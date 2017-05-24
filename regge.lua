#!/usr/bin/env luajit
require 'ext'
require 'vec'
local ImGuiApp = require 'imguiapp'
local gl = require 'gl'
local ig = require 'ffi.imgui'
local bit = bit32 or require 'bit'

local View = require 'glapp.view'
local Orbit = require 'glapp.orbit'

local App = class(Orbit(View.apply(ImGuiApp)))

App.title = 'Regge calculus demo'
App.viewDist = 2

local vtxs = table()
local edges = table()
local tris = table()
local tets = table()

local Triangle = class()

function Triangle:init(a,b,c)
	self[1] = a
	self[2] = b
	self[3] = c
end

function Triangle:angleForVertex(v)
	local a,b,c = table.unpack(self)
	local va,vb,vc = vtxs[a], vtxs[b], vtxs[c]
	local ab = (vb - va):normalize()
	local ac = (vc - va):normalize()
	local bc = (vc - vb):normalize()
	if v == a then return math.acos(ab:dot(ac)) end
	if v == b then return math.acos(-ab:dot(bc)) end
	if v == c then return math.acos(ac:dot(bc)) end
end

function Triangle:getvtxs(i)
	local a,b,c = table.unpack(self)
	if i == a then return a, b, c end
	if i == b then return b, c, a end
	if i == c then return c, a, b end
end

local Vtx = class(vec4)

function Vtx:init(t, i)
	for k,v in pairs(t) do self[k] = v end
	self.i = i
end

function Vtx:neighborTris()
	local i = self.i
	local tvs = table()
	for _,t in ipairs(tris) do
		if table.find(t, i) then
			tvs:insert(t)
		end
	end
	return tvs
end

function Vtx:totalAngle()
	local i = self.i
	-- tvs is all triangles that hold vertex i
	local tvs = self:neighborTris()
	local theta = 0
	for _,tv in ipairs(tvs) do
		theta = theta + tv:angleForVertex(i)	
	end
	return theta
end

function Vtx:curvature()
	return 2 * math.pi - self:totalAngle()
end

local function makevtx(x,y,z,t)
	for i,v in ipairs(vtxs) do
		if v[1] == x 
		and v[2] == y 
		and v[3] == z 
		and v[4] == t
		then 
			return i 
		end
	end
	vtxs:insert(Vtx({x,y,z,t}, #vtxs+1))
	return #vtxs
end

local function makeedge(a,b)
	for i,e in ipairs(edges) do
		if (e[1] == a and e[2] == b)
		or (e[2] == a and e[1] == b)
		then return i end
	end
	edges:insert{a,b}
	return #edges
end

local function maketri(a,b,c, args)
	local tri = Triangle(a,b,c)
	table.sort(tri)
	if args then
		for k,v in pairs(args) do
			tri[k] = v
		end
	end
	for i,t in ipairs(tris) do
		if t[1] == tri[1] and t[2] == tri[2] and t[3] == tri[3] then
			return i
		end
	end
	tris:insert(tri)
	
	makeedge(a,b)
	makeedge(a,c)
	makeedge(b,c)
	
	return #tris
end

local function maketet(a,b,c,d)
	local as = {a,b,c,d}
	table.sort(as)
	for i,t in ipairs(tets) do
		if t[1] == as[1] and t[2] == as[2] and t[3] == as[3] and t[4] == as[4] then
			return i
		end
	end
	tets:insert(as)
	
	maketri(a,b,c)
	maketri(a,b,d)
	maketri(a,c,d)
	maketri(b,c,d)
	
	return #tets
end

-- height: 1/2^2 + y^2 = 1 <=> y = sqrt(3/4) = sqrt(3) / 2
-- 2D center is x=1/2, y=x * tan(30) = 1/2 * 1/sqrt(3) = 1/(2*sqrt(3))
--  distance from any base point to base center: 1/sqrt(3)
-- depth: 1/3^2 + z^2 = 1 <=> z = sqrt(8/9) = 2 sqrt(2) / 3
-- inclination from base to top: atan(z/base) = atan(2 sqrt(2) / 3 / (1/sqrt(3))) = atan(2 sqrt(2/3))
local function vtxcoord(x,y,z)
	--[[
	return 
		x - .5 * y ,
		y * .5 * math.sqrt(3) ,
		z * math.sqrt(2) * 2/3,
		0
	--]]
	return x, y, z, 0
end

function App:init()
	App.super.init(self)
	local gridsize = vec3(10,10,10)

	--[[
	for k=1,gridsize[3] do
		for j=1,gridsize[2] do
			for i=1,gridsize[1] do
				local function makequad(i,j,k)
					local vi = makevtx(vtxcoord(i, j, k))
					local vx, vy, vxy
					if i < gridsize[1] then
						vx = makevtx(vtxcoord(i+1, j, k))
						makeedge(vi,vx)
						if j < gridsize[2] then
							vxy = makevtx(vtxcoord(i+1, j+1, k))
							makeedge(vi,vxy)
							maketri(vi,vx,vxy)
						end
					end
					if j < gridsize[2] then
						vy = makevtx(vtxcoord(i, j+1, k))
						makeedge(vi,vy)
						if i < gridsize[1] then
							maketri(vxy,vy,vi)
						end
					end
					return vi, vx, vy, vxy
				end

				local function maketetandtris(a,b,c,d)
					maketri(a,b,c)
					maketri(a,b,d)
					maketri(a,c,d)
					maketri(b,c,d)
					makeedge(a,b)
					makeedge(a,c)
					makeedge(a,d)
					makeedge(b,c)
					makeedge(b,d)
					makeedge(c,d)
					maketet(a,b,c,d)
				end

				local function makeprism(a,b,c, d,e,f)
					maketetandtris(d,e,f,c)
					maketetandtris(d,e,b,c)
					maketetandtris(d,a,b,c)
				end

				local vi, vx, vy, vxy = makequad(i,j,k)
				if k < gridsize[3] then
					local vz, vxz, vyz, vxyz = makequad(i,j,k+1)
					if vi and vx and vy and vxy
					and vz and vxz and vyz and vxyz
					then
						makeprism(vi, vx, vxy, vz, vxz, vxyz)
						makeprism(vxy, vy, vi, vxyz, vyz, vz)
					end
				end
			end
		end
	end
	--]]

	--[[
	deficit angle
	in 2D is 2 pi - sum of triangle angles meeting at a vertex
	in 3D is 2 pi - sum of tetrahedron angles meeting at an edge
	in 4D is 2 pi - sum of 4-simplex angles meeting at a triangle
	--]]

	local slices = table()
		
	local slicesize = 50
	local radius = slicesize / (2 * math.pi)

	local function vtxpos(i,z,dr)
		dr = dr or 0
		local r = radius + dr
		local th = 2 * math.pi * (i-1 - .5*z) / slicesize
		local c, s = r * math.cos(th), r * math.sin(th)
		return c, s, z * math.sqrt(3) / 2, 0
	end

	local function makeslice(z)
		local slice = table()
		for i=1,slicesize do
			slice:insert(makevtx(vtxpos(i,z)))
		end
		return slice
	end

	local function fuseslices(slice, slice2)
		for i=1,#slice do
			local a,b = slice[i], slice[i%slicesize+1]
			local x,y = slice2[i], slice2[i%slicesize+1]
			if a and x and y then
				maketri(y,x,a, {color={1,0,0}})
			end
			if a and b and y then
				maketri(a,b,y, {color={0,1,1}})
			end
		end
	end

	local function findedge(a,b)
		for _,e in ipairs(edges) do
			if (e[1] == a and e[2] == b)
			or (e[1] == b and e[2] == a)
			then
				return true
			end
		end
	end

	for i=1,20 do
		slices:insert(makeslice(i))
		if i > 1 then
			local pslice = slices[i-1]
			local slice = slices[i]
			fuseslices(pslice, slice)
			
			if i > 2 then
				-- now perturb radially to adjust to curvature
				-- phi = (R - R')^2
				-- dphi/dr = 2 (R - R') dR/dr
				-- R = 2 pi - sum angles
				-- dR/dr = -sum d/dr angle
				-- angle = acos(theta) 
				-- dangle / dr = -1 / sqrt(1 - theta^2) d/dr theta 
				-- theta = (v1-v0) dot (v2-v0) / (|v1-v0||v2-v0|)
				-- if v is v1 then dtheta/dr = dv1/dr dot (v2-v0) / (|v1-v0||v2-v0|) - (v1-v0) dot (v2-v0) dv1/dr / (|v1-v0||v2-v0|)^2
				-- if v is v2 then dtheta/dr = dv2/dr dot (v2-v0) / (|v1-v0||v2-v0|) - (v1-v0) dot (v2-v0) dv2/dr / (|v1-v0||v2-v0|)^2
				-- dv1/dr = v1/|v1|
				-- we want to perturb the latest slice
				-- such that the curvature of the previous slice
				-- is solved for
				local dphi_drs = table()
				for j=1,slicesize do
					dphi_drs[j] = 0
				end
				for j=1,slicesize do	-- for each of the latest slice vtxs
					local vi = slice[j]
					local v = vtxs[vi]
					for k=1,slicesize do	-- for each prev slice vtx that influences it ...
						local pvi = pslice[k]
						local pv = vtxs[pvi]
						if findedge(vi, pvi) then
							local tvs = pv:neighborTris()
							local R = pv:curvature()
						
							-- TODO make this a function of 'r'
							local theta = math.atan2(pv[2], pv[1])
							local desR = math.abs(theta) < (2 * math.pi / slicesize) and 3 or 0
							
							local dR_dr = 0
							for _,t in ipairs(tvs) do
								local a,b,c = t:getvtxs(vi)
								if a then
									local va, vb, vc = vtxs[a], vtxs[b], vtxs[c]
									local vba = vb - va
									local vca = vc - va
									local theta = vba:dot(vca) / (vba:length() * vca:length())
dtheta_dr = (b == pvi or c == pvi) and .1 or 0
									local dangle_dr = -1 / math.sqrt(1 - theta^2) * dtheta_dr 
									dR_dr = dR_dr - dangle_dr
								end
							end
							local dphi_dr = 2 * (R - desR) * dR_dr
							dphi_drs[j] = dphi_drs[j] + dphi_dr
						end
					end
				end
				for j=1,slicesize do
					local v = vtxs[slice[j]]
					-- r' = r + dr
					-- r'/r = 1 + dr/r
					local r = math.sqrt(v[1]^2 + v[2]^2)
					for l=1,2 do
						v[l] = v[l] * (1 + dphi_drs[j] / r)
					end
				end
			end
		end
	end

	-- recenter
	local com = vec4(0,0,0,0)
	for _,v in ipairs(vtxs) do
		com = com + v
	end
	com = com / #vtxs
	for i=1,#vtxs do
		for j=1,4 do
			vtxs[i][j] = vtxs[i][j] - com[j]
		end
	end	

	--[[
	for i,v in ipairs(vtxs) do
		print(v.i, v:totalAngle())
	end
	--]]
end

function App:initGL()
	App.super.initGL(self)
	gl.glEnable(gl.GL_BLEND)
	gl.glBlendFunc(gl.GL_SRC_ALPHA, gl.GL_ONE_MINUS_SRC_ALPHA)
	gl.glDisable(gl.GL_CULL_FACE)
	gl.glEnable(gl.GL_DEPTH_TEST)
end

function App:update()
	gl.glClear(bit.bor(gl.GL_COLOR_BUFFER_BIT, gl.GL_DEPTH_BUFFER_BIT))
	App.super.update(self)

	gl.glColor3f(1,1,1)
	gl.glPointSize(2)
	gl.glBegin(gl.GL_POINTS)
	for _,v in ipairs(vtxs) do
		gl.glVertex3d(v:unpack(1,3))
	end
	gl.glEnd()

	gl.glColor3f(1,1,0)
	gl.glBegin(gl.GL_LINES)
	for _,e in ipairs(edges) do
		for _,a in ipairs(e) do
			gl.glVertex3d(vtxs[a]:unpack(1,3))
		end
	end
	gl.glEnd()

	local alpha = 1
	gl.glBegin(gl.GL_TRIANGLES)
	for _,t in ipairs(tris) do
		if t.color then
			local c = t.color
			gl.glColor4d(c[1], c[2], c[3], alpha)
		else
			gl.glColor4d(1,1,0,alph)
		end
		for _,a in ipairs(t) do
			gl.glVertex3d(vtxs[a]:unpack(1,3))
		end
	end
	gl.glEnd()
end

function App:updateGUI()

end

App():run()
