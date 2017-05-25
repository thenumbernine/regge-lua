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

local Vertex = class()

function Vertex:init(args)
	for k,v in pairs(args) do self[k] = v end 
	self.edges = table()
	self.tris = table()
end

function Vertex:totalAngle()
	-- tvs is all triangles that hold this vertex
	local theta = 0
	for _,tv in ipairs(self.tris) do
		theta = theta + tv:angleForVertex(self)	
	end
	return theta
end

function Vertex:curvature()
	return 2 * math.pi - self:totalAngle()
end

function Vertex.__lt(a,b) return tostring(a) < tostring(b) end

function Vertex.__eq(a,b)
	return a.pos == b.pos 
end

local Edge = class()

function Edge:init(...)
	self.vtxs = table{...}
	self.tris = table()
end

function Edge:getVtxsByVtx(v)
	local a,b = self.vtxs:unpack()
	if v == a then return a,b end
	if v == b then return b,a end
end

local Triangle = class()

function Triangle:init(...)
	self.vtxs = table{...}
	self.edges = table()
end

function Triangle:angleForVertex(v)
	local a,b,c = self.vtxs:unpack()
	local ab = (b.pos - a.pos):normalize()
	local ac = (c.pos - a.pos):normalize()
	local bc = (c.pos - b.pos):normalize()
	if v == a then return math.acos(ab:dot(ac)) end
	if v == b then return math.acos(-ab:dot(bc)) end
	if v == c then return math.acos(ac:dot(bc)) end
end

function Triangle:getVtxsByVtx(v,w)
	local a,b,c = self.vtxs:unpack()
	if v and w then
		if v == a and w == b then return a, b, c end
		if v == a and w == c then return a, c, b end
		if v == b and w == a then return b, a, c end
		if v == b and w == c then return b, c, a end
		if v == c and w == a then return c, a, b end
		if v == c and w == b then return c, b, a end
	end
	if v then
		if v == a then return a, b, c end
		if v == b then return b, c, a end
		if v == c then return c, a, b end
	end
end

-- returns a,b,c where a,b are possessed by the edge
function Triangle:getVtxsByEdge(e)
	assert(#self.edges == 3)
	for _,te in ipairs(self.edges) do
		if te == e then
			local a,b = e.vtxs:unpack()
			local _, c = self.vtxs:find(nil, function(v) return v ~= a and v ~= b end)
			assert(a and b and c)
			return a,b,c
		end
	end
end

local function makevtx(x,y,z,t)
	for _,v in ipairs(vtxs) do
		if v.pos[1] == x 
		and v.pos[2] == y 
		and v.pos[3] == z 
		and v.pos[4] == t
		then 
			return v
		end
	end
	local v = Vertex{pos=vec4(x,y,z,t), index=#vtxs+1}
	vtxs:insert(v)
	return v
end

local function makeedge(a,b)
	for i,e in ipairs(edges) do
		if (e.vtxs[1] == a and e.vtxs[2] == b)
		or (e.vtxs[2] == a and e.vtxs[1] == b)
		then return e end
	end
	local e = Edge(a,b)
	edges:insert(e)
	a.edges:insert(e)
	b.edges:insert(e)
	return e
end

local function maketri(a,b,c, args)
	local vtxs = {a,b,c}
	table.sort(vtxs)
	for i,t in ipairs(tris) do
		if t.vtxs[1] == vtxs[1] 
		and t.vtxs[2] == vtxs[2] 
		and t.vtxs[3] == vtxs[3] 
		then
			return t
		end
	end
	
	local t = Triangle(table.unpack(vtxs))
	if args then
		for k,v in pairs(args) do
			t[k] = v
		end
	end
	tris:insert(t)
	
	t.edges:insert(makeedge(a,b))
	t.edges:insert(makeedge(a,c))
	t.edges:insert(makeedge(b,c))
	for _,e in ipairs(t.edges) do
		e.tris:insert(t)
	end
	for _,v in ipairs(t.vtxs) do
		v.tris:insert(t)
	end

	return t
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
		
	local slicesize = 10
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
		local ta = table()
		local tb = table()
		for i=1,#slice do
			local a,b = slice[i], slice[i%slicesize+1]
			local x,y = slice2[i], slice2[i%slicesize+1]
			if a and x and y then
				tb:insert(maketri(y,x,a, {color={1,0,0}}))
			end
			if a and b and y then
				ta:insert(maketri(a,b,y, {color={0,1,1}}))
			end
		end
		return ta, tb
	end

--[==[ grow randomly
	local sa = makeslice(0)
	local sb = makeslice(1)
	fuseslices(sa,sb)
	
	for i=1,20 do
		-- 1) pick a random external edge (with only one triangle)
		local openedges = edges:filter(function(edge) return #edge.tris == 1 end)
		-- 2) put a vertex opposite of it
		assert(#openedges > 0)
		local edge = openedges[math.random(#openedges)]
		local tri = edge.tris[1]
		local a,b,c = tri:getVtxsByEdge(edge)
		local v = makevtx((a.pos+b.pos-c.pos):unpack())
		maketri(b,v,c)
	end
--]==]
-- [==[ grow a new slice from the old, then minimize discrete curvature
	slices:insert(makeslice(1))
	slices:insert(makeslice(2))
	fuseslices(slices[1], slices[2])
	
	for i=3,3 do
		local pslice = slices[#slices]
		local slice = makeslice(i)
		slices:insert(slice)
		
		fuseslices(pslice, slice)
			
-- [=[ grad descent
		-- we want to perturb the latest slice
		-- such that the curvature of the previous slice
		-- is solved for
		-- phi = (R - R')^2
		-- dphi/dv = 2 (R - R') dR/dv
		-- R = 2 pi - sum angles
		-- dR/dv = -sum d/dv angle
		-- angle = acos(theta) 
		-- d/dv angle = -1 / sqrt(1 - theta^2) d/dv theta 
		-- theta = (v1-v0) dot (v2-v0) / (|v1-v0||v2-v0|)
		-- theta = (v1_i-v0_i) dot (v2_i-v0_i) / sqrt(sum_j (v1_j-v0_j)^2 sum_k (v2_k-v0_k)^2)
		-- d/dv1 theta = (v2 - v0) / (|v1-v0||v2-v0|)
		--		- ((v1-v0) dot (v2-v0)) / (|v1-v0|*|v2-v0|)^2
		--			 |v2-v0| * (v1-v0)/|v1-v0|
		-- = (|v1-v0| unit(v2-v0) - unit(v1-v0) dot unit(v2-v0)) / |v1-v0|^2
		-- d/dx x/|x| = I/|x| + x * (-1/|x|^2) * x/|x|
		--		= I / |x| - x*x / |x|^3)
		--		= (I (x.x) - x*x) / |x|^3
		-- d/dx x/|x| dot y/|y|
		--		= (y (x.x) - x (x.y)) / (|x|^3 |y|)
		--		= uy / lx - ux (ux . uy) / lx 
		local dphi_dvs = table()
		for j=1,slicesize do
			dphi_dvs[j] = vec4()
		end
		for k,pv in ipairs(pslice) do 		-- for each vertex in the previous slice 
			local R = pv:curvature()		-- look at the curvature at the vertex
			
			local thetacoord = math.atan2(pv.pos[2], pv.pos[1])
			local desR = math.abs(thetacoord) < (2 * math.pi / slicesize) and 3 or 0	-- look at the desired curvature

print('pslice',k,'R',R,'desR',desR)
			
			local a = pv
			for _,t in ipairs(pv.tris) do		-- for each triangle on that vertex ...
				for q,b in ipairs(t.vtxs) do	-- for each vertex on the triangle ...
					local j = slice:find(b)		-- if the vertex is in the next slice ...
					if j then
						local _,_,c = t:getVtxsByVtx(a,b)
						local vba = b.pos - a.pos
						local vca = c.pos - a.pos
						local lba = vba:length()
						local lca = vca:length()
						local uba = vba / lba
						local uca = vca / lca
						local theta = uba:dot(uca)
						local dtheta_dv = (uca - uba * theta) / lba
						local dangle_dv = dtheta_dv * (-1 / math.sqrt(1 - theta^2)) 
						local dR_dv = -dangle_dv
						local dphi_dv = dR_dv * (2 * (R - desR))
						b.pos = b.pos - dphi_dv
						--dphi_dvs[j] = dphi_dvs[j] + dphi_dv
					end
				end
			end
print('new R', pv:curvature())		
		end
		--[[
		for j=1,slicesize do
			local v = slice[j]
			local dt = 1
			v.pos = v.pos + dphi_dvs[j] * dt
		end
		--]]
--]=]			
	end
--]==]

	for i,slice in ipairs(slices) do
		print('slice',i)
		for j,v in ipairs(slice) do
			print('vtx',j,'curvature',v:curvature())
		end
	end

	-- recenter
	local com = vec4(0,0,0,0)
	for _,v in ipairs(vtxs) do
		com = com + v.pos
	end
	com = com / #vtxs
	for i=1,#vtxs do
		vtxs[i].pos = vtxs[i].pos - com
	end	

	--[[
	for i,v in ipairs(vtxs) do
		print(v.index, v:totalAngle())
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
		gl.glVertex3d(v.pos:unpack(1,3))
	end
	gl.glEnd()

	gl.glColor3f(1,1,0)
	gl.glBegin(gl.GL_LINES)
	for _,e in ipairs(edges) do
		for _,a in ipairs(e.vtxs) do
			gl.glVertex3d(a.pos:unpack(1,3))
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
			gl.glColor4d(1,0,0,alpha)
		end
		for _,a in ipairs(t.vtxs) do
			gl.glVertex3d(a.pos:unpack(1,3))
		end
	end
	gl.glEnd()
end

function App:updateGUI()

end

App():run()
