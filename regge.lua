#!/usr/bin/env luajit
require 'ext'
require 'vec'
local ffi = require 'ffi'
local ImGuiApp = require 'imguiapp'
local gl = require 'gl'
local ig = require 'ffi.imgui'
local bit = bit32 or require 'bit'
local vec4f = require 'ffi.vec.vec4f'

local View = require 'glapp.view'
local Orbit = require 'glapp.orbit'

require 'graph'

local App = class(Orbit(View.apply(ImGuiApp)))

App.title = 'Regge calculus demo'
App.viewDist = 2


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

--[=[ lattice
	local gridsize = vec3(3,3,3)
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
	
	for j,v in ipairs(vtxs) do
		print('vtx',j,'curvature',v:curvature())
	end
--]=]
-- [=[
	--[[
	deficit angle
	in 2D is 2 pi - sum of triangle angles meeting at a vertex
	in 3D is 2 pi - sum of tetrahedron angles meeting at an edge
	in 4D is 2 pi - sum of 4-simplex angles meeting at a triangle
	--]]

	local slices = table()
	
	local baselength = 1
	local slicesize = 20

	local function vtxpos(i,z,dr)
		dr = dr or 0
		local radius = baselength * slicesize / (2 * math.pi)
		local r = radius + dr
		local th = 2 * math.pi * (i-1 - .5*z) / slicesize
		local c, s = r * math.cos(th), r * math.sin(th)
		return c, s, z * math.sqrt(3) / 2 * baselength, 0
	end

	local function makeslice(z)
		local slice = table()
		for i=1,slicesize do
			local v = makevtx(vtxpos(i,z))
			local theta = math.atan2(v.pos[2], v.pos[1])
			v.R = -.3	 --math.abs(theta) < (2 * math.pi / slicesize) and .1 or -.1
			slice:insert(v)
		end
		return slice
	end

	local function extrudeslice(pslice)
		-- [[
		local slice = table()
		for i=1,#pslice do
			local a = pslice[i]
			local b = pslice[i%#pslice+1]
			local ts = a:getTrisByVtx(b)
			if #ts ~= 1 then
				error("expected 1 triangle but found "..#ts)
			end
			local t = ts[1]
			local _,_,c = t:getVtxsByVtx(a,b)
			--local vpos = (a.pos + b.pos - c.pos * 2):normalize()  + c.pos
			local vpos = ((a.pos + b.pos) * .5 - c.pos) + (a.pos + b.pos) * .5
			local v = makevtx(vpos:unpack())
			v.R = c.R or .5 * (a.R + b.R)
			slice:insert(v)
		end
		slice:insert(1, slice:remove())
		--]]
		--[[
		for i=#slice,1,-1 do
			local a = slice[i]
			local b = slice[i%#slice+1]
			if (a.pos-b.pos):length() > 1.5 * baselength then
				local v = makevtx(((a.pos + b.pos) * .5):unpack())
				v.R = (a.R + b.R) * .5
				slice:insert(i+1,v)
			end
		end
		--]]
		--slice:insert(1, slice:remove())
		return slice
	end

	-- TODO match nearest?
	local function fuseslices(slice, slice2)
		local ta = table()
		local tb = table()
		for i=1,#slice do
			local a,b = slice[i], slice[i%#slice+1]
			local x,y = slice2[i], slice2[i%#slice2+1]
			if a and x and y then
				tb:insert(maketri(y,x,a, {color={1,0,0}}))
			end
			if a and b and y then
				ta:insert(maketri(a,b,y, {color={0,1,1}}))
			end
		end
		return ta, tb
	end

	local function convergeSlice(slice, pslice)
		-- we want to perturb the latest slice
		-- such that the curvature of the previous slice
		-- is solved for
		-- phi = (R - R')^2
		-- dphi/dv = 2 (R - R') dR/dv
		-- R = 2 pi - sum angles
		-- dR/dv = -sum d/dv angle
		-- angle = acos(theta) 
		-- d/dv angle = -1 / sqrt(1 - theta^2) d/dv theta 
		-- theta = (b-a) dot (c-a) / (|b-a||c-a|)
		-- d/dx x/|x| = I/|x| + x * (-1/|x|^2) * x/|x|
		--		= I / |x| - x*x / |x|^3)
		--		= (I (x.x) - x*x) / |x|^3
		-- d/dx x/|x| dot y/|y|
		--		= (y (x.x) - x (x.y)) / (|x|^3 |y|)
		--		= (uy - ux (ux . uy)) / |x|
		local totalPhi = 0
		for _,v in ipairs(pslice) do
			local dK = v:curvature() - v.R
			totalPhi = totalPhi + dK * dK
		end
		print('total phi', totalPhi)		
		
		local dt = 1
		for try=1,100 do	
			local dphi_dvs = table()
			for j=1,#slice do
				dphi_dvs[j] = vec4()
			end
			for k,pv in ipairs(pslice) do 		-- for each vertex in the previous slice 
				local R = pv:curvature()		-- look at the curvature at the vertex
				
--io.write('pslice\t',k,'\tR\t',R,'\tdesR\t',pv.R,'\t|R-desR|\t', math.abs(R - pv.R))
				
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
							local dangle_dv = dtheta_dv * (-1 / math.sqrt(1 - theta * theta))
							local dR_dv = -dangle_dv
							local dphi_dv = dR_dv * (2 * (R - pv.R))
							dphi_dvs[j] = dphi_dvs[j] + dphi_dv
						end
					end
				end
--io.write('\tnew R\t', pv:curvature(),'\tnew |R-desR|\t',math.abs(R-pv.R),'\n')
			end
			for j=1,#slice do
				local v = slice[j]
				v.pos = v.pos - dphi_dvs[j] * dt
			end
		
			local totalPhi = 0
			for _,v in ipairs(pslice) do
				local dK = v:curvature() - v.R
				totalPhi = totalPhi + dK * dK
			end
			print('total phi', totalPhi)		
		
		end
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
	for i=1,11 do
		local pslice = slices[i-1]
		local slice
		if i > 2 then
			print('extruding slice '..i)
			slice = extrudeslice(slices[i-1])
		else
			slice = makeslice(i)
		end
		slices:insert(slice)
	
		if i >= 2 then
			local ta, tb = fuseslices(pslice, slice)
			if i >= 3 then
				convergeSlice(slice, pslice)
			end
			-- [[
			-- now split edges too big
			-- they're all about 1 for now
			for i=#slice,1,-1 do
				local va = slice[i]
				local vb = slice[(i-2)%#slice+1]
				local ts = va:getTrisByVtx(vb)
				if #ts ~= 1 then
					error("FOUND #ts ~= 1: ".. #ts)
					va.color = {1,0,0}
					vb.color = {1,0,0}
				else
					assert(#ts == 1)
					local t = ts[1]
					local edge = t:getEdgeForVtxs(va,vb)
					local len = (vb.pos - va.pos):length()
					
					if len > 1.2 then 
						
						-- notice this is destructive
						-- it adds a new vtx midway through edge
						-- removes edge
						-- removes all tris that reference edge
						-- and replaces each with two tris that instead use the new vtx
						local vn, ea, eb, newts = edge:split() 
						assert(#newts == 2)
						newts[1].color = {0,1,0}
						newts[2].color = {1,0,1}
						vn.R = (va.R + vb.R) * .5
						slice:insert(i, vn)
					end
				end
			end
			--]]
		end
	end
--]==]

	for i,slice in ipairs(slices) do
		print('slice',i)
		for j,v in ipairs(slice) do
			print('vtx',j,'curvature',v:curvature())
		end
	end
--]=]

	-- recenter
	local com = vec4(0,0,0,0)
	for _,v in ipairs(vtxs) do
		com = com + v.pos
	end
	com = com / #vtxs
	for i=1,#vtxs do
		vtxs[i].pos = vtxs[i].pos - com
	end	

	-- calculate triangle normals
	for _,t in ipairs(tris) do
		local a,b,c = t.vtxs:unpack()
		t.normal = vec3.cross(c.pos - b.pos, b.pos - a.pos):normalize()
	end

	-- calculate vertex normals
	for _,v in ipairs(vtxs) do v.normal = vec3() end
	for _,t in ipairs(tris) do
		for _,v in ipairs(t.vtxs) do
			if t.normal:dot(v.normal) > 0 then
				v.normal = v.normal + t.normal
			else
				v.normal = v.normal - t.normal
			end
		end
	end
	for _,v in ipairs(vtxs) do v.normal = v.normal:normalize() end

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
	App.super.update(self)
	gl.glClear(bit.bor(gl.GL_COLOR_BUFFER_BIT, gl.GL_DEPTH_BUFFER_BIT))

gl.glEnable(gl.GL_LIGHT0)
gl.glLightModelf(gl.GL_LIGHT_MODEL_LOCAL_VIEWER, gl.GL_FALSE)
gl.glLightModelf(gl.GL_LIGHT_MODEL_TWO_SIDE, gl.GL_FALSE)
gl.glLightModelfv(gl.GL_LIGHT_MODEL_AMBIENT, vec4f(0,0,0,0):ptr())
gl.glLightfv(gl.GL_LIGHT0, gl.GL_POSITION, vec4f(0,0,0,1):ptr())
gl.glLightfv(gl.GL_LIGHT0, gl.GL_AMBIENT, vec4f(.3,.3,.3,1):ptr())
gl.glLightfv(gl.GL_LIGHT0, gl.GL_DIFFUSE, vec4f(1,1,1,1):ptr())
gl.glColorMaterial(gl.GL_FRONT_AND_BACK, gl.GL_DIFFUSE)
gl.glEnable(gl.GL_COLOR_MATERIAL)

	if drawVertexes[0] then
		gl.glColor3f(1,1,1)
		gl.glPointSize(3)
		gl.glBegin(gl.GL_POINTS)
		for _,v in ipairs(vtxs) do
			if v.color then
				gl.glColor3f(table.unpack(v.color))
			else
				gl.glColor3f(1,1,1)
			end
			gl.glVertex3d(v.pos:unpack(1,3))
		end
		gl.glEnd()
	end

	if drawEdges[0] then
		gl.glColor3f(1,1,0)
		gl.glBegin(gl.GL_LINES)
		for _,e in ipairs(edges) do
			for _,a in ipairs(e.vtxs) do
				gl.glVertex3d(a.pos:unpack(1,3))
			end
		end
		gl.glEnd()
	end

	if drawNormals[0] then
		gl.glColor3f(1,0,1)
		gl.glBegin(gl.GL_LINES)
		for _,t in ipairs(tris) do
			local a,b,c = t.vtxs:unpack()
			local x = vec3(((a.pos + b.pos + c.pos) * (1/3)):unpack())
			gl.glVertex3d(x:unpack())
			gl.glVertex3d((x + t.normal):unpack())
		end
		gl.glEnd()
	end

	if useBlend[0] then
		gl.glBlendFunc(gl.GL_SRC_ALPHA, gl.GL_ONE)
		gl.glEnable(gl.GL_BLEND)
		gl.glDisable(gl.GL_DEPTH_TEST)
	end
	local alpha = 1
	if useLighting[0] then 
		gl.glEnable(gl.GL_LIGHTING)
	end
	if drawTriangles[0] then
		gl.glBegin(gl.GL_TRIANGLES)
		for _,t in ipairs(tris) do
			if t.color then
				local c = t.color
				gl.glColor4d(c[1], c[2], c[3], alpha)
			else
				gl.glColor4d(1,0,0,alpha)
			end
			for _,v in ipairs(t.vtxs) do
				gl.glNormal3d(v.normal:unpack(1,3))
				gl.glVertex3d(v.pos:unpack(1,3))
			end
		end
		gl.glEnd()
	end

	gl.glDisable(gl.GL_LIGHTING)
	gl.glDisable(gl.GL_BLEND)
	gl.glEnable(gl.GL_DEPTH_TEST)
end

drawVertexes = ffi.new('bool[1]', true)
drawEdges = ffi.new('bool[1]', true)
drawTriangles = ffi.new('bool[1]', true)
drawNormals = ffi.new('bool[1]', true)
useLighting = ffi.new('bool[1]', true)
useBlend = ffi.new('bool[1]', true)
function App:updateGUI()
	ig.igCheckbox('vertexes', drawVertexes)
	ig.igCheckbox('edges', drawEdges)
	ig.igCheckbox('triangles', drawTriangles)
	ig.igCheckbox('normals', drawNormals)
	ig.igCheckbox('light', useLighting)
	ig.igCheckbox('transparent', useBlend)
end

App():run()
