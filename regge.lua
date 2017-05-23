#!/usr/bin/env luajit
require 'ext'
local ImGuiApp = require 'imguiapp'
local gl = require 'gl'
local ig = require 'ffi.imgui'
local bit = bit32 or require 'bit'
local vec3 = require 'vec.vec3'

local View = require 'glapp.view'
local Orbit = require 'glapp.orbit'

local App = class(Orbit(View.apply(ImGuiApp)))

App.title = 'Regge calculus demo'
App.viewDist = 2

local vtxs = table()
local edges = table()
local faces = table()

local function makevtx(x,y,z)
	for i,v in ipairs(vtxs) do
		if v[1] == x and v[2] == y and v[3] == z then return i end
	end
	vtxs:insert(vec3(x,y,z))
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

function App:init()
	App.super.init(self)
	for j=1,10 do
		for i=1,10 do
			local a = makevtx(i,j,0)
			local b = makevtx(i+1,j,0)
			makeedge(a,b)
			local c = makevtx(i,j+1,0)
			makeedge(a,c)
		end
	end
end

function App:update()
	gl.glClear(bit.bor(gl.GL_COLOR_BUFFER_BIT, gl.GL_DEPTH_BUFFER_BIT))
	App.super.update(self)

	gl.glColor3f(1,1,1)
	gl.glPointSize(2)
	gl.glBegin(gl.GL_POINTS)
	for _,v in ipairs(vtxs) do
		gl.glVertex3d(v:unpack())
	end
	gl.glEnd()

	gl.glColor3f(1,1,0)
	gl.glBegin(gl.GL_LINES)
	for _,e in ipairs(edges) do
		for _,a in ipairs(e) do
			gl.glVertex3d(vtxs[a]:unpack())
		end
	end
	gl.glEnd()
end

function App:updateGUI()
	ig.igText('green = F')
	ig.igText('red = curl F')
	ig.igText('blue = curl^-1 curl F')

	if ig.igButton'apply' then
		apply(F2)
	end
end


App():run()
