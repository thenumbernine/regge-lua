-- TODO OOP this

vtxs = table()
edges = table()
tris = table()
tets = table()

Vertex = class()

function Vertex:init(args)
	for k,v in pairs(args) do self[k] = v end 
	self.edges = table()
	self.tris = table()
end

--[[
2D:
for a point,
for each triangle touching that point
look at the edges opposite that point
calculate the deficit angle 
cos(theta) = (c^2 - a^2 - b^2 / (2ab)
where 'c' is the opposite edge length
and 'a' and 'b' are lengths of the other triangle edges

3D:
for an edge
for each tetrad touching that edge
cos theta_ab = (cos theta_acd - cos theta_abc cos theta_abd) / (sin theta_abc sin theta_abd)
TODO implement this under Edge?
--]]
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

function Vertex:getTrisByVtx(v)
	local ts = table()
	for _,t in ipairs(self.tris) do
		if t.vtxs:find(v) then
			ts:insert(t)
		end
	end
	return ts
end

function Vertex.__lt(a,b) return tostring(a) < tostring(b) end

function Vertex.__eq(a,b)
	return a.pos == b.pos 
end

Edge = class()

function Edge:init(...)
	self.vtxs = table{...}
	self.tris = table()
end

function Edge:getVtxsByVtx(v)
	local a,b = self.vtxs:unpack()
	if v == a then return a,b end
	if v == b then return b,a end
end

function Edge:getTrisByVtxs(a,b)
	return self.tris:filter(function(t)
		return t.edges:find(self)
	end)
end

function Edge:remove()
	-- remove this edge from the main list
	for _,t in ipairs(self.tris) do
		t.edges:removeObject(self)
	end
	for _,v in ipairs(self.vtxs) do
		v.edges:removeObject(self)
	end
	edges:removeObject(self)
end

function Edge:split()
	local a,b = self.vtxs:unpack()
	local n = makevtx(((a.pos+b.pos)*.5):unpack())

	-- remove this edge from its vertices
	for _,v in ipairs(self.vtxs) do
		v.edges:removeObject(self)
	end

	-- remove this edge from its triangles, splitting and replacing them as we go
	local e1, e2
	local ts = table(self.tris)
	local newts = table()
	for _,t in ipairs(ts) do
		local _,_,c = t:getVtxsByEdge(self)

		-- remove this from the master list
		t:remove()

		-- make two new tris
		local nt = maketri(n,b,c)
		if t.color then nt.color = vec3(t.color:unpack()) end
		newts:insert(nt)
		
		local nt = maketri(c,a,n)
		if t.color then nt.color = vec3(t.color:unpack()) end
		newts:insert(nt)

		-- should be the same for all new tris created
		e1 = t:getEdgeForVtxs(a, n)
		e2 = t:getEdgeForVtxs(n, b)
	end
	
	self:remove()
	
	return n, e1, e2, newts
end

Triangle = class()

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

-- provide 2 vertexes, returns the 3rd
function Triangle:get3rdVtx(v,w)
	local a,b,c = self.vtxs:unpack()
	if v == a and w == b then return c end
	if v == a and w == c then return b end
	if v == b and w == a then return c end
	if v == b and w == c then return a end
	if v == c and w == a then return b end
	if v == c and w == b then return a end
end

-- provide 1 vertex, returns the 2nd and 3rd vertexes
function Triangle:get23Vtxs(v)
	local a,b,c = self.vtxs:unpack()
	if v == a then return b, c end
	if v == b then return c, a end
	if v == c then return a, b end
end

function Triangle:getVtxsByVtx(v,w)
	local a,b,c = self.vtxs:unpack()
	if v and w then
		return v, w, self:get3rdVtx(v, x)
	end
	if v then
		return v, self:get23Vtxs(v)
	end
	return a,b,c
end

function Triangle:remove()
	tris:removeObject(self)
	for _,e in ipairs(self.edges) do
		e.tris:removeObject(self)
	end
	for _,v in ipairs(self.vtxs) do
		v.tris:removeObject(self)
	end
	-- TODO remove any edges / vertices that are alone
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

function Triangle:getEdgeForVtxs(a,b)
	for _,e in ipairs(self.edges) do
		if e.vtxs:find(a) and e.vtxs:find(b) then return e end
	end
end

function makevtx(x,y,z,t)
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

function makeedge(a,b)
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

function maketri(a,b,c, args)
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

function maketet(a,b,c,d)
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
