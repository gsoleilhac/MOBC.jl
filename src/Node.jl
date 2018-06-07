mutable struct Node{T}
	m::JuMP.Model
	z::Float64
	x::Vector{Float64}
	zparent::Float64
	nbcover::Int
	f0::Set{Int}
	f1::Set{Int}
	bound1::Float64
	bound2::Float64
	cstr1::T
	cstr2::T
end

function Node(m, b1, b2, c1::T, c2::T) where T 
	Node(m, NaN, zeros(m.numCols), Inf, 0, Set{Int}(), Set{Int}(), Float64(b1), Float64(b2), c1, c2)
end
# Base.copy(n::Node) = begin warn("deepcopy of a node") ; return Node(copy(n.m), NaN, n.x, n.z, n.nbcover, copy(n.f0), copy(n.f1), n.bound1, n.bound2, n.cstr1, n.cstr2) end
shallowcopy(n::Node) = Node(n.m, NaN, n.x, n.z, n.nbcover, copy(n.f0), copy(n.f1), n.bound1, n.bound2, n.cstr1, n.cstr2)


mutable struct NodeParragh{T}
	m::JuMP.Model
	x::Vector{Vector{Float64}}
	nbcover::Int
	f0::Set{Int}
	f1::Set{Int}
	bound1::Float64
	bound2::Float64
	cstr1::T
	cstr2::T
end

function NodeParragh(m, b1, b2, c1::T, c2::T) where T
	NodeParragh(m, Vector{Float64}[], 0, Set{Int}(), Set{Int}(), Float64(b1), Float64(b2), c1, c2)
end
# Base.copy(n::NodeParragh) = NodeParragh(copy(n.m), Vector{Float64}[], n.nbcover, copy(n.f0), copy(n.f1))
shallowcopy(n::NodeParragh) = NodeParragh(n.m, Vector{Float64}[], n.nbcover, copy(n.f0), copy(n.f1), n.bound1, n.bound2, n.cstr1, n.cstr2)

Base.show(io::IO, n::NodeParragh) = print(io, "Node[f0=$(n.f0), f1=$(n.f1), <$(n.bound1)-$(n.bound2)>]")