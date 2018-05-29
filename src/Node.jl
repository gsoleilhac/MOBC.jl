mutable struct Node
	m::JuMP.Model
	z::Float64
	x::Vector{Float64}
	zparent::Float64
	nbcover::Int
	f0::Set{Int}
	f1::Set{Int}
end
Node(m) = Node(copy(m), NaN, Float64[], Inf, 0, Set{Int}(), Set{Int}())
Base.copy(n::Node, copy_model=true) = Node(copy_model ? copy(n.m) : n.m, NaN, Float64[], n.z, n.nbcover, copy(n.f0), copy(n.f1))

Base.show(io::IO, n::Node) = begin
	print(io, "Node(")
	if isnan(n.z)
		print(io, "unsolved")
	else
		print(io, "z=$(n.z)")
		print(io, "x=$(n.x)")
	end
	print(io, ")")
end

mutable struct NodeParragh
	m::JuMP.Model
	x::Vector{Vector{Float64}}
	nbcover::Int
	f0::Set{Int}
	f1::Set{Int}
end
NodeParragh(m) = NodeParragh(copy(m), Vector{Float64}[], 0, Set{Int}(), Set{Int}())
Base.copy(n::NodeParragh, copy_model=true) = NodeParragh(copy_model ? copy(n.m) : n.m, Vector{Float64}[], n.nbcover-1, copy(n.f0), copy(n.f1))