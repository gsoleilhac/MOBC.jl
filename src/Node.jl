mutable struct Node
	m::JuMP.Model
	z::Float64
	x::Vector{Float64}
	zparent::Float64
	nbcover::Int
end
Node(m) = Node(copy(m), NaN, [], Inf, 0)
Base.copy(n::Node) = Node(copy(n.m), NaN, [], n.z, n.nbcover-1)

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