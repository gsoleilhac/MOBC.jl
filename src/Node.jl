mutable struct Node
	m::JuMP.Model
	z::Float64
	x::Vector{Float64}
	zparent::Float64
end
Node(m) = Node(m, NaN, [], Inf)
Base.copy(n::Node) = Node(copy(n.m), NaN, [], n.z)

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