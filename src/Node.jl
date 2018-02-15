mutable struct Node
	m::JuMP.Model
	z::Float64
	x::Vector{Float64}
	f0::Vector{Int}
	f1::Vector{Int}
	seen::Vector{Vector{Any}}
	ubz1::Float64
	ubz2::Float64
	Node(m,f0,f1,seen) = new(m,NaN,[],f0,f1,seen,Inf,Inf)
	Node(m, z, x, f0, f1, seen) = new(m,z,x,f0,f1,seen,Inf,Inf)
	Node(m, z, x, f0, f1, seen, ubz1, ubz2) = new(m,z,x,f0,f1,seen,ubz1, ubz2)
end


Base.copy(n::Node) = Node(copy(n.m), NaN, [], copy(n.f0), copy(n.f1), copy(n.seen), n.ubz1, n.ubz2)

Base.show(io::IO, n::Node) = begin
	if isnan(n.z)
		println(io, "\tunsolved")
	else
		println(io, "\t z=$(n.z)")
		println(io, "\t x=$(n.x)")
	end
	println(io, "\tf0=$(n.f0)")
	println(io, "\tf1=$(n.f1)")
	print(io, ")")
end