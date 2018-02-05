module MOBC
using JuMP, vOptGeneric

@enum NODE_STATUS SOLVED BOUNDED PENDING

mutable struct Node
	m::JuMP.Model
	status::NODE_STATUS
	isinteger::Bool
	z::Float64
	x::Vector{Float64}
	f0::Vector{Int}
	f1::Vector{Int}
	feasible::Bool
	seen::Vector{Vector{Int}}
	left::Node
	right::Node
	Node(m,z,x,f0,f1,seen)=new(m,SOLVED,z,x,f0,f1,true,seen)
	Node(m,f0,f1,seen)=new(m,PENDING,0,[],f0,f1,true,seen)
end

Base.show(io::IO, n::Node) = begin
	print(io, "Node(\t$(n.status), \n")
	n.status != PENDING && println(io, "\tz=$(n.z)")
	n.status != PENDING && println(io, "\tx=$(n.x)")
	println(io, "\tf0=$(n.f0)")
	println(io, "\tf1=$(n.f1)")
	n.status != PENDING && println(io, "\tfeasible:$(n.feasible)")
	isdefined(n, :left) && println(io, "\thas left child")
	isdefined(n, :right) && println(io, "\thas right child")
	print(io, ")")
end

isonlybinary(m, obj) = all(m.colCat[map(x->getfield(x, :col), obj.aff.vars)] .== :Bin)
hasinteger(m, obj) = any(m.colCat[map(x->getfield(x, :col), obj.aff.vars)] .== :Int)
isinteger(x, ind_INT) = all(x->x==1.||x==0., x[ind_INT])

function solve_BC(vm)
	vd = getvOptData(vm)
	@assert length(vd.objs) == 2
	@assert isonlybinary(vm, vd.objs[1])
	@assert !hasinteger(vm, vd.objs[1]) && !hasinteger(vm, vd.objs[2])

	if vd.objSenses[1] == :Max
		vd.objSenses[1] = :Min
		vd.objs[1] = -1 * vd.objs[1]
	end
	if vd.objSenses[2] == :Max
		vd.objSenses[2] = :Min
		vd.objs[2] = -1 * vd.objs[2]
	end


	z1, z2 = vd.objs
	ind_INT = find(x->x==:Bin, vm.colCat)
	ind_CONT = find(x->x==:Cont, vm.colCat)

	vm_lex = copy(vm)
	status = solve(vm_lex, method=:lex)
	status == :Optimal || return status
	zUL, zLR = getY_N(vm_lex)
	@show zUL, zLR

	z = (z1 / (zLR[1]-zUL[1])) + (z2 / (zUL[2]-zLR[2])) - (zUL[1] / (zLR[1] - zUL[1])) - (zLR[2] / (zUL[2] - zLR[2]))
	
	JuMP.setobjective(vm, :Min, z.aff)

	res = solve(vm, ignore_solve_hook=true, relaxation=true)

	zlp = getobjectivevalue(vm)
	@show zlp
	xlp = vm.colVal
	@show xlp

	root = Node(vm, zlp, xlp, [], [], [])

	return vm, root
end

end # module
