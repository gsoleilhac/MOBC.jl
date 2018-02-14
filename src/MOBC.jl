module MOBC
using JuMP, vOptGeneric, PyPlot, NSGAII, Suppressor

@enum NODE_STATUS SOLVED BOUNDED PENDING

mutable struct Node
	m::JuMP.Model
	z::Float64
	x::Vector{Float64}
	f0::Vector{Int}
	f1::Vector{Int}
	seen::Vector{Vector{Any}}
	Node(m,f0,f1,seen) = new(m,NaN,[],f0,f1,seen)
	Node(m, z, x, f0, f1, seen) = new(m,z,x,f0,f1,seen)
end

Base.copy(n::Node) = Node(copy(n.m), NaN, [], copy(n.f0), copy(n.f1), copy(n.seen))

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

isonlybinary(m, obj) = all(m.colCat[map(x->getfield(x, :col), obj.aff.vars)] .== :Bin)
hasinteger(m, obj) = any(m.colCat[map(x->getfield(x, :col), obj.aff.vars)] .== :Int)
isbinary(x, inds) = all(y-> y==1. || y==0., view(x, inds))

#Check against all local nadirs, ideally : should just check against the worst
function isfathomable(x, incumbent::AbstractVector)
	res = true
	for i = 1:length(incumbent)-1
		nadir = incumbent[i][2] + incumbent[i+1][1]
		x <= nadir && (res = false)
	end
	res
end

function isdominated(z1, z2, LN)
	for (v1, v2) in LN
		if v1 < z1 && v2 <= z2 || v1 <= z1 && v2 < z2 
			return true
		end
	end
	false
end

function integerbranch(n, inds)

	newnode = copy(n)

	inds0 = find(x->x==0., view(n.x, inds))
	inds1 = find(x->x==1., view(n.x, inds))

	#NO-GOOD CSTR
	@constraint(newnode.m, sum(JuMP.Variable(newnode.m, j) for j in inds0) + sum(1 - JuMP.Variable(newnode.m, j) for j in inds1) >= 1)

	return newnode
end

function evaluate(x, cst, coeffs, vars)::Float64
	res = cst
	for i = 1:length(coeffs)
	    @inbounds res += coeffs[i] * x[vars[i]]
	end
	res
end
evaluate(x, obj) = evaluate(x, obj.aff.constant, obj.aff.coeffs, map(x->getfield(x, :col), obj.aff.vars))

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
	inds_BIN = find(x->x==:Bin, vm.colCat)
	inds_CONT = find(x->x==:Cont, vm.colCat)

	vm_lex = copy(vm)
	status = solve(vm_lex, method=:lex)
	status == :Optimal || return status
	zUL, zLR = getY_N(vm_lex)
	@show zUL, zLR

	Ƶ1 = (z1 - zUL[1]) / (zLR[1] - zUL[1])
	Ƶ2 = (z2 - zLR[2]) / (zUL[2] - zLR[2])
	Ƶ = Ƶ1 + Ƶ2
	
	JuMP.setobjective(vm, :Min, Ƶ.aff)



	# x1 = [getvalue(JuMP.Variable(vm_lex, i), 1) for i = 1:vm.numCols]
	# x2 = [getvalue(JuMP.Variable(vm_lex, i), 2) for i = 1:vm.numCols]
	# modelnsga = copy(vm)
	# modelnsga.ext[:vOpt].objs = [Ƶ1, Ƶ2]
	# modelnsga.ext[:vOpt].objSenses = [:Min, :Min]
	# ns = nsga(150, 500, modelnsga, seed=[x1, x2], pmut=0.3)
	# LN = sort!([Tuple(x.y) for x in ns], by = x -> x[1])



	LN = [(0., 1.), (1., 0.)]
	S = [Node(vm, [], [], [])]

	cpt = 0
	while !isempty(S)
		branch(S, LN, inds_BIN, Ƶ1, Ƶ2)
		cpt += 1
	end

	println()
	@show cpt
	return LN, S, inds_BIN, Ƶ1, Ƶ2
end


function branch(S::AbstractVector{Node}, LN::Vector{Tuple{Float64, Float64}}, inds_BIN, obj1, obj2)

	n = pop!(S)

	res = @suppress solve(n.m, ignore_solve_hook=true, relaxation=true)
	res != :Optimal && return

	n.z = getobjectivevalue(n.m)
	n.x = n.m.colVal
	isfathomable(n.z, LN) && return

	if isbinary(n.x, inds_BIN)
		z1, z2 = evaluate(n.x, obj1), evaluate(n.x, obj2)
		if !isdominated(z1, z2, LN)
			push!(LN, (z1, z2))
			inds_delete = find(x->isdominated(x[1],x[2],LN), LN)
			deleteat!(LN, inds_delete)
			sort!(LN, by = x -> first(x))
			clf()
			plot(map(x->x[1], LN), map(x->x[2], LN), "rx")
		end
		newnode = integerbranch(n, inds_BIN)
		push!(S, newnode)
	else
		n1, n2 = branch(n, inds_BIN)
		push!(S, n1)
		push!(S, n2)

	end
end

function branch(n, inds)
	i = findfirst(x -> x!=0. && x!=1., view(n.x, inds))
	@assert i!= 0
	n1, n2 = copy(n), copy(n)
	n1.m.colUpper[i] = 0. ; @assert n1.m.colLower[i] == 0. "i=$i, colLower=$(n1.m.colLower[i])"
	n2.m.colLower[i] = 1. ; @assert n2.m.colUpper[i] == 1. "i=$i, colUpper=$(n1.m.colUpper[i])"
	push!(n1.f0, i)
	push!(n2.f1, i)
	n1, n2
end


end # module
