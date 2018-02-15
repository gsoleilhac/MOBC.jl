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


nadirs(incumbent) = [(incumbent[i+1][1], incumbent[i][2]) for i = 1:length(incumbent)-1]

#Check against all local nadirs, ideally : should just check against the worst
function isfathomable(x, incumbent::AbstractVector)
	res = true
	for n in nadirs(incumbent)
		nadir = n[1]+n[2]
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

dominates(x, y) = x[1] < y[1] && x[2] <= y[2] || x[1] <= y[1] && x[2] < y[2]

function paretobranch(n, z1, z2, LN, obj1, obj2)
	i = findfirst(x -> dominates(x, (z1,z2)), LN)
	left = i-1
	while left > 0 && sum(nadirs(LN)[left]) <= z1+z2
		left -= 1
	end
	right = i
	while right < length(LN) && sum(nadirs(LN)[right]) <= z1+z2
		right += 1
	end

	boundz1 = left > 0 ? nadirs(LN)[left][1] : 0.
	boundz2 = right < length(LN) ? nadirs(LN)[right][2] : 0.


	### PLOT ###
	clf()
	LNPLOT=[]
	for i = 1:length(LN)
		push!(LNPLOT, LN[i])
		i!=length(LN) && push!(LNPLOT, nadirs(LN)[i])
	end
	plot(map(x->x[1], LNPLOT), map(x->x[2], LNPLOT), "r--")
	plot(map(x->x[1], LN), map(x->x[2], LN), "ro")
	plot(map(x->x[1], nadirs(LN)), map(x->x[2], nadirs(LN)), "rs", markersize="2")
	plot([0, z1+z2], [z1+z2, 0], "b--")
	plot([z1], [z2], "bo")
	plot([boundz1, boundz1], [0, 1], "g--")
	plot([0, 1], [boundz2, boundz2], "g--")
	# readline(STDIN)
	### PLOT ### 

	res = Node[]

	if boundz1 != 0.
		n1 = copy(n)
		@constraint(n1.m, copy(obj1, n1.m).aff <= boundz1 - 1e-6)
		if boundz1 ≈ z1
			#Add a no-good constraint to avoid cycling
			inds0 = find(x->x==0., n.x)
			inds1 = find(x->x==1., n.x)
			@constraint(n1.m, sum(JuMP.Variable(n1.m, j) for j in inds0) + sum(1 - JuMP.Variable(n1.m, j) for j in inds1) >= 1)
		end
		push!(res, n1)
	end
	if boundz2 != 0.
		n2 = copy(n)
		@constraint(n2.m, copy(obj2, n2.m).aff <= boundz2 - 1e-6)
		if boundz2 ≈ z2
			#Add a no-good constraint to avoid cycling
			inds0 = find(x->x==0., n.x)
			inds1 = find(x->x==1., n.x)
			@constraint(n2.m, sum(JuMP.Variable(n2.m, j) for j in inds0) + sum(1 - JuMP.Variable(n2.m, j) for j in inds1) >= 1)
		end
		push!(res, n2)
	end

	res
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


	if (objSense1 = vd.objSenses[1]) == :Max
		vd.objSenses[1] = :Min
		vd.objs[1] = -1 * vd.objs[1]
	end
	if (objSense2 = vd.objSenses[2]) == :Max
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



	x1 = [getvalue(JuMP.Variable(vm_lex, i), 1) for i = 1:vm.numCols]
	x2 = [getvalue(JuMP.Variable(vm_lex, i), 2) for i = 1:vm.numCols]
	modelnsga = copy(vm)
	modelnsga.ext[:vOpt].objs = [Ƶ1, Ƶ2]
	modelnsga.ext[:vOpt].objSenses = [:Min, :Min]
	ns = nsga(200, 2000, modelnsga, 3, seed=[x1, x2], pmut=0.3)
	LN = sort!([Tuple(x.y) for x in ns], by = x -> x[1])



	# LN = [(0., 1.), (1., 0.)]
	S = [Node(vm, [], [], [])]

	cpt = 0
	while !isempty(S) && cpt < 10000
		process_node(S, LN, inds_BIN, Ƶ1, Ƶ2)
		cpt += 1
		cpt % 50 == 1 && @show cpt
	end

	println()
	@show cpt
	return LN, S, inds_BIN, Ƶ1, Ƶ2
end


function process_node(S::AbstractVector{Node}, LN::Vector{Tuple{Float64, Float64}}, inds_BIN, obj1, obj2)

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
			inds_delete = find(x->dominates((z1,z2), x), LN)
			deleteat!(LN, inds_delete)
			sort!(LN, by = x -> first(x))

			### PLOT ###
			clf()
			LNPLOT=[]
			for i = 1:length(LN)
				push!(LNPLOT, LN[i])
				i!=length(LN) && push!(LNPLOT, nadirs(LN)[i])
			end
			plot(map(x->x[1], LNPLOT), map(x->x[2], LNPLOT), "r--")
			plot(map(x->x[1], LN), map(x->x[2], LN), "ro")
			plot(map(x->x[1], nadirs(LN)), map(x->x[2], nadirs(LN)), "rs", markersize="2")
			plot([0, z1+z2], [z1+z2, 0], "b--")
			plot([z1], [z2], "bo")
			# readline(STDIN)
			### PLOT ### 

		else

			for node in paretobranch(n, z1, z2, LN, obj1, obj2)
				push!(S, node)
			end

		end
			# newnode = integerbranch(n, inds_BIN)
			# push!(S, newnode)
	else
		z1, z2 = evaluate(n.x, obj1), evaluate(n.x, obj2)
		if isdominated(z1, z2, LN)
			for node in paretobranch(n, z1, z2, LN, obj1, obj2)
				push!(S, node)
			end
		else
			n1, n2 = branch(n, inds_BIN)
			push!(S, n1)
			push!(S, n2)
		end
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
