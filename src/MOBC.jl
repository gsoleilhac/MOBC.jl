module MOBC
using JuMP, vOptGeneric, PyPlot, NSGAII, Suppressor, NamedTuples

include("Node.jl")
include("NonDomPoints.jl")

isonlybinary(m, obj) = all(m.colCat[map(x->getfield(x, :col), obj.aff.vars)] .== :Bin)
hasinteger(m, obj) = any(m.colCat[map(x->getfield(x, :col), obj.aff.vars)] .== :Int)
isbinary(x, inds) = all(y-> y==1. || y==0., view(x, inds))
dominates(x, y) = x[1] < y[1] && x[2] <= y[2] || x[1] <= y[1] && x[2] < y[2]
evaluate(x, obj) = evaluate(x, obj.aff.constant, obj.aff.coeffs, map(x->getfield(x, :col), obj.aff.vars))
function evaluate(x, cst, coeffs, vars)::Float64
	res = cst
	for i = 1:length(coeffs)
	    @inbounds res += coeffs[i] * x[vars[i]]
	end
	res
end

function solve_BC(vm, limit=5000)

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
	xUL = [getvalue(JuMP.Variable(vm_lex, i), 1) for i = 1:vm.numCols]
	xLR = [getvalue(JuMP.Variable(vm_lex, i), 2) for i = 1:vm.numCols]
	zUL, zLR = getY_N(vm_lex)

	Ƶ1 = (z1 - zUL[1]) / (zLR[1] - zUL[1])
	Ƶ2 = (z2 - zLR[2]) / (zUL[2] - zLR[2])
	Ƶ = Ƶ1 + Ƶ2
	
	JuMP.setobjective(vm, :Min, Ƶ.aff)

	modelnsga = copy(vm)
	modelnsga.ext[:vOpt].objs = [Ƶ1, Ƶ2]
	modelnsga.ext[:vOpt].objSenses = [:Min, :Min]
	ns = unique(x->x.y, nsga(100, 500, modelnsga, 3, seed=[xUL, xLR], pmut=0.3))
	sort!(ns, by=x->x.y[1])

	LN = NonDomPoints(map(x->x.pheno, ns), map(x->x.y, ns))
	# LN = [(0., 1.), (1., 0.)]

	S = [Node(vm)]

	cpt = 0
	while !isempty(S) && cpt < limit
		cpt % 1000 == 0 && @show cpt
		sort!(S, by = x->x.zparent, rev=true)
		process_node(S, LN, inds_BIN, Ƶ1, Ƶ2)
		cpt += 1
	end

	#Transform results back to the original problem :
	if objSense1 == :Max
		vd.objSenses[1] = :Max
		vd.objs[1] = -1 * vd.objs[1]
	end
	if objSense2 == :Max
		vd.objSenses[2] = :Max
		vd.objs[2] = -1 * vd.objs[2]
	end
	YN = [(evaluate(x, vd.objs[1]), evaluate(x, vd.objs[2])) for x in LN.xe]

	@show cpt
	return @NT(YN = YN, XE = LN.xe)
end


function process_node(S::AbstractVector{Node}, LN, inds_BIN, obj1, obj2)

	n = pop!(S)
	n.zparent != Inf && isfathomable(n.zparent, LN) && return
	res = @suppress solve(n.m, ignore_solve_hook=true, relaxation=true)
	res != :Optimal && return
	n.z = getobjectivevalue(n.m)
	n.x = n.m.colVal
	isfathomable(n.z, LN) && return

	if isbinary(n.x, inds_BIN)
		z1, z2 = evaluate(n.x, obj1), evaluate(n.x, obj2)
		if !isdominated(z1, z2, LN)
			push!(LN, n.x, (z1, z2))
			push!(S, integerbranch(n, inds_BIN))

			## PLOT ###
			# clf() ; LNPLOT=[]
			# for i = 1:length(LN)
			# 	push!(LNPLOT, LN.yn[i])
			# 	i!=length(LN) && push!(LNPLOT, nadirs(LN)[i])
			# end
			# plot(map(x->x[1], LNPLOT), map(x->x[2], LNPLOT), "r--")
			# plot(map(x->x[1], LN.yn), map(x->x[2], LN.yn), "ro")
			# plot(map(x->x[1], nadirs(LN)), map(x->x[2], nadirs(LN)), "rs", markersize="2")
			# plot([0, z1+z2], [z1+z2, 0], "b--")
			# plot([z1], [z2], "bo")
			## PLOT ### 

		else
			append!(S, paretobranch(n, z1, z2, LN, obj1, obj2))
		end
	else
		z1, z2 = evaluate(n.x, obj1), evaluate(n.x, obj2)
		if isdominated(z1, z2, LN)
			for node in paretobranch(n, z1, z2, LN, obj1, obj2)
				push!(S, node)
			end
		else
			n1, n2 = basicbranch(n, inds_BIN)
			push!(S, n1)
			push!(S, n2)
		end
	end
end


function integerbranch(n, inds)

	newnode = copy(n)
	inds0 = find(x->x==0., view(n.x, inds))
	inds1 = find(x->x==1., view(n.x, inds))
	#NO-GOOD CSTR
	@constraint(newnode.m, sum(JuMP.Variable(newnode.m, j) for j in inds0) + sum(1 - JuMP.Variable(newnode.m, j) for j in inds1) >= 1)
	return newnode
end


function paretobranch(n, z1, z2, LN, obj1, obj2)

	i = findfirst(x -> dominates(x, (z1,z2)), LN.yn)
	#Find the first dominated nadir going left from z1/z2
	left = i-1
	while left > 0 && sum(nadirs(LN)[left]) <= z1+z2
		left -= 1
	end
	#Find the first dominated nadir going right from z1/z2
	right = i
	while right < length(LN) && sum(nadirs(LN)[right]) <= z1+z2
		right += 1
	end

	boundz1 = left > 0 ? nadirs(LN)[left][1] : 0.
	boundz2 = right < length(LN) ? nadirs(LN)[right][2] : 0.


	## PLOT ###
	# clf() ; LNPLOT=[]
	# for i = 1:length(LN)
	# 	push!(LNPLOT, LN.yn[i])
	# 	i!=length(LN) && push!(LNPLOT, nadirs(LN)[i])
	# end
	# plot(map(x->x[1], LNPLOT), map(x->x[2], LNPLOT), "r--")
	# plot(map(x->x[1], LN.yn), map(x->x[2], LN.yn), "ro")
	# plot(map(x->x[1], nadirs(LN)), map(x->x[2], nadirs(LN)), "rs", markersize="2")
	# plot([0, z1+z2], [z1+z2, 0], "b--") ; plot([z1], [z2], "bo")
	# plot([boundz1, boundz1], [0, 1], "g--")
	# plot([0, 1], [boundz2, boundz2], "g--")
	## PLOT ### 

	res = Node[]

	if boundz1 != 0.
		n1 = copy(n)
		@constraint(n1.m, copy(obj1, n1.m).aff <= boundz1 - 1e-4)
		push!(res, n1)
	end
	if boundz2 != 0.
		n2 = copy(n)
		@constraint(n2.m, copy(obj2, n2.m).aff <= boundz2 - 1e-4)
		push!(res, n2)
	end
	res
end


function basicbranch(n, inds)
	i = findfirst(x -> x!=0. && x!=1., view(n.x, inds))
	n1, n2 = copy(n), copy(n)
	setlowerbound(JuMP.Variable(n1.m, i), 1.)
	setupperbound(JuMP.Variable(n2.m, i), 0.)
	n1, n2
end


end # module
