module MOBC
using JuMP, vOptGeneric, PyPlot, NSGAII, Suppressor, NamedTuples, Nullables

export solve_BC

include("Node.jl")
include("Plots.jl")
include("NonDomPoints.jl")

isonlybinary(m, obj) = all(m.colCat[map(x->getfield(x, :col), obj.aff.vars)] .== :Bin)
isbinary(x) = all(y-> y==1. || y==0., x)
dominates(::Type{Min}, x, y) = x[1] < y[1] && x[2] <= y[2] || x[1] <= y[1] && x[2] < y[2]
dominates(::Type{Max}, x, y) = x[1] > y[1] && x[2] >= y[2] || x[1] >= y[1] && x[2] > y[2]

evaluate(x, obj) = evaluate(x, obj.aff.constant, obj.aff.coeffs, map(x->getfield(x, :col), obj.aff.vars))
function evaluate(x, cst, coeffs, vars)::Float64
	res = cst
	for i = 1:length(coeffs)
	    @inbounds res += coeffs[i] * x[vars[i]]
	end
	res
end

function solve_BC(vm, limit=500)

	#Asserts
	vd = getvOptData(vm)
	@assert length(vd.objs) == 2
	@assert isonlybinary(vm, vd.objs[1]) && isonlybinary(vm, vd.objs[2])

	@assert vd.objSenses[1] == vd.objSenses[2]
	sense = vd.objSenses[1] == :Max ? Max : Min
	JuMP.setobjectivesense(vm, vd.objSenses[1])

	z1, z2 = vd.objs
	
	#Calculate the convex relaxation
	status = solve(vm, method=:dicho)
	status == :Optimal || return status

	YN_convex = getY_N(vm)
	XE_convex = [[getvalue(JuMP.Variable(vm, i), j) for i = 1:vm.numCols] for j = 1:length(YN_convex)]
	xUL, zUL = XE_convex[1], YN_convex[1]
	xLR, zUL = XE_convex[end], YN_convex[end]
	length(XE_convex) > 1 || return @NT(YN = YN_convex, XE = XE_convex)

	# Use NSGAII to get a better incumbent
	modelnsga = copy(vm)
	modelnsga.ext[:vOpt].objs = [z1, z2]
	# ns = unique(x->x.y, nsga(Val(2), 100, 500, modelnsga, 3, seed=[x_1_2, x_2_1], pmut=0.3))
	ns = unique(x->x.y, nsga(Val(2), 200, 500, modelnsga, seed=XE_convex, pmut=0.3))
	sort!(ns, by=x->x.y[1])
	LN_NSGA = NonDomPoints(sense, map(x->x.pheno, ns), map(x->Tuple(x.y), ns))

	yn = getY_N(vm)
	@assert issorted(yn, by = x->x[1])
	xe = [[getvalue(JuMP.Variable(vm, i), j) for i = 1:vm.numCols] for j = 1:length(yn)]

	resxe = []
	resyn = []

	for i = 1:length(yn)-1

		# println("################################")
		m = copy(vm)
		# LNGlobal = NonDomPoints(sense, xe, Tuple.(yn))
		LNGlobal = LN_NSGA
		LN = NonDomPoints(sense, [xe[i], xe[i+1]], [Tuple(yn[i]), Tuple(yn[i+1])])

		# @show LN.yn
		for j = 1:length(LN_NSGA.yn)
			if LN_NSGA.yn[j][1] > LN.yn[1][1] && LN_NSGA.yn[j][2] > LN.yn[end][2]
				push!(LN, LN_NSGA.xe[j], LN_NSGA.yn[j])
			end
		end
		# @show LN.yn
		Ƶ1, Ƶ2 = copy(z1, m), copy(z2, m)
		Ƶ = LN.λ[1]*Ƶ1 + LN.λ[2]*Ƶ2
		JuMP.setobjective(m, vd.objSenses[1], Ƶ.aff)
		if sense == Max
			@constraint(m, Ƶ1.aff >= LN.yn[1][1] + 1)
			@constraint(m, Ƶ2.aff >= LN.yn[end][2] + 1)
		else
			@constraint(m, Ƶ1.aff <= LN.yn[end][1] - 1)
			@constraint(m, Ƶ2.aff <= LN.yn[1][2] - 1)
		end

		#Stack of nodes to evaluate
		S = [Node(m)]

		#Solve while there are nodes to process
		cpt = 0
		while !isempty(S) && cpt < limit
			sort!(S, by = x->x.zparent, rev=true)
			process_node(S, sense, LN, Ƶ1, Ƶ2, LNGlobal)
			cpt += 1
		end

		cpt == limit && println("cpt limit reached")

		append!(resxe, LN.xe)
		append!(resyn, LN.yn)

		# return [1,2],[3,4]
	end


	resxe = unique(resxe)
	YN = [(evaluate(x, vd.objs[1]), evaluate(x, vd.objs[2])) for x in resxe]

	# @show cpt
	return @NT(YN = YN, XE = resxe)
end


function process_node(S::AbstractVector{Node}, sense, LN, obj1, obj2, LNGlobal)

	n = pop!(S)
	n.zparent != Inf && isfathomable(n.zparent, LN) && return
	res = @suppress solve(n.m, ignore_solve_hook=true, relaxation=true)
	res != :Optimal && return
	n.z = getobjectivevalue(n.m)
	n.x = n.m.colVal

	z1, z2 = evaluate(n.x, obj1), evaluate(n.x, obj2)

	if isfathomable(n.z, LN)
		# println("fathomable")
		# plot_int_found(LN, LNGlobal, z1, z2, sleeptime=0.01)
		return
	end
	
	if isbinary(n.x)
		# println("isbinary")
		if !isdominated(z1, z2, LN)
			# println("new int solution found")
			push!(LN, n.x, (z1, z2))
			push!(S, integerbranch(n))

			# plot_int_found(LN, LNGlobal, z1, z2)
		else
			# println("dominated but not fathomable : paretobranch")
			append!(S, paretobranch(sense, n, z1, z2, LN, obj1, obj2))
		end
	else
		# println("not binary")
		if isdominated(z1, z2, LN)
			# println("dominated : paretobranch")
			# plot_int_found(LN, LNGlobal, z1, z2)
			for node in paretobranch(sense, n, z1, z2, LN, obj1, obj2)
				push!(S, node)
			end
		else
			# println("not dominated : basic branch")
			# plot_int_found(LN, LNGlobal, z1, z2, sleeptime=0.01)
			n1, n2 = basicbranch(n)
			push!(S, n1)
			push!(S, n2)
		end
	end

end


function integerbranch(n)

	newnode = copy(n)
	inds0 = find(x->x==0., n.x)
	inds1 = find(x->x==1., n.x)
	#NO-GOOD CSTR
	@constraint(newnode.m, sum(JuMP.Variable(newnode.m, j) for j in inds0) + sum(1 - JuMP.Variable(newnode.m, j) for j in inds1) >= 1)
	return newnode
end


function paretobranch(sense::Type{Min}, n, z1, z2, LN, obj1, obj2)
	@show "paretomin"
	λ = LN.λ
	WS = sum((z1,z2).*λ)
	i = findfirst(x -> dominates(sense, x, (z1,z2)), LN.yn)
	#Find the first dominated nadir going left from (z1,z2)
	left = i-1
	while left > 0 && sum(nadirs(LN)[left].*λ) <= WS
		left -= 1
	end

	#Find the first dominated nadir going right from (z1,z2)
	right = i
	while right < length(LN) && sum(nadirs(LN)[right].*λ) <= WS
		right += 1
	end

	boundz1 = left > 0 ? nadirs(LN)[left][1] : NaN
	boundz2 = right < length(LN) ? nadirs(LN)[right][2] : NaN

	res = Node[]

	if !isnan(boundz1)
		n1 = copy(n)
		@constraint(n1.m, copy(obj1, n1.m).aff <= boundz1 - 1)
		push!(res, n1)
	end
	if !isnan(boundz2)
		n2 = copy(n)
		@constraint(n2.m, copy(obj2, n2.m).aff <= boundz2 - 1)
		push!(res, n2)
	end
	@assert length(res) >= 1
	res
end

function paretobranch(sense::Type{Max}, n, z1, z2, LN, obj1, obj2)
	λ = LN.λ
	WS = sum((z1,z2).*λ)
	i = findfirst(x -> dominates(sense, x, (z1,z2)), LN.yn)
	#Find the first dominated nadir going left from (z1,z2)
	left = i-1
	while left > 0 && sum(nadirs(LN)[left].*λ) >= WS
		left -= 1
	end

	#Find the first dominated nadir going right from (z1,z2)
	right = i
	while right < length(LN) && sum(nadirs(LN)[right].*λ) >= WS
		right += 1
	end

	boundz1 = right < length(LN) ? nadirs(LN)[right][1] : NaN
	boundz2 = left > 0 ? nadirs(LN)[left][2] : NaN

	res = Node[]

	if !isnan(boundz1)
		n1 = copy(n)
		@constraint(n1.m, copy(obj1, n1.m).aff >= boundz1 + 1)
		push!(res, n1)
	end
	if !isnan(boundz2)
		n2 = copy(n)
		@constraint(n2.m, copy(obj2, n2.m).aff >= boundz2 + 1)
		push!(res, n2)
	end
	# plot_pareto_branch(LN, z1, z2, boundz1, boundz2, sleeptime=0.01)
	@assert length(res) >= 1
	res
end


function basicbranch(n)
	i = findfirst(x -> !isinteger(x), n.x)
	n1, n2 = copy(n), copy(n)
	setlowerbound(JuMP.Variable(n1.m, i), 1.)
	setupperbound(JuMP.Variable(n2.m, i), 0.)
	n1, n2
end


end # module
