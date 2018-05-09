module MOBC
using JuMP, vOptGeneric, PyPlot, NSGAII, Suppressor, NamedTuples, Nullables, Compat, ProgressMeter

export solve_BC

include("Node.jl")
include("Plots.jl")
include("NonDomPoints.jl")
include("CoverCuts.jl")

isonlybinary(m, obj) = all(m.colCat[map(x->getfield(x, :col), obj.aff.vars)] .== :Bin)
isbinary(x) = all(y-> y==1. || y==0., x)
dominates(::Type{Min}, x, y) = x[1] < y[1] && x[2] <= y[2] || x[1] <= y[1] && x[2] < y[2]
dominates(::Type{Max}, x, y) = x[1] > y[1] && x[2] >= y[2] || x[1] >= y[1] && x[2] > y[2]
weakly_dominates(::Type{Min}, x, y) = x[1] <= y[1] && x[2] <= y[2]
weakly_dominates(::Type{Max}, x, y) = x[1] >= y[1] && x[2] >= y[2]

evaluate(x, obj) = evaluate(x, obj.aff.constant, obj.aff.coeffs, map(x->getfield(x, :col), obj.aff.vars))
function evaluate(x, cst, coeffs, vars)::Float64
	res = cst
	for i = 1:length(coeffs)
	    @inbounds res += coeffs[i] * x[vars[i]]
	end
	res
end

function solve_BC(vm, limit=500 ;  showplot = false, docovercuts = true ; global_nsga = true)

	#Asserts
	vd = getvOptData(vm)
	@assert length(vd.objs) == 2
	@assert all(isonlybinary.(vm, vd.objs))

	@assert vd.objSenses[1] == vd.objSenses[2]
	sense = vd.objSenses[1] == :Max ? Max : Min
	JuMP.setobjectivesense(vm, vd.objSenses[1])

	cstrData = ConstraintData(vm)

	z1, z2 = vd.objs
	
	#Calculate the convex relaxation
	status = solve(vm, method=:dicho)
	status == :Optimal || return status

	YN_convex = map(x->round.(Int, x), getY_N(vm))
	XE_convex = [[getvalue(JuMP.Variable(vm, i), j) for i = 1:vm.numCols] for j = 1:length(YN_convex)]
	@assert issorted(YN_convex, by = x->x[1])

	# Filter X_E and Y_N (to fix a case with values becoming dominated when cast into int)  
	# e.g. : ([683.0, 616.0], [682.9999999999435, 619.0000000000904]) => ([683, 616], [683, 619])
	inds = Int[]
	dominates = sense==Max ? (a,b) -> a[1] >= b[1] && a[2] >= b[2] : (a,b) -> a[1] <= b[1] && a[2] <= b[2]
	for i = 1:length(YN_convex)-1
		if dominates(YN_convex[i], YN_convex[i+1])
			push!(inds, i+1)
		elseif dominates(YN_convex[i+1], YN_convex[i])
			push!(inds, i)
		end
	end
	deleteat!(YN_convex, inds)
	deleteat!(XE_convex, inds)

	xUL, zUL = XE_convex[1], YN_convex[1]
	xLR, zUL = XE_convex[end], YN_convex[end]
	length(XE_convex) > 1 || return @NT(YN = YN_convex, XE = XE_convex, nodes=0)

	# Use NSGAII to get a better incumbent
	modelnsga = copy(vm)
	modelnsga.ext[:vOpt].objs = [z1, z2]

	if global_nsga
		# 1 global nsga
		ns = nsga(100, 1000, modelnsga, seed=XE_convex, pmut=0.3, showprogress=false)
	else
		# 1 nsga per triangle
		ns = union((nsga(100, 200, modelnsga, seed=[XE_convex[i], XE_convex[i+1]], pmut=0.3, showprogress=false) for i = 1:length(XE_convex)-1)...)
		#if we do one nsga per triangle, we can have dominated solutions here.
		NSGAII.fast_non_dominated_sort!(ns, sense==Max ? NSGAII.Max() : NSGAII.Min()) ; filter!(x->x.rank == 1, ns) 
	end
	
	ns = unique(x->x.y, ns)
	
	
	sort!(ns, by=x->x.y[1])
	LN_NSGA = NonDomPoints(sense, map(x->x.pheno, ns), map(x->Tuple(x.y), ns))


	resxe = []
	resyn = []

	nbNodesTotal = 0

	# @showprogress 0.1 
	for i = 1:length(YN_convex)-1

		m = copy(vm)
		LNGlobal = LN_NSGA #for plots
		LN = NonDomPoints(sense, [XE_convex[i], XE_convex[i+1]], [Tuple(YN_convex[i]), Tuple(YN_convex[i+1])])
		
		if abs(LN.yn[1][1]-LN.yn[end][1]) == 1 || abs(LN.yn[1][2]-LN.yn[end][2]) == 1
			append!(resxe, LN.xe)
			append!(resyn, LN.yn)
			continue
		end

		for j = 1:length(LN_NSGA.yn)
			if LN.yn[1][1] < LN_NSGA.yn[j][1] < LN.yn[end][1] 
				push!(LN, LN_NSGA.xe[j], LN_NSGA.yn[j])
			end
		end

		Ƶ1, Ƶ2 = copy(z1, m), copy(z2, m)
		Ƶ = LN.λ[1]*Ƶ1 + LN.λ[2]*Ƶ2
		JuMP.setobjective(m, vd.objSenses[1], Ƶ.aff)
		if sense == Max
			@constraint(m, Ƶ1.aff >= LN.yn[1][1] + 1)
			@constraint(m, Ƶ2.aff >= LN.yn[end][2] + 1)
			# @constraint(m, Ƶ.aff <= LN.λ[1]*LN.yn[1][1] + LN.λ[2]*LN.yn[1][2])
		else
			@constraint(m, Ƶ1.aff <= LN.yn[end][1] - 1)
			@constraint(m, Ƶ2.aff <= LN.yn[1][2] - 1)
		end

		#Stack of nodes to evaluate
		S = [Node(m)]

		
		#Solve while there are nodes to process
		cpt = 0
		while !isempty(S) && cpt < limit
			# sort!(S, by = x->x.zparent, rev=true)
			process_node(pop!(S), S, sense, LN, Ƶ1, Ƶ2, LNGlobal, cstrData, showplot, docovercuts)
			nbNodesTotal += 1
			cpt += 1
		end
		
		cpt == limit && println("cpt limit reached")
		
		append!(resxe, LN.xe)
		append!(resyn, LN.yn)

	end


	resxe = unique(resxe)
	YN = [(evaluate(x, vd.objs[1]), evaluate(x, vd.objs[2])) for x in resxe]

	# @show cpt
	return @NT(YN = YN, XE = resxe, nodes = nbNodesTotal)
end


function process_node(n::Node, S, sense, LN, obj1, obj2, LNGlobal, cstrData, showplot, docovercuts)

	n.zparent != Inf && isfathomable(n.zparent, LN) && return
	res = @suppress solve(n.m, ignore_solve_hook=true, relaxation=true)
	res != :Optimal && return
	n.z = getobjectivevalue(n.m)
	n.x = n.m.colVal

	z1, z2 = evaluate(n.x, obj1), evaluate(n.x, obj2)

	if isfathomable(n.z, LN)
		# println("fathomable")
		# @show n.x, n.z
		# @show sum(worst_nadir(LN).*LN.λ)
		# @show worst_nadir(LN)
		showplot && plot_int_found(LN, LNGlobal, z1, z2, marker="k.", sleeptime=0.01)
		return
	end
	
	if isbinary(n.x)
		if isdominated(z1, z2, LN)
			# println("int, dominated but not fathomable : paretobranch")
			showplot && plot_int_found(LN, LNGlobal, z1, z2, marker = "kx")
			append!(S, paretobranch(sense, n, z1, z2, LN, obj1, obj2, showplot))
		else
			# println("new int solution found")
			push!(LN, n.x, (z1, z2))
			push!(S, integerbranch(n))
			showplot && plot_int_found(LN, LNGlobal, z1, z2)
		end
	else
		if isdominated(z1, z2, LN)
			# println("not binary, dominated : paretobranch")
			showplot && plot_int_found(LN, LNGlobal, z1, z2, marker="r.")
			for node in paretobranch(sense, n, z1, z2, LN, obj1, obj2, showplot)
				push!(S, node)
			end
		else
			# println("not binary, not dominated : basic branch")
			showplot && plot_int_found(LN, LNGlobal, z1, z2, sleeptime=0.01, marker="k.")
			if docovercuts && n.nbcover <= 10 && find_cover_cuts(n, cstrData)
				n.nbcover += 1
				process_node(n, S, sense, LN, obj1, obj2, LNGlobal, cstrData, showplot, docovercuts)
			else
				n1, n2 = basicbranch(n)
				push!(S, n1)
				push!(S, n2)
			end
		end
	end
end


function integerbranch(n)
	newnode = copy(n)
	inds0 = find(x->x==0., n.x)
	inds1 = find(x->x==1., n.x)
	#NO-GOOD CSTR
	@constraint(newnode.m, sum(JuMP.Variable(newnode.m, j) for j in inds0) + sum(1. - JuMP.Variable(newnode.m, j) for j in inds1) >= 1.)
	return newnode
end


function paretobranch(sense::Type{Min}, n, z1, z2, LN, obj1, obj2, showplot)
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

	boundz1 = left > 0 ? nadirs(LN)[left][1] - 1. : NaN
	boundz2 = right < length(LN) ? nadirs(LN)[right][2] - 1. : NaN

	res = Node[]

	if !isnan(boundz1)
		n1 = copy(n)
		@constraint(n1.m, copy(obj1, n1.m).aff <= boundz1)
		push!(res, n1)
	end
	if !isnan(boundz2)
		n2 = copy(n)
		@constraint(n2.m, copy(obj2, n2.m).aff <= boundz2)
		push!(res, n2)
	end
	showplot && plot_pareto_branch(LN, z1, z2, boundz1, boundz2, sleeptime=0.1)
	#@assert length(res) >= 1
	res
end

function paretobranch(sense::Type{Max}, n, z1, z2, LN, obj1, obj2, showplot)
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

	boundz1 = right < length(LN) ? nadirs(LN)[right][1] + 1. : NaN
	boundz2 = left > 0 ? nadirs(LN)[left][2] + 1. : NaN

	res = Node[]

	if !isnan(boundz1)
		n1 = copy(n)
		@constraint(n1.m, copy(obj1, n1.m).aff >= boundz1)
		push!(res, n1)
	end
	if !isnan(boundz2)
		n2 = copy(n)
		@constraint(n2.m, copy(obj2, n2.m).aff >= boundz2)
		push!(res, n2)
	end
	showplot && plot_pareto_branch(LN, z1, z2, boundz1, boundz2, sleeptime=0.01)
	#@assert length(res) >= 1
	res
end


function basicbranch(n)
	i = findfirst(x -> !isinteger(x), n.x)
	n1, n2 = copy(n), copy(n)
	setlowerbound(JuMP.Variable(n1.m, i), 1.) ; push!(n1.f1, i)
	setupperbound(JuMP.Variable(n2.m, i), 0.) ; push!(n2.f0, i)
	n1, n2
end


end # module
