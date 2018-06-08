function solve_stidsen(model, limit=Inf ;  showplot = false, docovercuts = true , global_branch = false, use_nsga = true, global_nsga = true, lift_covers = false)

	vm = copy(model)
	vd = getvOptData(vm)
	@assert length(vd.objs) == 2
	@assert all(isonlybinary.(vm, vd.objs))
	@assert vd.objSenses[1] == vd.objSenses[2]

	sense = vd.objSenses[1] == :Max ? Max : Min
	JuMP.setobjectivesense(vm, vd.objSenses[1])
	cstrData = ConstraintData(vm)

	z1, z2 = vd.objs
	status = global_branch ? solve(vm, method=:lexico, verbose=false) : solve(vm, method=:dicho, round_results=true)
	status == :Optimal || return status

	YN_convex = map(x->round.(Int, x), getY_N(vm))
	XE_convex = [[round(getvalue(JuMP.Variable(vm, i), j)) for i = 1:vm.numCols] for j = 1:length(YN_convex)]

	if length(YN_convex) == 1 || YN_convex[1] == YN_convex[end]
		return @NT(YN = YN_convex, XE = XE_convex, nodes=0)
	end

	if !issorted(YN_convex, by = first)
		s = sortperm(YN_convex, by = first)
        YN_convex = YN_convex[s]
        XE_convex = XE_convex[s]
	end

	if use_nsga
		modelnsga = copy(vm)
		modelnsga.ext[:vOpt].objs = [z1, z2]

		if global_nsga
			ns = nsga_binary(50, 500, modelnsga, seed=XE_convex, pmut=0.3, showprogress=false)
		else
			ns = union((nsga_binary(50, 200, modelnsga, seed=[XE_convex[i], XE_convex[i+1]], pmut=0.3, showprogress=false) for i = 1:length(XE_convex)-1)...)
			#if we do one nsga per triangle, we can have dominated solutions here.
			NSGAII.fast_non_dominated_sort!(ns, sense==Max ? NSGAII.Max() : NSGAII.Min())
			filter!(x->x.rank == 1, ns)
		end
		
		ns = unique(x->x.y, ns)
		sort!(ns, by=x->x.y[1])
		LN_NSGA = NonDomPoints(sense, map(x->x.pheno, ns), map(x->Tuple(x.y), ns))
	end

	resxe = Vector{Int}[]
	resyn = Tuple{Int,Int}[]

	nbNodesTotal = 0

	for i = 1:length(YN_convex)-1

		m = copy(vm)
		LN = NonDomPoints(sense, [XE_convex[i], XE_convex[i+1]], [Tuple(YN_convex[i]), Tuple(YN_convex[i+1])])
		LNGlobal = use_nsga ? LN_NSGA : LN #for plots
		
		if abs(LN.yn[1][1]-LN.yn[end][1]) == 1 || abs(LN.yn[1][2]-LN.yn[end][2]) == 1
			append!(resxe, LN.xe)
			append!(resyn, LN.yn)
			continue
		end

		if use_nsga
			for j = 1:length(LN_NSGA.yn)
				if LN.yn[1][1] < LN_NSGA.yn[j][1] < LN.yn[end][1] 
					push!(LN, LN_NSGA.xe[j], LN_NSGA.yn[j])
				end
			end
		end

		Ƶ1, Ƶ2 = copy(z1, m), copy(z2, m)
		Ƶ = LN.λ[1]*Ƶ1 + LN.λ[2]*Ƶ2
		JuMP.setobjective(m, vd.objSenses[1], Ƶ.aff)

		if sense == Max
			bound1 = LN.yn[1][1] + 1
			bound2 = LN.yn[end][2] + 1
			cstr1 = @constraint(m, Ƶ1.aff >= bound1)
			cstr2 = @constraint(m, Ƶ2.aff >= bound2)
			# @constraint(m, Ƶ.aff <= LN.λ[1]*LN.yn[1][1] + LN.λ[2]*LN.yn[1][2])
		else
			bound1 = LN.yn[end][1] - 1
			bound2 = LN.yn[1][2] - 1
			cstr1 = @constraint(m, Ƶ1.aff <= bound1)
			cstr2 = @constraint(m, Ƶ2.aff <= bound2)
		end

		#Stack of nodes to evaluate
		S = [Node(m, bound1, bound2, cstr1, cstr2)]

		#Solve while there are nodes to process
		cpt = 0
		while !isempty(S) && cpt < limit
			sort!(S, by = x->x.zparent, rev=(sense==Min))
			process_node_stidsen(pop!(S), S, sense, LN, Ƶ1, Ƶ2, LNGlobal, cstrData, showplot, docovercuts, lift_covers)
			nbNodesTotal += 1
			cpt += 1
		end
		
		cpt == limit && println("node limit reached")

		append!(resxe, LN.xe)
		append!(resyn, LN.yn)

	end


	resxe = unique(resxe)
	YN = [(evaluate(x, vd.objs[1]), evaluate(x, vd.objs[2])) for x in resxe]

	# @show cpt
	return @NT(YN = YN, XE = resxe, nodes = nbNodesTotal)
end

function setRHS!(n::Union{Node, NodeParragh})
	JuMP.setRHS(n.cstr1, n.bound1)
	JuMP.setRHS(n.cstr2, n.bound2)
end

function fix_variables!(n::Union{Node, NodeParragh})
	for i in n.f1
		setlowerbound(JuMP.Variable(n.m, i), 1.)
	end
	for i in n.f0
		setupperbound(JuMP.Variable(n.m, i), 0.)
	end
end

function unfix_variables!(n::Union{Node, NodeParragh})
	for i in n.f1
		setlowerbound(JuMP.Variable(n.m, i), 0.)
	end
	for i in n.f0
		setupperbound(JuMP.Variable(n.m, i), 1.)
	end
end


function process_node_stidsen(n::Node, S, sense, LN, obj1, obj2, LNGlobal, cstrData, showplot, docovercuts, lift_covers)

	n.zparent != Inf && isfathomable(n.zparent, LN) && return

	fix_variables!(n)
	setRHS!(n)
	res = @suppress solve(n.m, ignore_solve_hook=true, relaxation=true)
	unfix_variables!(n)

	res != :Optimal && return
	n.z = round(getobjectivevalue(n.m), 8)
	n.x = round.(n.m.colVal, 8)
	
	z1, z2 = evaluate(n.x, obj1), evaluate(n.x, obj2)

	if isfathomable(n.z, LN)
		showplot && plot_int_found(LN, LNGlobal, z1, z2, marker="k.", sleeptime=0.01)
		return
	end
	
	if isbinary(n.x)
		if isweaklydominated(z1, z2, LN)
			# println("int, dominated but not fathomable : paretobranch")
			showplot && plot_int_found(LN, LNGlobal, z1, z2, marker = "kx")
			append!(S, paretobranch(sense, n, z1, z2, LN, obj1, obj2, showplot))
		else
			# println("new int solution found")
			push!(LN, n.x, (z1, z2))
			# push!(S, integerbranch(n))
			append!(S, paretobranch(sense, n, z1, z2, LN, obj1, obj2, showplot))
			showplot && plot_int_found(LN, LNGlobal, z1, z2, sleeptime=1)
		end
	else
		# println("not binary")
		if isweaklydominated(z1, z2, LN)
			# println("not binary, dominated : paretobranch")
			showplot && plot_int_found(LN, LNGlobal, z1, z2, marker="r.")
			for node in paretobranch(sense, n, z1, z2, LN, obj1, obj2, showplot)
				push!(S, node)
			end
		else
			# println("not binary, not dominated : basic branch")
			showplot && plot_int_found(LN, LNGlobal, z1, z2, sleeptime=0.01, marker="k.")
			if docovercuts && n.nbcover <= 7
				if find_cover_cuts(n, cstrData, lift_covers)
					n.nbcover += 1
				else
					n.nbcover = 8
				end
				return process_node_stidsen(n, S, sense, LN, obj1, obj2, LNGlobal, cstrData, showplot, docovercuts, lift_covers)
			else
				n1, n2 = basicbranch(n)
				push!(S, n1)
				push!(S, n2)
			end
		end
	end
end

function paretobranch(sense::Type{Min}, n, z1, z2, LN, obj1, obj2, showplot)
	λ = LN.λ
	WS = sum((z1,z2).*λ)
	i = findfirst(x -> weakly_dominates(sense, x, (z1,z2)), LN.yn)
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
		n1 = shallowcopy(n)
		n1.bound1 = boundz1
		# @constraint(n1.m, copy(obj1, n1.m).aff <= boundz1)
		push!(res, n1)
	end
	if !isnan(boundz2)
		n2 = shallowcopy(n)
		n2.bound2 = boundz2
		# @constraint(n2.m, copy(obj2, n2.m).aff <= boundz2)
		push!(res, n2)
	end
	showplot && plot_pareto_branch(LN, z1, z2, boundz1, boundz2, sleeptime=0.01)
	#@assert length(res) >= 1
	res
end

function paretobranch(sense::Type{Max}, n, z1, z2, LN, obj1, obj2, showplot)
	λ = LN.λ
	WS = sum((z1,z2).*λ)
	i = findfirst(x -> weakly_dominates(sense, x, (z1,z2)), LN.yn)
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

	boundz1 = right < length(LN) ? nadirs(LN)[right][1] + 0.5 : NaN
	boundz2 = left > 0 ? nadirs(LN)[left][2] + 0.5 : NaN

	res = Node[]

	if !isnan(boundz1)
		n1 = shallowcopy(n)
		n1.bound1 = boundz1
		# @constraint(n1.m, copy(obj1, n1.m).aff >= boundz1)
		push!(res, n1)
	end
	if !isnan(boundz2)
		n2 = shallowcopy(n)
		n2.bound2 = boundz2
		# @constraint(n2.m, copy(obj2, n2.m).aff >= boundz2)
		push!(res, n2)
	end
	showplot && plot_pareto_branch(LN, z1, z2, boundz1, boundz2, sleeptime=0.01)
	#@assert length(res) >= 1
	res
end

function basicbranch(n)
	i = indmax(map(x->x>0.5 ? 1-x : x, n.x))
	n1, n2 = shallowcopy(n), shallowcopy(n)
	push!(n2.f0, i)
	push!(n1.f1, i)
	n1, n2
end