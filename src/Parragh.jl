
struct Point
    x::Float64
    y::Float64
end

import Base.first ; first(p::Point) = p.x
import Base.last ; last(p::Point) = p.y

struct Segment
    p1::Point
    p2::Point
    c::Point
    a::Float64
    b::Float64
end

function Segment(p1::Point, p2::Point, c::Point)
	a = (p2.y - p1.y) / (p2.x - p1.x)
	b = p1.y - a*p1.x
	Segment(p1, p2, c, a, b)
end
Segment(p1::Point, p2::Point, s::Type{Max}) = Segment(p1, p2, Point(p1.x, p2.y))
Segment(p1::Point, p2::Point, s::Type{Min}) = Segment(p1, p2, Point(p2.x, p1.y))

function filterSegment(::Type{Min}, s::Segment, u::Point)
	S = Segment[]
	p1, p2, c = s.p1, s.p2, s.c
	if u.y >= s.a * u.x + s.b || u.y >= p1.y || u.x >= p2.x
		if u.x <= p1.x
			push!(S, Segment(p1, p2, Point(c.x, min(c.y, u.y))))
		elseif u.y <= p2.y
			push!(S, Segment(p1, p2, Point(min(c.x, u.x), c.y)))
		else
			push!(S, s)
		end
	else
		if u.x > p1.x
			push!(S, Segment(p1, Point(u.x,s.a * u.x + s.b), Point(u.x,c.y)))
		end
		if u.y > p2.y
			push!(S, Segment(Point((u.y - s.b) / s.a, u.y), p2,Point(c.x, u.y)))
		end
	end
	return S
end

function filterLB(LB, UB, debug=false)
	S = LB
	for u in UB
		@show u
		S_prime = Segment[]
		for s in S
			append!(S_prime, filterSegment(Min, s, u))
		end
		S = S_prime
		if debug
			@show S
			clf()
			plotdualbound(S)
			plot(first.(UB), last.(UB), "ko")
			sleep(1)
			println()
		end
	end
	return S
end




function toSegmentList(lb, nadir)
	res = Segment[]
	for i = 1:length(lb)-1
		c = Point(i==length(lb)-1 ? nadir[1] : lb[i+1][1], nadir[2])
		push!(res, Segment(Point(lb[i]...), Point(lb[i+1]...), c))
	end
	res
end

toPointList(ub) = [Point(u...) for u in ub]



function solve_parragh(vm, limit=500 ;  showplot = false, docovercuts = true , global_branch = false, use_nsga = true, global_nsga = true, lift_covers = false)

	vd = getvOptData(vm)
	@assert length(vd.objs) == 2
	@assert all(isonlybinary.(vm, vd.objs))
	@assert vd.objSenses[1] == vd.objSenses[2]

	sense = vd.objSenses[1] == :Max ? Max : Min
	JuMP.setobjectivesense(vm, vd.objSenses[1])
	cstrData = ConstraintData(vm)

	z1, z2 = vd.objs
	status = global_branch ? solve(vm, method=:lexico) : solve(vm, method=:dicho, round_results=true)
	status == :Optimal || return status

	YN_convex = map(x->round.(Int, x), getY_N(vm))
	XE_convex = [[getvalue(JuMP.Variable(vm, i), j) for i = 1:vm.numCols] for j = 1:length(YN_convex)]

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
			ns = nsga(100, 1000, modelnsga, seed=XE_convex, pmut=0.3, showprogress=false)
		else
			ns = union((nsga(100, 200, modelnsga, seed=[XE_convex[i], XE_convex[i+1]], pmut=0.3, showprogress=false) for i = 1:length(XE_convex)-1)...)
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
		# Ƶ = LN.λ[1]*Ƶ1 + LN.λ[2]*Ƶ2
		# JuMP.setobjective(m, vd.objSenses[1], Ƶ.aff)

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
			sort!(S, by = x->x.zparent, rev=(sense==Min))
			process_node_parragh(pop!(S), S, sense, LN, Ƶ1, Ƶ2, LNGlobal, cstrData, showplot, docovercuts, lift_covers)
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

function process_node_parragh(n::Node, S, sense, LN, obj1, obj2, LNGlobal, cstrData, showplot, docovercuts, lift_covers)
	#n.zparent != Inf && isfathomable(n.zparent, LN) && return
	res = @suppress solve(n.m, method=:dicho, relaxation=true)
	res != :Optimal && return

	YN_relax = getY_N(n.m)
    XE_convex = [[getvalue(JuMP.Variable(vm, i), j) for i = 1:vm.numCols] for j = 1:length(YN_convex)]
    

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
			showplot && plot_int_found(LN, LNGlobal, z1, z2)
		end
	else
		if isweaklydominated(z1, z2, LN)
			# println("not binary, dominated : paretobranch")
			showplot && plot_int_found(LN, LNGlobal, z1, z2, marker="r.")
			for node in paretobranch(sense, n, z1, z2, LN, obj1, obj2, showplot)
				push!(S, node)
			end
		else
			# println("not binary, not dominated : basic branch")
			showplot && plot_int_found(LN, LNGlobal, z1, z2, sleeptime=0.01, marker="k.")
			if docovercuts && n.nbcover <= 10 && find_cover_cuts(n, cstrData, lift_covers)
				n.nbcover += 1
				process_node_parragh(n, S, sense, LN, obj1, obj2, LNGlobal, cstrData, showplot, docovercuts, lift_covers)
			else
				n1, n2 = basicbranch(n)
				push!(S, n1)
				push!(S, n2)
			end
		end
	end
end