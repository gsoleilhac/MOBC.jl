struct Point
    x::Float64
    y::Float64
end

import Base.show ; Base.show(io::IO, p::Point) = print(io, "P($(round(p.x, 2)), $(round(p.y, 2)))")
import Base.first ; first(p::Point) = p.x
import Base.last ; last(p::Point) = p.y

struct Segment{T}
    p1::Point
    p2::Point
    c::Point
    a::Float64
	b::Float64
	Segment{T}(p1::P, p2::P, c::P, a, b) where {T<:Sense, P<:Point} = new(p1, p2, c, a, b)
	Segment{T}(p1, p2, c) where T = Segment{T}(Point(p1), Point(p2), Point(c))
	function Segment{T}(p1::Point, p2::Point, c::Point) where T
		a = (p2.y - p1.y) / (p2.x - p1.x)
		b = p1.y - a*p1.x
		Segment{T}(p1, p2, c, a, b)
	end
end
Base.show(io::IO, s::Segment) = print(io, "Sgm[$(s.p1), $(s.p2)]")

function filterSegment(s::Segment{Min}, u0::Point)
	u = Point(u0.x - 0.999, u0.y - 0.999)
	S = Segment{Min}[]
	p1, p2, c = s.p1, s.p2, s.c
	if u.y >= s.a * u.x + s.b || u.y >= p1.y || u.x >= p2.x
		if u.x <= p1.x
			push!(S, Segment{Min}(p1, p2, Point(c.x, min(c.y, u.y))))
		elseif u.y <= p2.y
			push!(S, Segment{Min}(p1, p2, Point(min(c.x, u.x), c.y)))
		else
			push!(S, s)
		end
	else
		if u.x > p1.x
			push!(S, Segment{Min}(p1, Point(u.x,s.a * u.x + s.b), Point(u.x,c.y)))
		end
		if u.y > p2.y
			push!(S, Segment{Min}(Point((u.y - s.b) / s.a, u.y), p2,Point(c.x, u.y)))
		end
	end
	return S
end

function coversInteger(s::Segment{Min})
	p1, p2, a, b, c = s.p1, s.p2, s.a, s.b, s.c
	return (p1 == p2 || floor(c.y) >= a * floor(c.x) + b) && floor(c.y) >= ceil(p2.y) && floor(c.x) >= ceil(p1.x)
end


function filterDualBound(dual::Vector{Segment{T}}, primal, debug=false) where T<:Sense
	S = dual
	for u in primal
		# @show u
		S_prime = Segment{T}[]
		for s in S
			append!(S_prime, filterSegment(s, u))
		end
		S = S_prime 
		if debug && u == primal[end]
			# @show S
			clf()
			plotdualbound(S)
			plot(first.(primal), last.(primal), "ko", markersize=4)
			# plot(first.(primal), last.(primal), "k-", linewidth=1)
			nadirsline_x=[primal[1].x]
			nadirsline_y=[primal[1].y]
			for i = 2:length(primal)
				append!(nadirsline_x, [primal[i].x, primal[i].x])
				append!(nadirsline_y, [primal[i-1].y, primal[i].y])
			end
			plot(nadirsline_x, nadirsline_y, "k-", linewidth=1)
			plot(nadirsline_x, nadirsline_y, "ks", markersize=2)
			sleep(0.2)
			# println()
		end
	end
	res = groupby_continuous(S)
	for r in res
		filter!(coversInteger, r)
	end
	filter!(!isempty, res)
	res
end

function toSegmentList(lb, nadir, sense)
	res = Segment{sense}[]
	for i = 1:length(lb)-1
		c = Point(i==length(lb)-1 ? nadir[1] : lb[i+1][1], nadir[2])
		push!(res, Segment{sense}(Point(lb[i]...), Point(lb[i+1]...), c))
	end
	res
end

toPointList(ub) = [Point(u...) for u in ub]

function groupby_continuous(v::Vector{T})::Vector{Vector{T}} where T<:Segment
	isempty(v) && return Vector{Segment}[]
	res = [[v[1]]]
	for i = 2:length(v)
		if v[i].p1 == res[end][end].p2
			push!(res[end], v[i])
		else
			push!(res, [v[i]])
		end
	end
	res
end



function solve_parragh(vm, limit=500 ;  showplot = false, docovercuts = false , global_branch = false, use_nsga = true, global_nsga = true, lift_covers = false)

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
		S = [NodeParragh(m)]
		
		#Solve while there are nodes to process
		cpt = 0
		while !isempty(S) && cpt < limit
			# sort!(S, by = x->x.zparent, rev=(sense==Min))
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

function process_node_parragh(n::NodeParragh, S, sense, LN, obj1, obj2, LNGlobal, cstrData, showplot, docovercuts, lift_covers)

	res = @suppress solve(n.m, method=:dicho, relax=true)
	res != :Optimal && return

	YN_relax = getY_N(n.m)
	# XE_convex = [[getvalue(JuMP.Variable(n.m, i), j) for i = 1:n.m.numCols] for j = 1:length(YN_relax)]
	
	n.x = [[getvalue(JuMP.Variable(n.m, i), j) for i = 1:n.m.numCols] for j = 1:length(YN_relax)]
    
	inds = find(x -> all(isinteger, x), YN_relax)
	for i in inds
		x = [getvalue(JuMP.Variable(n.m, col), i) for col = 1:n.m.numCols]
		z1, z2 = YN_relax[i]
		if all(isinteger, x) && !isweaklydominated(z1, z2, LN)
			push!(LN, x, (z1, z2))
		end
	end

	
	dual = toSegmentList(YN_relax, local_nadir(LN), sense);
	primal = toPointList(LN.yn)
	
	if showplot
		clf() ; 
		for s in dual
			plot([s.p1.x, s.p2.x], [s.p1.y, s.p2.y], "g-", linewidth=1)
			plot([s.p1.x, s.p2.x], [s.p1.y, s.p2.y], "gs", markersize=4)
		end
		plot(first.(primal), last.(primal), "ko", markersize=4)
		nadirsline_x=[primal[1].x]
		nadirsline_y=[primal[1].y]
		for i = 2:length(primal)
			append!(nadirsline_x, [primal[i].x, primal[i].x])
			append!(nadirsline_y, [primal[i-1].y, primal[i].y])
		end
		plot(nadirsline_x, nadirsline_y, "k-", linewidth=1)
		plot(nadirsline_x, nadirsline_y, "ks", markersize=2)
		sleep(0.2)
	end
	
	filteredDual = MOBC.filterDualBound(dual, primal, showplot)

	if isempty(filteredDual)
		#Dual bound set is dominated by the primal set, discard that node
		return
	elseif length(filteredDual) == 1
		if docovercuts && n.nbcover <= 10 && find_cover_cuts(n, cstrData, lift_covers)
			n.nbcover += 1
			process_node_parragh(n, S, sense, LN, obj1, obj2, LNGlobal, cstrData, showplot, docovercuts, lift_covers)
		else
			#branch classique
			n1, n2 = basicbranch_parragh(n)
			push!(S, n1)
			push!(S, n2)
		end
	else
		for segmentList in filteredDual
			push!(S, parragh_branch(n, segmentList, obj1, obj2))
		end
	end

end

function parragh_branch(n, segmentList::Vector{Segment{Min}}, obj1, obj2)
	node = copy(n)
	boundz1 = maximum(x-> x.c.x, segmentList)
	boundz2 = maximum(x-> x.c.y, segmentList)
	@constraint(node.m, copy(obj1, node.m).aff <= boundz1)
	@constraint(node.m, copy(obj2, node.m).aff <= boundz2)
	node
end

function basicbranch_parragh(n)
	i = findfirst(x-> !(x in n.f0) && !(x in n.f1), 1:n.m.numCols)
	n1, n2 = copy(n), copy(n, false)
	setlowerbound(JuMP.Variable(n1.m, i), 1.) ; push!(n1.f1, i)
	setupperbound(JuMP.Variable(n2.m, i), 0.) ; push!(n2.f0, i)
	n1, n2
end