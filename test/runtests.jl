using MOBC, vOptGeneric, CPLEX, PyPlot, Suppressor
using Base.Test

function random_instance(range, n, rng)
	MT = MersenneTwister(rng)
	m = vModel(solver = CplexSolver(CPX_PARAM_SCRIND = 0))
	p1 = rand(MT, range, n)
	p2 = rand(MT, range, n)
	w = [rand(MT, range, n) for i = 1:3]
	c = sum.(w).÷2

	@variable(m, x[1:length(p1)], Bin)
	@addobjective(m, Max, dot(x, p1))
	@addobjective(m, Max, dot(x, p2))
	@constraint(m, dot(x, w[1]) <= c[1])
	@constraint(m, dot(x, w[2]) <= c[2])
	@constraint(m, dot(x, w[3]) <= c[3])
	m
end

function hard_instance(solver=CplexSolver(CPX_PARAM_SCRIND = 0))
    n = 30
    p1 = [85 ,88 ,82 ,64 ,55 ,63 ,77 ,92 ,52 ,64 ,91 ,59 ,90 ,74 ,81 ,80 ,75 ,95 ,52 ,83 ,88 ,50 ,94 ,67 ,91 ,95 ,50 ,74 ,56 ,76]
    p2 = [56 ,58 ,78 ,55 ,69 ,86 ,92 ,84 ,88 ,60 ,62 ,99 ,91 ,62 ,56 ,96 ,75 ,71 ,77 ,87 ,80 ,87 ,72 ,99 ,91 ,74 ,94 ,50 ,98 ,100]
    w1 = [79 ,72 ,60 ,72 ,53 ,53 ,56 ,69 ,57 ,73 ,55 ,66 ,51 ,67 ,74 ,57 ,85 ,98 ,73 ,77 ,58 ,69 ,70 ,76 ,72 ,71 ,87 ,55 ,51 ,83]
    w2 = [67 ,67 ,96 ,88 ,83 ,52 ,60 ,77 ,60 ,75 ,79 ,66 ,71 ,73 ,65 ,58 ,69 ,100 ,67 ,99 ,59 ,69 ,81 ,95 ,59 ,79 ,86 ,56 ,65 ,73]
    w3 = [55 ,72 ,77 ,62 ,51 ,58 ,85 ,76 ,57 ,52 ,83 ,74 ,69 ,89 ,61 ,66 ,66 ,71 ,94 ,61 ,59 ,60 ,69 ,51 ,64 ,84 ,92 ,62 ,67 ,50]
    b = [1019, 1097, 1018]

    m = vModel(solver = solver)
    @variable(m, x[1:n], Bin)
    @addobjective(m, Max, dot(x,p1))
    @addobjective(m, Max, dot(x,p2))
    @constraint(m, dot(x, w1) <= b[1])
    @constraint(m, dot(x, w2) <= b[2])
    @constraint(m, dot(x, w3) <= b[3])

    m
end

function benchmark(itemrange, nrange, runpersize ; args...)
	rng = 0
	res_stidsen=Dict{Int, Vector{Float64}}() ; nb_nodes_stidsen=Dict{Int, Vector{Int}}()
	res_parragh=Dict{Int, Vector{Float64}}() ; nb_nodes_parragh=Dict{Int, Vector{Int}}()
	res_eps=Dict{Int, Vector{Float64}}()
	for n in nrange
		res_stidsen[n]=Float64[] ; res_parragh[n]=Float64[] ; res_eps[n]=Float64[]
		nb_nodes_stidsen[n]=Int[] ; nb_nodes_parragh[n]=Int[]
		for i = 1:runpersize
			m = random_instance(itemrange, n, rng)
			rng += 1
			global m_current = m
			valBC, tBC, _ = @timed solve_stidsen(m, Inf ; args...);
			#valPAR, tPAR, _ = @timed solve_parragh(m, Inf ; args...);
			valPAR, tPAR = valBC, tBC
			status,tEPS,_ = @suppress @timed solve(m, method=:epsilon, round_results=true);
			if !(length(getY_N(m)) == length(unique(valBC.YN)) == length(unique(valPAR.YN)))
				print(m)
				@show valBC.YN
				@show valPAR.YN
				@show getY_N(m)
				figure(3)
				clf()
				plot(map(first, valPAR.YN), map(last,valPAR.YN), "b>", markersize="5", label="parragh")
				plot(map(first, valBC.YN), map(last, valBC.YN), "r<", markersize="5", label="stidsen")
				plot(map(first, getY_N(m)), map(last, getY_N(m)), "go", markersize="2", label="epsilon")
				legend()
				error()
			end

			push!(res_stidsen[n], tBC) ; push!(res_parragh[n], tPAR) ; push!(res_eps[n], tEPS)
			push!(nb_nodes_parragh[n], valPAR.nodes) ; push!(nb_nodes_stidsen[n], valBC.nodes)

			figure(1) ; clf()
			xlabel("n") ; ylabel("t(s)") ; title("time per instance")
			plot([],[],"gx",markersize=3,label="stidsen") ; plot([],[],"bx",markersize=3,label="parragh") ; plot([],[],"k.",markersize=3,label="epsilon")
			for k in keys(res_stidsen)
				plot(fill(k, length(res_stidsen[k])), res_stidsen[k], "gx", markersize=3)
				plot(fill(k, length(res_stidsen[k])), res_parragh[k], "bx", markersize=3)
				plot(fill(k, length(res_stidsen[k])), res_eps[k], "k.", markersize=3)
			end
			legend()

			figure(2) ; clf()
			xlabel("n") ; ylabel("#nodes") ; title("average time and #nodes")
			k = sort(collect(keys(res_stidsen)))
			plot(k, [mean(nb_nodes_stidsen[a]) for a in k], "gs", markersize=3) ; plot(k, [mean(nb_nodes_stidsen[a]) for a in k], "g--", linewidth=1)
			plot(k, [mean(nb_nodes_parragh[a]) for a in k], "bs", markersize=3) ; plot(k, [mean(nb_nodes_parragh[a]) for a in k], "b--", linewidth=1)
			twinx() ; ylabel("t(s)")
			plot([],[],"g--",markersize=3,label="#nodes stidsen") #labels dont work after twinx(), had to put a dummy plot 
			plot([],[],"b--",markersize=3,label="#nodes parragh")
			plot(k, [mean(res_stidsen[a]) for a in k], "gx", markersize=3) ; plot(k, [mean(res_stidsen[a]) for a in k], "g-", linewidth=1, label="tmoy stidsen")
			plot(k, [mean(res_parragh[a]) for a in k], "bx", markersize=3) ; plot(k, [mean(res_parragh[a]) for a in k], "b-", linewidth=1, label="tmoy parragh")
			plot(k, [mean(res_eps[a])for a in k], "k-", linewidth=1, label="tmoy epsilon")
			legend()
		end
		println("\n####################################\nn = $n")
		println("tmoy eps : $(mean(res_eps[n]))") ; println("tmoy stidsen : $(mean(res_stidsen[n]))") ; println("tmoy parragh : $(mean(res_parragh[n]))")
		println("moy : $(round(mean(nb_nodes_stidsen[n]),1)) noeuds par BC Stidsen") ; println("moy : $(round(mean(nb_nodes_parragh[n]),1)) noeuds par BC Parragh")
	end
	res_stidsen, res_parragh, res_eps, nb_nodes_stidsen, nb_nodes_parragh
end

benchmark(1:15, 10:11, 1, use_nsga=true, global_branch=false, docovercuts=true, lift_covers=true)
# inst = hard_instance2(CplexSolver(CPX_PARAM_SCRIND = 0))
srand(0)
benchmark(1:15, 10:5:40, 5, use_nsga=true, global_branch=false, docovercuts=true, lift_covers=true)
