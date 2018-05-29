using MOBC, vOptGeneric, CPLEX, PyPlot, Suppressor
using Base.Test

function random_instance(range, n)
	m = vModel(solver = CplexSolver(CPX_PARAM_SCRIND = 0))
	p1 = rand(range, n)
	p2 = rand(range, n)
	w = [rand(range, n) for i = 1:3]
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
    n = 25
    p1 = [90 ,100 ,97 ,82 ,84 ,87 ,99 ,87 ,84 ,89 ,90 ,93 ,86 ,86 ,82 ,96 ,83 ,82 ,97 ,85 ,84 ,89 ,82 ,90 ,81]
    p2 = [96 ,81 ,87 ,97 ,80 ,85 ,95 ,100 ,85 ,86 ,98 ,95 ,80 ,98 ,98 ,86 ,80 ,91 ,98 ,87 ,96 ,94 ,98 ,87 ,84]
    w1 = [100 ,95 ,80 ,98 ,98 ,99 ,81 ,83 ,86 ,98 ,86 ,81 ,82 ,91 ,94 ,86 ,92 ,100 ,83 ,87 ,91 ,93 ,92 ,80 ,93]
    w2 = [88 ,91 ,89 ,93 ,100 ,99 ,92 ,98 ,100 ,99 ,99 ,96 ,80 ,94 ,92 ,85 ,99 ,94 ,83 ,96 ,84 ,99 ,94 ,82 ,86]
    w3 = [90 ,94 ,80 ,93 ,92 ,95 ,82 ,82 ,98 ,80 ,84 ,84 ,89 ,94 ,97 ,96 ,82 ,83 ,82 ,93 ,94 ,83 ,97 ,83 ,96]
    b = [1123, 1156, 1111]

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
	figure(1) ; clf()
	figure(2) ; clf()
	res_stidsen=Dict{Int, Vector{Float64}}() ; nb_nodes_stidsen=Dict{Int, Vector{Int}}()
	res_parragh=Dict{Int, Vector{Float64}}() ; nb_nodes_parragh=Dict{Int, Vector{Int}}()
	res_eps=Dict{Int, Vector{Float64}}()
	for n in nrange
		res_stidsen[n]=Float64[] ; res_parragh[n]=Float64[] ; res_eps[n]=Float64[]
		nb_nodes_stidsen[n]=Int[] ; nb_nodes_parragh[n]=Int[]
		for i = 1:runpersize
			m = random_instance(itemrange, n)
			print(m)
			global m_current = m
			valBC, tBC, _ = @timed solve_stidsen(m, Inf ; args...);
			valPAR, tPAR, _ = @timed solve_parragh(m, Inf ; args...);
			status,tEPS,_ = @suppress @timed solve(m, method=:epsilon, round_results=true);
			if !(length(getY_N(m)) == length(unique(valBC.YN)) == length(unique(valPAR.YN)))
				@suppress solve(m, method=:epsilon);
				@show sort(valBC.YN, by = x->x[1])
				@show getY_N(m)
				@show unique(cast_to_int_and_filter(getY_N(m), :Max))
				figure(3)
				clf()
				plot(map(x->x[1], getY_N(m)), map(x->x[2], getY_N(m)), "b>", markersize="5")
				plot(map(x->x[1], valBC.YN), map(x->x[2], valBC.YN), "r<", markersize="5")
				error()
			end

			push!(res_stidsen[n], tBC)
			push!(res_parragh[n], tPAR)
			push!(res_eps[n], tEPS)
			push!(nb_nodes_parragh[n], valPAR.nodes)
			push!(nb_nodes_stidsen[n], valBC.nodes)

			figure(1)
			clf()
			for k in keys(res_stidsen)
				plot(fill(k, length(res_stidsen[k])), res_stidsen[k], "gx", markersize="3")
				plot(fill(k, length(res_stidsen[k])), res_parragh[k], "bx", markersize="3")
				plot(fill(k, length(res_stidsen[k])), res_eps[k], "b.", markersize="3")
			end

			figure(2)
			clf()
			k = sort(collect(keys(res_stidsen)))
			plot(k, [mean(nb_nodes_stidsen[a]) for a in k], "gs", markersize="3")
			plot(k, [mean(nb_nodes_parragh[a]) for a in k], "bs", markersize="3")
			twinx()
			plot([], [], "gs", markersize="3", label="#nodes_stidsen") #just for legend
			plot([], [], "bs", markersize="3", label="#nodes_parragh") #just for legend
			plot(k, [mean(res_stidsen[a]) for a in k], "gx", markersize="3", label="tmoy stidsen")
			plot(k, [mean(res_parragh[a]) for a in k], "bx", markersize="3", label="tmoy parragh")
			plot(k, [mean(res_eps[a])for a in k], "b.", markersize="3", label="tmoy Eps")
			legend()
		end
		println("n = $n")
		println("tmoy eps : $(mean(res_eps[n]))")
		println("tmoy stidsen : $(mean(res_stidsen[n]))")
		println("tmoy parragh : $(mean(res_parragh[n]))")
		println("moy : $(mean(nb_nodes_stidsen[n])) noeuds par BC Stidsen")
		println("moy : $(mean(nb_nodes_parragh[n])) noeuds par BC Parragh")
	end
	res_stidsen, res_parragh, res_eps, nb_nodes_stidsen, nb_nodes_parragh
end

benchmark(50:100, 10:10, 1)
# inst = hard_instance2(CplexSolver(CPX_PARAM_SCRIND = 0))
srand(0)
benchmark(50:100, 5:5:40, 5, use_nsga=true, global_branch=true, docovercuts=false)
# solve_BC(inst, Inf)