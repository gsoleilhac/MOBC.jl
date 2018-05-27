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

#tres long sans cover cuts
function hard_instance(solver)
    n = 15
    p1 = [95 ,85 ,80 ,82 ,100 ,95 ,81 ,94 ,82 ,84, 88 ,99 ,89 ,91 ,98]
    p2 = [100 ,88 ,85 ,91 ,94 ,99 ,97 ,94 ,97 ,87, 96 ,92 ,91 ,99 ,80]
    w1 = [98 ,95 ,91 ,94 ,92 ,95 ,88 ,94 ,86 ,84, 86 ,85 ,82 ,84 ,92]
    w2 = [84 ,100 ,88 ,97 ,97 ,85 ,100 ,98 ,98 ,82, 81 ,80 ,81 ,88 ,89]
    w3 = [84 ,83 ,94 ,84 ,83 ,85 ,82 ,99 ,92 ,97, 91 ,92 ,96 ,80 ,100]
    b = [673, 674, 671]

    m = vModel(solver = solver)
    @variable(m, x[1:n], Bin)
    @addobjective(m, Max, dot(x,p1))
    @addobjective(m, Max, dot(x,p2))
    @constraint(m, dot(x, w1) <= b[1])
    @constraint(m, dot(x, w2) <= b[2])
    @constraint(m, dot(x, w3) <= b[3])

    m
end

#tres long avec cover cuts
function hard_instance2(solver)
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

function cast_to_int_and_filter(yn::Vector{Vector{Float64}}, sense=:Max)
	yn = [round.(Int, x) for x in yn]
	inds = Int[]
	dominates = sense==:Max ? (a,b) -> a[1] >= b[1] && a[2] >= b[2] : (a,b) -> a[1] <= b[1] && a[2] <= b[2]
	for i = 1:length(yn)-1
		if dominates(yn[i], yn[i+1])
			push!(inds, i+1)
		elseif dominates(yn[i+1], yn[i])
			push!(inds, i)
		end
	end
	deleteat!(yn, inds)
	yn
end

function benchmark(itemrange, nrange, runpersize ; args...)
	figure(1) ; clf()
	figure(2) ; clf()
	res_bc=Dict{Int, Vector{Float64}}()
	res_eps=Dict{Int, Vector{Float64}}()
	nb_nodes=Dict{Int, Vector{Int}}()
	for n in nrange
		res_bc[n]=Float64[]
		res_eps[n]=Float64[]
		nb_nodes[n]=Int[]
		for i = 1:runpersize
			m = random_instance(itemrange, n)
			print(m)
			global m_problem = m
			valBC, tBC, _ = @timed solve_stidsen(m, Inf ; args...);
			status,tEPS,_ = @suppress @timed solve(m, method=:epsilon, round_results=true);
			if length(getY_N(m)) != length(unique(valBC.YN))
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

			push!(res_bc[n], tBC)
			push!(res_eps[n], tEPS)
			push!(nb_nodes[n], valBC.nodes)

			figure(1)
			clf()
			for k in keys(res_bc)
				plot(fill(k, length(res_bc[k])), res_bc[k], "rx", markersize="3")
				plot(fill(k, length(res_bc[k])), res_eps[k], "b.", markersize="3")
			end

			figure(2)
			clf()
			k = sort(collect(keys(res_bc)))
			plot(k, [mean(nb_nodes[a]) for a in k], "gs", markersize="3")
			twinx()
			plot([], [], "gs", markersize="3", label="#nodes") #just for legend
			plot(k, [mean(res_bc[a]) for a in k], "rx", markersize="3", label="tmoy BC")
			plot(k, [mean(res_eps[a])for a in k], "b.", markersize="3", label="tmoy Eps")
			legend()
		end
		println("n = $n")
		println("tmoy eps : $(mean(res_eps[n]))")
		println("tmoy bc : $(mean(res_bc[n]))")
		println("moy : $(mean(nb_nodes[n])) noeuds par BC")
	end
	res_bc, res_eps, nb_nodes
end

benchmark(50:100, 10:10, 1)
# inst = hard_instance2(CplexSolver(CPX_PARAM_SCRIND = 0))
srand(0)
benchmark(10:100, 5:5:40, 10, use_nsga=false)
# solve_BC(inst, Inf)