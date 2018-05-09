using MOBC, vOptGeneric, CPLEX, PyPlot, Suppressor
using Base.Test

function random_instance(range, n, nbcstr=10, density=0.2)
	m = vModel(solver = CplexSolver(CPX_PARAM_SCRIND = 0))
	p1 = rand(range, n)
	p2 = rand(range, n)

	@variable(m, x[1:n], Bin)
	@addobjective(m, Max, dot(x, p1))
	@addobjective(m, Max, dot(x, p2))
    for i = 1:nbcstr
        mask = rand(n) .< density
        while count(mask) < 2
            mask[rand(1:n)] = true
        end
        @constraint(m, sum(x[mask]) <= 1)
    end
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

function benchmark(itemrange, nrange, runpersize)
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
			valBC, tBC, _ = @timed solve_BC(m, Inf);
			@suppress valEPS,tEPS,_ = @timed solve(m, method=:epsilon);
			if length(unique(cast_to_int_and_filter(getY_N(m), :Max))) != length(unique(valBC.YN))
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
benchmark(10:100, 5:3:35, 25)
# solve_BC(inst, Inf)