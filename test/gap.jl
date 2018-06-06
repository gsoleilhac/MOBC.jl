using MOBC, vOptGeneric, CPLEX, PyPlot, Suppressor
using Base.Test

function benchmark( ; fig1 = 1, fig2 = 2, args...)
	res_stidsen=Dict{Int, Vector{Float64}}() ; nb_nodes_stidsen=Dict{Int, Vector{Int}}()
	res_parragh=Dict{Int, Vector{Float64}}() ; nb_nodes_parragh=Dict{Int, Vector{Int}}()
	res_eps=Dict{Int, Vector{Float64}}()
	for filename in ("gap1.txt", "gap2.txt")
		filepath = joinpath(Pkg.dir("MOBC"), "test", "Instances_gap", "gap_orlib", filename)
		instances = parseGAP_orlib(filepath)
		n = instances[1].m.numCols
		res_stidsen[n]=Float64[] ; res_parragh[n]=Float64[] ; res_eps[n]=Float64[]
		nb_nodes_stidsen[n]=Int[] ; nb_nodes_parragh[n]=Int[]
		for inst in instances
			m = inst.m
			global m_current = m
			valBC, tBC, _ = @timed solve_stidsen(m, Inf ; args...);
			sleep(0.5)
			valPAR, tPAR, _ = @timed solve_parragh(m, Inf ; args...);
			sleep(0.5)
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

			figure(fig1) ; clf()
			xlabel("n") ; ylabel("t(s)") ; title("time per instance")
			plot([],[],"gx",markersize=3,label="stidsen") ; plot([],[],"bx",markersize=3,label="parragh") ; plot([],[],"k.",markersize=3,label="epsilon")
			for k in keys(res_stidsen)
				plot(fill(k, length(res_stidsen[k])), res_stidsen[k], "gx", markersize="3")
				plot(fill(k, length(res_stidsen[k])), res_parragh[k], "bx", markersize="3")
				plot(fill(k, length(res_stidsen[k])), res_eps[k], "k.", markersize="3")
			end
			legend()

			figure(fig2) ; clf()
			xlabel("n") ; ylabel("#nodes") ; title("average time and #nodes")
			k = sort(collect(keys(res_stidsen)))
			plot(k, [mean(nb_nodes_stidsen[a]) for a in k], "gs", markersize="3") ; plot(k, [mean(nb_nodes_stidsen[a]) for a in k], "g--", linewidth=1)
			plot(k, [mean(nb_nodes_parragh[a]) for a in k], "bs", markersize="3") ; plot(k, [mean(nb_nodes_parragh[a]) for a in k], "b--", linewidth=1)
			twinx() ; ylabel("t(s)")
			plot([],[],"gs",markersize="3",label="#nodes stidsen") #labels dont work after twinx(), had to put a dummy plot 
			plot([],[],"bs",markersize="3",label="#nodes parragh")
			plot(k, [mean(res_stidsen[a]) for a in k], "gx", markersize="3", label="tmoy stidsen") ; plot(k, [mean(res_stidsen[a]) for a in k], "g-", linewidth=1)
			plot(k, [mean(res_parragh[a]) for a in k], "bx", markersize="3", label="tmoy parragh") ; plot(k, [mean(res_parragh[a]) for a in k], "b-", linewidth=1)
			plot(k, [mean(res_eps[a])for a in k], "b.", markersize="3", label="tmoy epsilon")
			legend()
		end
		println("\n####################################\nn = $n")
		println("tmoy eps : $(mean(res_eps[n]))") ; println("tmoy stidsen : $(mean(res_stidsen[n]))") ; println("tmoy parragh : $(mean(res_parragh[n]))")
		println("moy : $(round(mean(nb_nodes_stidsen[n]),1)) noeuds par BC Stidsen") ; println("moy : $(round(mean(nb_nodes_parragh[n]),1)) noeuds par BC Parragh")
	end
	res_stidsen, res_parragh, res_eps, nb_nodes_stidsen, nb_nodes_parragh
end

srand(0)
benchmark(use_nsga=true, global_branch=false, docovercuts=true, lift_covers=true)
