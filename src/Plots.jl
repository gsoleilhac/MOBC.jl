function plot_int_found(LN, lnglobal, z1, z2 ; sleeptime = 0.01)

	# @show LN.yn
	# @show lnglobal.yn

	figure(1)
	clf()
	LNPLOT=[]
	for i = 1:length(LN)
		push!(LNPLOT, LN.yn[i])
		i!=length(LN) && push!(LNPLOT, nadirs(LN)[i])
	end
	plot(map(x->x[1], LNPLOT), map(x->x[2], LNPLOT), "r--")
	plot(map(x->x[1], LN.yn), map(x->x[2], LN.yn), "r^")
	plot(map(x->x[1], nadirs(LN)), map(x->x[2], nadirs(LN)), "rs", markersize="4")
	# plot([0, z1+(LN.λ[2]/LN.λ[1])*z2 ], [z2+(LN.λ[1]/LN.λ[2])*z1, 0], "b--")
	plot([z1], [z2], "bo")

	λ = LN.λ
	WS = sum((z1,z2).*λ)
	plot([LN.yn[1][1], (WS-λ[2]*LN.yn[end][2])/λ[1]], [(WS-λ[1]*LN.yn[1][1])/λ[2], LN.yn[end][2]], "k--")



	figure(2)
	
	LNPLOT=[]
	for i = 1:length(lnglobal)
		push!(LNPLOT, lnglobal.yn[i])
		i!=length(lnglobal) && push!(LNPLOT, nadirs(lnglobal)[i])
	end
	plot(map(x->x[1], LNPLOT), map(x->x[2], LNPLOT), "r--")
	plot(map(x->x[1], lnglobal.yn), map(x->x[2], lnglobal.yn), "r^")
	plot(map(x->x[1], nadirs(lnglobal)), map(x->x[2], nadirs(lnglobal)), "rs", markersize="4")
	# plot([0, z1+(lnglobal.λ[2]/lnglobal.λ[1])*z2 ], [z2+(lnglobal.λ[1]/lnglobal.λ[2])*z1, 0], "b--")
	plot([z1], [z2], "bo")

	sleep(sleeptime)
end

function plot_pareto_branch(LN, z1, z2, bound1, bound2 ; sleeptime = 0.01)


	figure(1)
	clf()


	LNPLOT=[]
	for i = 1:length(LN)
		push!(LNPLOT, LN.yn[i])
		i!=length(LN) && push!(LNPLOT, nadirs(LN)[i])
	end
	plot(map(x->x[1], LNPLOT), map(x->x[2], LNPLOT), "r--")
	plot(map(x->x[1], LN.yn), map(x->x[2], LN.yn), "r^")
	plot(map(x->x[1], nadirs(LN)), map(x->x[2], nadirs(LN)), "rs", markersize="4")
	# plot([0, z1+(LN.λ[2]/LN.λ[1])*z2 ], [z2+(LN.λ[1]/LN.λ[2])*z1, 0], "b--")
	!isnan(bound1) && plot([bound1, bound1], [LN.yn[1][2], LN.yn[end][2]])
	!isnan(bound2) && plot([LN.yn[1][1], LN.yn[end][1]], [bound2, bound2])
	plot([z1], [z2], "bo")

	λ = LN.λ
	WS = sum((z1,z2).*λ)
	plot([LN.yn[1][1], (WS-λ[2]*LN.yn[end][2])/λ[1]], [(WS-λ[1]*LN.yn[1][1])/λ[2], LN.yn[end][2]], "k--")
	sleep(sleeptime)

end