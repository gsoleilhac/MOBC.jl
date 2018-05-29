function plot_int_found(LN, lnglobal, z1, z2 ; marker = "b.", sleeptime = 0.01)

	figure(1) ; clf()
	LNPLOT=[]
	for i = 1:length(LN)
		push!(LNPLOT, LN.yn[i])
		i!=length(LN) && push!(LNPLOT, nadirs(LN)[i])
	end
	plot(map(x->x[1], LNPLOT), map(x->x[2], LNPLOT), "r--")
	plot(map(x->x[1], LN.yn), map(x->x[2], LN.yn), "r^")
	plot(map(x->x[1], nadirs(LN)), map(x->x[2], nadirs(LN)), "rs", markersize="4")
	plot([z1], [z2], marker)

	λ = LN.λ
	WS = sum((z1,z2).*λ)
	plot([LN.yn[1][1], (WS-λ[2]*LN.yn[end][2])/λ[1]], [(WS-λ[1]*LN.yn[1][1])/λ[2], LN.yn[end][2]], "k--")

	ax = gca()
    ax[:set_xlim]([-1+min(LN.yn[1][1], z1), 1+max(LN.yn[end][1], z1+1)])
	ax[:set_ylim]([-1+min(LN.yn[end][2]), 1+max(LN.yn[1][2], z2+1)])
	show()

	# figure(2)
	# LNPLOT=[]
	# for i = 1:length(lnglobal)
	# 	push!(LNPLOT, lnglobal.yn[i])
	# 	i!=length(lnglobal) && push!(LNPLOT, nadirs(lnglobal)[i])
	# end
	# plot(map(x->x[1], LNPLOT), map(x->x[2], LNPLOT), "r--")
	# plot(map(x->x[1], lnglobal.yn), map(x->x[2], lnglobal.yn), "r^", markersize="2")
	# plot(map(x->x[1], nadirs(lnglobal)), map(x->x[2], nadirs(lnglobal)), "rs", markersize="2")
	# # plot([0, z1+(lnglobal.λ[2]/lnglobal.λ[1])*z2 ], [z2+(lnglobal.λ[1]/lnglobal.λ[2])*z1, 0], "b--")
	# plot([z1], [z2], marker, markersize="2")
	# show()

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

	ax = gca()
    ax[:set_xlim]([-1+min(LN.yn[1][1], z1), 1+max(LN.yn[end][1], z1)])
	ax[:set_ylim]([-1+min(LN.yn[end][2]), 1+max(LN.yn[1][2])])
	show()
	
	sleep(sleeptime)
end


function plotdualbound(sl::Vector{Segment{T}}) where T
	for s in sl
		if T == Min
			plot([s.p1.x, s.p2.x, s.c.x, s.c.x, s.p1.x, s.p1.x], [s.p1.y, s.p2.y, s.p2.y, s.c.y, s.c.y, s.p1.y], "g-", linewidth=1)
		else
			plot([s.p1.x, s.p2.x, s.p2.x, s.c.x, s.c.x, s.p1.x], [s.p1.y, s.p2.y, s.c.y, s.c.y, s.p1.y, s.p1.y], "g-", linewidth=1)
		end
		plot([s.p1.x, s.p2.x], [s.p1.y, s.p2.y], "gs", markersize=3)
		plot([s.c.x, s.c.x], [s.c.y, s.c.y], "kx", markersize=3)
	end
end