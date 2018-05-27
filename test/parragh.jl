using MOBC, vOptGeneric, CPLEX, vOptGeneric, PyPlot, Suppressor
using Base.Test

function random_instance(range, n)
	m = vModel(solver = CplexSolver(CPX_PARAM_SCRIND = 0))
	p1 = rand(range, n)
	p2 = rand(range, n)
	w = [rand(range, n) for i = 1:3]
	c = 3*sum.(w).÷5

	@variable(m, x[1:length(p1)], Bin)
	@addobjective(m, Min, dot(x, p1))
	@addobjective(m, Min, dot(x, p2))
	@constraint(m, dot(x, w[1]) >= c[1])
	@constraint(m, dot(x, w[2]) >= c[2])
	@constraint(m, dot(x, w[3]) >= c[3])
	m
end

m = random_instance(50:100, 20)

solve(m, method=:dicho)
ub = getY_N(m);
ub = MOBC.toPointList(ub);

solve(m, method=:dicho, relax=true)
lb = getY_N(m);
nadir = (maximum(map(first, ub)), maximum(map(last, ub)));
lb = MOBC.toSegmentList(lb, nadir);

clf()
MOBC.plotdualbound(lb)
plot(first.(ub), last.(ub), "ko")

MOBC.filterLB(lb, ub, true)



ub = MOBC.toPointList([[2,13], [4,10], [6,6], [10,4], [14,1]])
lb = MOBC.toSegmentList([[1,15], [3,11], [5,8], [9,6], [14,5]], (14,13))