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
	m, p1, p2
end

# srand(14)

m, p1, p2 = random_instance(50:100, 20)

solve(m, method=:dicho)
ub = getY_N(m);
ub = MOBC.toPointList(ub);
nadir = (maximum(first.(ub)), maximum(last.(ub)));

@constraint(m, dot(m[:x], p1) <= nadir[1])
@constraint(m, dot(m[:x], p2) <= nadir[2])


@constraint(m, m[:x][rand(1:7)] == rand(0:1))
@constraint(m, m[:x][rand(8:14)] == rand(0:1))
@constraint(m, m[:x][rand(15:20)] == rand(0:1))

solve(m, method=:dicho, relax=true)
lb = getY_N(m);
lb = MOBC.toSegmentList(lb, nadir, MOBC.Min);

clf() ; MOBC.plotdualbound(lb) ; plot(first.(ub), last.(ub), "ko")

MOBC.filterDualBound(lb, ub, true)



# ub = MOBC.toPointList([[2,13], [4,10], [6,6], [10,4], [14,1]]);
# lb = MOBC.toSegmentList([[1,15], [3,11], [5,8], [9,6], [14,5]], (14,13));