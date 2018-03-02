using MOBC, vOptGeneric, GLPKMathProgInterface, PyPlot
using Base.Test

m = vModel(solver = GLPKSolverMIP())

nbint = 80
nbcont = 100
range_p = 50:100
range_w = 50:100

p1 = rand(range_p, nbint)
p2, p2cont = rand(range_p, nbint), rand(range_p, nbcont)
w, wcont = rand(range_w, nbint), rand(range_w, nbcont)
c = (sum(w) + sum(wcont)) ÷ 2

@variable(m, x[1:length(p1)], Bin)
@variable(m, 0 <= y[1:length(p2cont)] <= 1)
@addobjective(m, Max, dot(x, p1))
@addobjective(m, Max, dot(x, p2) + dot(y, p2cont))
@constraint(m, dot(x, w) + dot(y, wcont)<= c)

YN, XE = MOBC.solve_BC(m, 5000)

solve(m, method=:dichotomy)
voptYN=getY_N(m)


sleep(2)
clf()
plot(map(x->x[1], voptYN), map(x->x[2], voptYN), "bo", markersize="3")
plot(map(x->x[1], YN), map(x->x[2], YN), "rx", markersize="3")
