using MOBC, vOptGeneric, GLPKMathProgInterface, PyPlot
using Base.Test

m = vModel(solver = GLPKSolverMIP())

p1 = [65,90,90,77,95,84,70,94,66,92,74,97,60,60,65,97,93]
p2, p2cont = [77,94,71,63,96,82,85,75,72,91,99,63,84,87,79,94,90], [74,23,98,45,100,49,42]
w, wcont = [80,87,68,72,66,77,99,85,70,93,98,72,100,89,67,86,91], [45,98,15,65,78,12,29]
c = 900

@variable(m, x[1:length(p1)], Bin)
@variable(m, 0 <= y[1:length(p2cont)] <= 1)
@addobjective(m, Max, dot(x, p1))
@addobjective(m, Max, dot(x, p2) + dot(y, p2cont))
@constraint(m, dot(x, w) + dot(y, wcont)<= c)

YN, XE = MOBC.solve_BC(m, 5000)

solve(m, method=:dichotomy)
voptYN=getY_N(m)

plot(map(x->x[1], voptYN), map(x->x[2], voptYN), "bo", markersize="3")
plot(map(x->x[1], YN), map(x->x[2], YN), "rx", markersize="3")