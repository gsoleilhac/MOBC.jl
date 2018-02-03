using MOBC, vOptGeneric, GLPKMathProgInterface
using Base.Test

m = vModel(solver = GLPKSolverMIP())

p1 = [65,90,90,77,95,84,70,94,66,92,74,97,60,60,65,97,93,60,69,74]
p2, p2cont = [77,94,71,63,96,82,85,75,72,91,99,63,84,87,79,94,90,60,69,62], [74,23,98,45,100,49,42]
w, wcont = [80,87,68,72,66,77,99,85,70,93,98,72,100,89,67,86,91,79,71,99], [45,98,15,65,78,12,29]
c = 900
size = length(p1)

@variable(m, x[1:size], Bin)
@variable(m, 0 <= y[1:7] <= 1)
@addobjective(m, Max, dot(x, p1))
@addobjective(m, Max, dot(x, p2) + dot(y, p2cont))
@constraint(m, dot(x, w) + dot(y, wcont)<= c)

n = MOBC.solve_BC(m)