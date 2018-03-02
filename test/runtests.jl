using MOBC, vOptGeneric, GLPKMathProgInterface, PyPlot
using Base.Test

m = vModel(solver = GLPKSolverMIP())

p1 = [65,90,90,77,95,84,70,94,66,92,74,97,60,60,65,97,93]
p2 = [77,94,71,63,96,82,85,75,72,91,99,63,84,87,79,94,90]
w = [80,87,68,72,66,77,99,85,70,93,98,72,100,89,67,86,91]
c = 900

@variable(m, x[1:length(p1)], Bin)
@addobjective(m, Max, dot(x, p1))
@addobjective(m, Max, dot(x, p2))
@constraint(m, dot(x, w) <= c)

YN, XE = solve_BC(m, 50000)
solve(m, method=:epsilon)
voptYN=getY_N(m)
clf()
plot(map(x->x[1], voptYN), map(x->x[2], voptYN), "b>", markersize="6")
plot(map(x->x[1], YN), map(x->x[2], YN), "r<", markersize="6")

#Random instance

function random_instance(range, n)
	m = vModel(solver = GLPKSolverMIP())
	p1 = rand(range, n)
	p2 = rand(range, n)
	w = rand(range, n)
	c = sum(w)÷2

	@variable(m, x[1:length(p1)], Bin)
	@addobjective(m, Max, dot(x, p1))
	@addobjective(m, Max, dot(x, p2))
	@constraint(m, dot(x, w) <= c)
	m
end


colors = ["b", "g", "r", "k", "y"]
range = 50:100
figure(1) ; clf()
for n = 10:10:100
	for i = 1:5
		m = random_instance(range, n)
		valBC, tBC, _ = @timed solve_BC(m, 5000000);
		valEPS,tEPS,_ = @timed solve(m, method=:epsilon);
		# @show valBC.YN |> unique |> length
		# @show getY_N(m) |> unique |> length
		@show length(unique(getY_N(m))) == length(unique(valBC.YN))
		plot([n], [tBC], "$(colors[i])x", markersize="4")
		plot([n], [tEPS], "$(colors[i])s", markersize="3")

	end
end

