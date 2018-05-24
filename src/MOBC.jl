module MOBC
using JuMP, vOptGeneric, PyPlot, NSGAII, Suppressor, NamedTuples, Nullables, Compat, ProgressMeter

export solve_stidsen, solve_parragh, parseGAP, parseGAP_orlib

include("Node.jl")
include("Plots.jl")
include("NonDomPoints.jl")
include("CoverCuts.jl")
include("ParserGAP.jl")
include("Stidsen.jl")
include("Parragh.jl")

isonlybinary(m, obj) = all(m.colCat[map(x->getfield(x, :col), obj.aff.vars)] .== :Bin)
isbinary(x) = all(y-> y==1. || y==0., x)
dominates(::Type{Min}, x, y) = x[1] < y[1] && x[2] <= y[2] || x[1] <= y[1] && x[2] < y[2]
dominates(::Type{Max}, x, y) = x[1] > y[1] && x[2] >= y[2] || x[1] >= y[1] && x[2] > y[2]
weakly_dominates(::Type{Min}, x, y) = x[1] <= y[1] && x[2] <= y[2]
weakly_dominates(::Type{Max}, x, y) = x[1] >= y[1] && x[2] >= y[2]

evaluate(x, obj) = evaluate(x, obj.aff.constant, obj.aff.coeffs, map(x->getfield(x, :col), obj.aff.vars))
function evaluate(x, cst, coeffs, vars)::Float64
	res = cst
	for i = 1:length(coeffs)
	    @inbounds res += coeffs[i] * x[vars[i]]
	end
	res
end

# function integerbranch(n)
# 	newnode = copy(n, false)
# 	inds0 = find(x->x==0., n.x)
# 	inds1 = find(x->x==1., n.x)
# 	#NO-GOOD CSTR
# 	@constraint(newnode.m, sum(JuMP.Variable(newnode.m, j) for j in inds0) + sum(1. - JuMP.Variable(newnode.m, j) for j in inds1) >= 1.)
# 	return newnode
# end

end # module
