struct ConstraintData
    nb_cstr::Int
    indices::Vector{Vector{Int}}
    coeffs::Vector{Vector{Float64}}
    lb::Vector{Float64}
    ub::Vector{Float64}
end

ConstraintData(m) = ConstraintData(
    length(m.linconstr),
    [map(v-> getfield(v, :col), CSTR.terms.vars) for CSTR in m.linconstr], 
    [CSTR.terms.coeffs for CSTR in m.linconstr], 
    [CSTR.lb for CSTR in m.linconstr], 
    [CSTR.ub for CSTR in m.linconstr]
)

isacover(C, a, b) = sum(a[C]) > b
function extend(C, a)
    m = maximum(a[C])
    res = copy(C)
    for i in eachindex(a)
        if a[i] >= m && !(i in res)
            push!(res, i)
        end
    end
    sort!(res)
end
function isviolated(C, x, card)
    # println("##########")
    # @show x[C]
    # @show sum(x[C])
    # @show card
    # println("##########")
    sum(x[C]) > card + 1e-3
end

function separateHeur(x, a, b)
    order = sortperm((1 .- x)./ a)
    C = Int[]
    astar = b
    while !isempty(order)
        p = popfirst!(order)
        a[p] <= astar && push!(C, p)
        isempty(C) && return (Set{Int}(), 0)
        if isacover(C, a, b)
            ECI = extend(C, a)
            if isviolated(ECI, x, length(C)-1)
                # @show C, sortperm((1 .- x)./ a)
                return (ECI, length(C)-1)
            else
                kstar = argmax(a[C])
                astar = a[kstar]
                deleteat!(C, kstar)
            end
        end
    end
    Int[], 0
end

function find_cover_cuts(n, cstrData)
    # @show length(n.m.linconstr)
    cuts = Tuple{Vector{Int}, Int}[]
    for i = 1:cstrData.nb_cstr
        if cstrData.ub[i] > 0 && cstrData.lb[i] == -Inf && cstrData.ub[i] != Inf
            C, card = separateHeur(view(n.m.colVal, cstrData.indices[i]), cstrData.coeffs[i], cstrData.ub[i])
            if !isempty(C)
                @assert isviolated(C, view(n.m.colVal, cstrData.indices[i]), card)
                push!(cuts, (cstrData.indices[i][C], card))
            end
        end
    end
    isempty(cuts) && return false
    # @show cuts
    for c in cuts
        vars, card = c
        # @show vars, card
        @constraint(n.m, sum(JuMP.Variable(n.m, j) for j in vars) <= card)
    end
    # @show n.m.colVal
    
    true
end