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

function isminimal(C, a, b)
    sum(a[C]) - minimum(a[C]) <= b
end

function isastrongcover(C, inds, a, b)
    !isminimal(C, a, b) && return false
    EC = extend(C, a)
    length(EC) == length(a) && return true
    aij1 = maximum(a[C])
    aij2 = maximum(a[setdiff(1:length(a), EC)])
    return sum(a[C]) - aij1 + aij2 <= b
end

function isviolated(C, x, card)
    sum(x[C]) > card + 1e-3
end

function separateHeur(x, a, b)
    order = sortperm((1 .- x)./ a)
    C = Int[]
    astar = b
    while !isempty(order)
        p = popfirst!(order)
        a[p] <= astar && push!(C, p)
        isempty(C) && return (Int[], 0)
        if isacover(C, a, b)
            ECI = extend(C, a)
            if isviolated(ECI, x, length(C)-1)
                return (C, length(C)-1)
            else
                kstar = argmax(a[C])
                astar = a[kstar]
                deleteat!(C, kstar)
            end
        end
    end
    Int[], 0
end

function find_cover_cuts(n, cstrData, lift_covers)
    cuts = Tuple{Vector{Int}, Int}[]
    success = false
    for i = 1:cstrData.nb_cstr
        if cstrData.ub[i] > 0 && cstrData.lb[i] == -Inf && cstrData.ub[i] != Inf
            inds = setdiff(cstrData.indices[i], union(n.f0, n.f1)) #indices of unfixed variables involved in the constraint
            x = view(n.m.colVal, inds) #value of those variables in the current solution
            a = cstrData.coeffs[i][in.(cstrData.indices[i], (inds,))] #coefficients of those variables in the current constraint
            b = cstrData.ub[i] #RHS of the constraint

            
            # Remove the sum of coeffs of variables already fixed to 1 from the RHS
            for j = 1:length(cstrData.indices[i])
                if cstrData.indices[i][j] in n.f1
                    b -= cstrData.coeffs[i][j]
                end
            end
            
            if dot(x, a) > b - 1e-4 #if the constraint is activated
                @assert length(a) == length(x) == length(inds)
                
                C, card = separateHeur(x, a, b) #try to find a cover
                if !isempty(C) #if success, add it to the list of cover cuts
                    # println("\n")
                    # @show inds
                    # @show a, b
                    # success = true

                    # println(C, " | " , inds[C])
                    # println(n.m.colVal[inds[C]],"> $card")
                    @assert isviolated(inds[extend(C, a)], n.m.colVal, card)
                    @assert isempty(intersect(inds[C], union(n.f0, n.f1))) "C : $C, f0 : $(n.f0), f1 : $(n.f1)"

                    # @show isminimal(C, a, b)
                    # @show isastrongcover(C, inds, a, b)

                    if lift_covers && isastrongcover(C, inds, a, b)
                        res = lifting_procedure(C, inds, a, b, card)
                        # println(res)
                        @constraint(n.m, dot([JuMP.Variable(n.m, j) for j in inds], res) <= card)

                    else
                        @constraint(n.m, sum(JuMP.Variable(n.m, j) for j in inds[extend(C, a)]) <= card)
                    end
                end
            end
        end
    end
   success
end

function lifting_procedure(C, inds, a, b, card)
    sorted_ECi = sort(extend(C, a), by = x->a[x], rev=true)
    q = card
    n0 = setdiff(1:length(a), sorted_ECi)
    aj = ones(Int, length(a))
    aj[n0] .= 0

    lb, ub = 0, a[sorted_ECi[1]]
    for h = 2:q
        lb += a[sorted_ECi[h]]
        ub += a[sorted_ECi[h+1]]
        nh = find(x-> lb <= x <= ub, a)
        aj[nh] .= h
    end

    return aj

end