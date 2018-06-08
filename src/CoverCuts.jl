mutable struct ConstraintData
    nb_cstr::Int
    indices::Vector{Vector{Int}}
    coeffs::Vector{Vector{Float64}}
    lb::Vector{Float64}
    ub::Vector{Float64}
end

function ConstraintData(m)
    res = ConstraintData(
            length(m.linconstr),
            [map(v-> getfield(v, :col), CSTR.terms.vars) for CSTR in m.linconstr], 
            [CSTR.terms.coeffs for CSTR in m.linconstr], 
            [CSTR.lb for CSTR in m.linconstr], 
            [CSTR.ub for CSTR in m.linconstr])

    to_discard = find(x-> x<=0 || x==Inf, res.ub)
    
    if !isempty(to_discard)
        res.nb_cstr = res.nb_cstr - length(to_discard)
        deleteat!(res.indices, to_discard)
        deleteat!(res.coeffs, to_discard)
        deleteat!(res.lb, to_discard)
        deleteat!(res.ub, to_discard)
    end

    res
end


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

function find_cover_cuts(n::Node, cstrData, lift_covers)
    success = false
    for i = 1:cstrData.nb_cstr
        inds = cstrData.indices[i] #indices of variables involved in the constraint
        x = view(n.m.colVal, inds) #value of the variables in the current solution
        a = cstrData.coeffs[i] #coefficients of the variables in the current constraint
        b = cstrData.ub[i] #RHS of the constraint
        
        if dot(x, a) > b - 1e-4 #if the constraint is activated                
            C, card = separateHeur(x, a, b) #try to find a cover
            if !isempty(C) #if success, add it to the list of cover cuts
                success = true
                if lift_covers && isastrongcover(C, inds, a, b)
                    res = lifting_procedure(C, inds, a, b, card)
                    @constraint(n.m, dot([JuMP.Variable(n.m, j) for j in inds], res) <= card)
                else
                    @constraint(n.m, sum(JuMP.Variable(n.m, j) for j in inds[extend(C, a)]) <= card)
                end
            end
        end
    end
   success
end

function find_cover_cuts(n::NodeParragh, cstrData, lift_covers)
    cuts = Set{Tuple{Vector{Int}, Int}}()
    lifted_cuts = Set{Tuple{Vector{Int}, Vector{Int}, Int}}()
    success = false
    for i = 1:cstrData.nb_cstr
        inds = cstrData.indices[i] #indices of variables involved in the constraint
        a = cstrData.coeffs[i] #coefficients of thothese variables in the current constraint
        b = cstrData.ub[i] #RHS of the constraint

        for ind_x in 1:length(n.x) #for each solution of the convex linear relaxation
            x =  view(n.x[ind_x], inds) #value of the variables in the solution
            
            if dot(x, a) > b - 1e-4 #if the constraint is activated                
                C, card = separateHeur(x, a, b) #try to find a cover
                if !isempty(C) #if success, add it to the list of cover cuts
                    success = true
                    if lift_covers && isastrongcover(C, inds, a, b)
                        res = lifting_procedure(C, inds, a, b, card)
                        @show maximum(res)
                        push!(lifted_cuts, (inds, res, card))
                        # @constraint(n.m, dot([JuMP.Variable(n.m, j) for j in inds], res) <= card)
                    else
                        push!(cuts, (inds[extend(C, a)], card))
                        # @constraint(n.m, sum(JuMP.Variable(n.m, j) for j in inds[extend(C, a)]) <= card)
                    end
                end
            end
        end
    end
    if success
        for (cut, card) in cuts
            @constraint(n.m, sum(JuMP.Variable(n.m, k) for k in cut) <= card)
        end
        for (inds, cut, card) in lifted_cuts
            @constraint(n.m, dot([JuMP.Variable(n.m, k) for k in inds], cut) <= card)
        end
    end
   success
end

function lifting_procedure(C, inds, a, b, card)
    sorted_Ci = sort(C, by = x->a[x], rev=true)
    sorted_ECi = sort(extend(C, a), by = x->a[x], rev=true)
    q = card
    n0 = setdiff(1:length(a), sorted_ECi)
    aj = ones(Int, length(a))
    aj[n0] .= 0

    lb = a[sorted_Ci[1]]
    ub = lb + a[sorted_Ci[2]]
    for h = 2:q
        lb += a[sorted_Ci[h]]
        ub += a[sorted_Ci[h+1]]
        maximum(a) < lb && break
        nh = find(x-> lb <= x < ub, a)
        aj[nh] .= h
    end

    return aj

end