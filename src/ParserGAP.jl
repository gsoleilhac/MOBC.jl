using CPLEX

struct MO_GAP_Inst
    m::Model
    a::Matrix{Int}
    b::Vector{Int}
    c1::Matrix{Int}
    c2::Matrix{Int}
end


function parseGAP(fname, solver=CplexSolver())
    f = open(fname)
    n,m = parse.(Int, split(readline(f)))

    a = zeros(Int, n, m)
    c = zeros(Int, n, m)

    for i = 1:n
        cij = Int[]
        while length(cij) < m
            append!(cij, parse.(Int, split(readline(f))))
        end
        @assert length(cij) == m
        c[i, :] .= cij
    end
    
    for i = 1:n
        aij = Int[]
        while length(aij) < m
            append!(aij, parse.(Int, split(readline(f))))
        end
        @assert length(aij) == m
        a[i, :] .= aij
    end

    b = Int[]
    while length(b) < n
        append!(b, parse.(Int, split(readline(f))))
    end

    c2 = reshape(shuffle(MersenneTwister(0), vec(c)), n, m)

    #m = m÷4

    mod = vModel(solver=solver)
    @variable(mod, x[1:n, 1:m], Bin)
    @addobjective(mod, Min, sum(c[i,j]*x[i,j] for i = 1:n, j = 1:m))
    @addobjective(mod, Min, sum(c2[i,j]*x[i,j] for i = 1:n, j = 1:m))
    @constraint(mod, [i=1:n], sum(a[i,j]*x[i,j] for j = 1:m) <= b[i])
    @constraint(mod, [j=1:m], sum(x[i,j] for i = 1:n) == 1)

    mod,a,b,c,c2
end

function parseGAP_orlib(fname, solver=CplexSolver(CPX_PARAM_SCRIND = 0))
    f = open(fname)
    nbinst = parse(Int, chomp(readline(f)))
    res = MO_GAP_Inst[]

    for _ = 1:nbinst

        n,m = parse.(Int, split(readline(f)))
        a = zeros(Int, n, m)
        c = zeros(Int, n, m)

        for i = 1:n
            cij = Int[]
            while length(cij) < m
                append!(cij, parse.(Int, split(readline(f))))
            end
            @assert length(cij) == m
            c[i, :] .= cij
        end

        for i = 1:n
            aij = Int[]
            while length(aij) < m
                append!(aij, parse.(Int, split(readline(f))))
            end
            @assert length(aij) == m
            a[i, :] .= aij
        end

        b = Int[]
        while length(b) < n
            append!(b, parse.(Int, split(readline(f))))
        end

        c2 = reshape(shuffle(MersenneTwister(0), vec(c)), n, m)
        #m = m÷4

        mod = vModel(solver=solver)
        @variable(mod, x[1:n, 1:m], Bin)
        @addobjective(mod, Min, sum(c[i,j]*x[i,j] for i = 1:n, j = 1:m))
        @addobjective(mod, Min, sum(c2[i,j]*x[i,j] for i = 1:n, j = 1:m))
        @constraint(mod, [i=1:n], sum(a[i,j]*x[i,j] for j = 1:m) <= b[i])
        @constraint(mod, [j=1:m], sum(x[i,j] for i = 1:n) == 1)

        push!(res, MO_GAP_Inst(mod,a,b,c,c2))
    end
    close(f)
    res
end
