using vOptGeneric, MOBC, PyPlot, CPLEX

function parseSSCFLP(fname, solver = CplexSolver(CPX_PARAM_SCRIND = 0))
    line_to_vec(f) = parse.(Int, split(chomp(readline(f))))
    file = open(fname)
    nbfac, nbclient = line_to_vec(file)
    I = 1:nbfac
    J = 1:nbclient
    readline(file)
    s, f = Int[], Int[]
    for i = 1:nbfac
        si, fi = line_to_vec(file)
        push!(s, si)
        push!(f, fi)
    end
    readline(file)
    d = line_to_vec(file)
    readline(file)
    c = readdlm(file, Int)
    close(file)
    
    m = vModel(solver=solver)
    @variable(m, x[1:nbfac, 1:nbclient], Bin)
    @variable(m, y[1:nbfac], Bin)
    @addobjective(m, Min, vecdot(c, x))
    @addobjective(m, Min, dot(f, y))
    @constraint(m, [j=J], sum(x[i,j] for i ∈ I) == 1)
    @constraint(m, [i=I], sum(d[j]*x[i,j] for j ∈ J) <= s[i]*y[i])

    m
end

for folder in ["size_5_10", "size_10_20"]#, "size_15_30"]
    cd(joinpath(Pkg.dir("MOBC"), "test", "Instances_SSCFLP", folder)) do
        for file in readdir()
            m = parseSSCFLP(file)
            println(file)
            @time solve_stidsen(m, use_nsga=true, global_branch=false, docovercuts=true, lift_covers=true)
            @time solve_parragh(m, use_nsga=true, global_branch=false, docovercuts=true, lift_covers=true)
            @time solve(m, method=:epsilon, round_results=true, verbose=false, suppress_warnings=true)
        end
    end
end