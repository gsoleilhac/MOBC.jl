slices(n) = 1./tan.(linspace(0,π/2,n+2)[2:end-1])

function RCNP(::Type{Min}, m, α, Ƶ1, Ƶ2)

	@show α
	mSUP = copy(m); mINF = copy(m)
	@constraint(mSUP, copy(Ƶ2, mSUP).aff >= α*copy(Ƶ1, mSUP).aff)
	@constraint(mINF, copy(Ƶ2, mINF).aff <= α*copy(Ƶ1, mINF).aff)
	s1 = @suppress solve(mSUP, method=:lexico)
	s2 = @suppress solve(mINF, method=:lexico)

	zA = last(getY_N(mSUP))
	xA = [getvalue(JuMP.Variable(m, i), 2) for i = 1:vm.numCols]

	zB = first(getY_N(mINF))
	xB = [getvalue(JuMP.Variable(m, i), 1) for i = 1:vm.numCols]

	@show zA, zB

	if zA[1] <= zB[1] && zA[2] <= zB[2]
		@constraint(mINF, copy(Ƶ2,mINF).aff <= zA[2] - 1e-4)
		s2 = @suppress solve(mINF, method=:lexico)
		@show s2
		if s2 == :Infeasible
			return Nullable((xA, zA)), Nullable{typeof((xA,zA))}()
		end
		zB = first(getY_N(mINF))
		xB = [getvalue(JuMP.Variable(m, i), 1) for i = 1:vm.numCols]
	elseif zB[1] <= zA[1] && zB[2] <= zA[2]
		@constraint(mSUP, copy(Ƶ1,mSUP).aff <= zB[1] - 1e-4)
		s1 = @suppress solve(mSUP, method=:lexico)
		@show s1
		if s1 == :Infeasible
			return Nullable{typeof((xB,zB))}(), Nullable((xB, zB))
		end
		zA = last(getY_N(mSUP))
		xA = [getvalue(JuMP.Variable(m, i), 2) for i = 1:vm.numCols]
	end
	println("###")
	return Nullable((xA,zA)), Nullable((xB,zB))
end


function RCNP(::Type{Max}, m, α, Ƶ1, Ƶ2)

	@show α
	mSUP = copy(m); mINF = copy(m)
	@constraint(mSUP, copy(Ƶ2, mSUP).aff >= α*copy(Ƶ1, mSUP).aff)
	@constraint(mINF, copy(Ƶ2, mINF).aff <= α*copy(Ƶ1, mINF).aff)
	s1 = @suppress solve(mSUP, method=:lexico)
	s2 = @suppress solve(mINF, method=:lexico)

	zA = first(getY_N(mSUP))
	xA = [getvalue(JuMP.Variable(m, i), 1) for i = 1:m.numCols]

	zB = last(getY_N(mINF))
	xB = [getvalue(JuMP.Variable(m, i), 1) for i = 1:m.numCols]

	@show zA, zB

	if zA[1] >= zB[1] && zA[2] >= zB[2]
		@constraint(mINF, copy(Ƶ2,mINF).aff >= zA[2] + 1e-4)
		s2 = @suppress solve(mINF, method=:lexico)
		@show s2
		if s2 == :Infeasible
			return Nullable((xA, zA)), Nullable{typeof((xA,zA))}()
		end
		zB = last(getY_N(mINF))
		xB = [getvalue(JuMP.Variable(m, i), 1) for i = 1:m.numCols]
	elseif zB[1] >= zA[1] && zB[2] >= zA[2]
		@constraint(mSUP, copy(Ƶ1,mSUP).aff >= zB[1] + 1e-4)
		s1 = @suppress solve(mSUP, method=:lexico)
		@show s1
		if s1 == :Infeasible
			return Nullable{typeof((xB,zB))}(), Nullable((xB, zB))
		end
		zA = first(getY_N(mSUP))
		xA = [getvalue(JuMP.Variable(m, i), 1) for i = 1:m.numCols]
	end
	@show zA, zB
	println("###")
	return Nullable((xA,zA)), Nullable((xB,zB))
end


function select_slices(YN, k)
	K = length(YN)
	@assert k <= K÷2
	G = OffsetArray(Vector{Int}, 0:k)
	G[0] = [0]
	G[k] = [K-1]
	for i = 1:k-1
		G[i] = collect(2i-1:K-2k+2i-1)
	end
	@show G
	
	d = Dict{Tuple{Int,Int}, Int}()
	for col = 0:k-1
		for i in G[col] .+ 1
			for j in G[col+1] .+ 1
				if col == 0
					d[i,j] = abs((YN[i][1]-YN[j][1]) * (YN[i][2]-YN[j][2]))
				else	
					if j >= i+2
						d[i,j] = abs((YN[i+1][1]-YN[j][1]) * (YN[i+1][2]-YN[j][2]))
					end
				end
			end
		end
	end
	
	m = Model(solver = CplexSolver(CPX_PARAM_SCRIND = 0))
	@variable(m, x[keys(d)], Bin)
	@variable(m, wmax, Int) 

	@objective(m, Min, wmax)

	edges_1_x = Iterators.filter(x->first(x)==1, keys(d))
	edges_x_K = Iterators.filter(x->last(x)==K, keys(d))

	@constraint(m, sum(x[e] for e in edges_1_x) == 1)
	@constraint(m, sum(x[e] for e in edges_x_K) == 1)
	@constraint(m, [i = keys(d)], wmax >= x[i]*d[i])
	@constraint(m, sum(x) == k)
	for i = 2:K-1
		in_edges = Iterators.filter(x->last(x)==i, keys(d))
		out_edges = Iterators.filter(x->first(x)==i, keys(d))
		@constraint(m, sum(x[j] for j in in_edges) == sum(x[l] for l in out_edges))
	end
	
	@show solve(m)
	@show getobjectivevalue(m)

	x_res = getvalue(x)

	res=Int[]
	i = 1
	j = 2
	while j != K
		if i==1 && j==2 && x_res[(i,j)] == 1.
			push!(res, 2)
			i=2
		elseif length(res) == k-1
			@assert x_res[(i,K)] == 1 res
			push!(res, K)
			j = K
		else
			j = findfirst(x->x_res[(i,x)]==1, (i+2):K) + (i+1)
			push!(res, j)
			i = j
		end
	end
	res


end