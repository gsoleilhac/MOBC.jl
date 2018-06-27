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
