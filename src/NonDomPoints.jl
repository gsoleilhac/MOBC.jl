
mutable struct NonDomPoints
	xe::Vector{Vector{Float64}}
	yn::Vector{Tuple{Float64,Float64}}
	nadirs::Vector{Tuple{Float64,Float64}}
	worst_nadir::Tuple{Float64,Float64}
end

function NonDomPoints(xe, yn)
	@assert issorted(yn, by = x-> x[1])
	res = NonDomPoints(xe, Tuple.(yn), [], (1.,1.))
	update_ln!(res)
end

function update_ln!(p)
	p.nadirs = [(p.yn[i+1][1], p.yn[i][2]) for i = 1:length(p.yn)-1]
	p.worst_nadir = p.nadirs[1]
	for n in p.nadirs
		if sum(n) > sum(p.worst_nadir)
			p.worst_nadir = n
		end
	end
	p
end

nadirs(p) = p.nadirs
worst_nadir(p) = p.worst_nadir
isfathomable(x, p::NonDomPoints) = x >= sum(p.worst_nadir)

function isdominated(z1, z2, p)
	for (v1, v2) in p.yn
		if v1 < z1 && v2 <= z2 || v1 <= z1 && v2 < z2 
			return true
		end
	end
	false
end

Base.length(p::NonDomPoints) = length(p.yn)

function Base.push!(p::NonDomPoints, x, y)
	@assert !isdominated(y..., p)
	ind = searchsortedfirst(p.yn, y, by = x -> first(x))
	insert!(p.xe, ind, x)
	insert!(p.yn, ind, y)
	inds_delete = find(x->dominates(y, x), p.yn)
	deleteat!(p.xe, inds_delete)
	deleteat!(p.yn, inds_delete)
	update_ln!(p)
end