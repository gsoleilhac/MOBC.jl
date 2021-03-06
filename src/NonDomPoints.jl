abstract type Sense end
struct Max<:Sense end
struct Min<:Sense end

mutable struct NonDomPoints{T}
	xe::Vector{Vector{Int}}
	yn::Vector{Tuple{Int,Int}}
	λ::Tuple{Int,Int}
	nadirs::Vector{Tuple{Int,Int}}
	worst_nadir::Tuple{Int,Int}
	function NonDomPoints(t::Type{T}, xe, yn) where T<:Sense
		@assert issorted(yn, by = y-> y[1])
		λ = round(Int, yn[1][2] - yn[end][2]), round(Int, yn[end][1] - yn[1][1])
		λ = λ[1]÷gcd(λ[1], λ[2]), λ[2]÷gcd(λ[1], λ[2])
		res = new{T}([round.(Int, x) for x in xe], [round.(Int, x) for x in yn], λ, [])
		update_ln!(res)
		res
	end
end

function update_ln!(p::NonDomPoints{Max})
	p.nadirs = [(p.yn[i][1], p.yn[i+1][2]) for i = 1:length(p.yn)-1]
	if isempty(p.nadirs)
		p.worst_nadir = p.yn[1].+1
		return p
	end
	p.worst_nadir = p.nadirs[1].+1
	for n in p.nadirs
		if sum((n.+1).*p.λ) < sum(p.worst_nadir.*p.λ)
			p.worst_nadir = n.+1
		end
	end
	p
end

function update_ln!(p::NonDomPoints{Min})
	p.nadirs = [(p.yn[i+1][1], p.yn[i][2]) for i = 1:length(p.yn)-1]
	if isempty(p.nadirs)
		p.worst_nadir = p.yn[1].-1
		return p
	end
	p.worst_nadir = p.nadirs[1].-1
	for n in p.nadirs
		if sum((n.-1).*p.λ) > sum(p.worst_nadir.*p.λ)
			p.worst_nadir = n.-1
		end
	end
	p
end

nadirs(p::NonDomPoints) = p.nadirs
worst_nadir(p::NonDomPoints) = p.worst_nadir
isfathomable(z, p::NonDomPoints{Min}) = z > sum(worst_nadir(p).*p.λ)
isfathomable(z, p::NonDomPoints{Max}) = z < sum(worst_nadir(p).*p.λ)
local_nadir(p::NonDomPoints{Min}) = (p.yn[end][1], p.yn[1][2])
local_nadir(p::NonDomPoints{Max}) = (p.yn[1][1], p.yn[end][2])

function isweaklydominated(z1, z2, p::NonDomPoints{Min})
	for (v1, v2) in p.yn
		if v1 <= z1 && v2 <= z2
			return true
		end
	end
	false
end

function isweaklydominated(z1, z2, p::NonDomPoints{Max})
	for (v1, v2) in p.yn
		if v1 >= z1 && v2 >= z2
			return true
		end
	end
	false
end

Base.length(p::NonDomPoints) = length(p.yn)

function Base.push!(p::NonDomPoints{S}, x, y) where S<:Sense
	#@assert !isweaklydominated(y..., p) "y : $y \n p : $p"
	if isweaklydominated(y..., p)
		return p
	end
	ind = searchsortedfirst(p.yn, y, by = y -> first(y))
	insert!(p.xe, ind, x)
	insert!(p.yn, ind, y)
	inds_delete = find(x->dominates(S, y, x), p.yn)
	
	deleteat!(p.xe, inds_delete)
	deleteat!(p.yn, inds_delete)
	update_ln!(p)
end