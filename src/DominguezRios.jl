struct Box
	l::Vector{Float64}
    u::Vector{Float64}
    dim::Int
	priority::Float64 
	function Box(l::Vector{Float64}, u::Vector{Float64}, priority::Float64 = 0.0)
        @assert length(l) == length(u) "Dimension mismatch"
		new(l, u, length(l), priority)
	end
end

function volume(b::Box)
    return prod(b.u - b.l)
end

function volume(b::Box, l::Vector{Float64})
    return prod(b.u - l)
end

function Base.isempty(b::Box)
	return volume(b) == 0
end

function scaled_priority(l::Vector{Float64}, u::Vector{Float64}, zᵢ::Vector{Float64}, zₙ::Vector{Float64})
	return prod((u - l) ./ (zₙ - zᵢ))
end

function reduced_scaled_priority(l::Vector{Float64}, u::Vector{Float64}, i::Int, z::Vector{Float64}, zᵢ::Vector{Float64}, zₙ::Vector{Float64})
	if i ≠ length(z)
		return scaled_priority(l, u, zᵢ, zₙ)
	else
		return scaled_priority(l, u, zᵢ, zₙ) - scaled_priority(l, z, zᵢ, zₙ)
	end
end

function p_split(u::Vector{Float64}, z::Vector{Float64})
	if length(u) != length(z)
		error("Dimension mismatch!")
	end
	return [[i == j ? z[i] : u[i] for i in 1:length(u)] for j in 1:length(u)]
end

function p_partition(B::Box, z::Vector{Float64}, zᵢ::Vector{Float64}, zₙ::Vector{Float64})
	ẑ = max.(z, B.l)
    boxes = Box[]
    for i in 1:length(z)
        l_new = vcat(B.l[1:i], ẑ[i+1:end])
        u_new = vcat(B.u[1:i-1], ẑ[i], B.u[i+1:end])
        priority = reduced_scaled_priority(l_new, u_new, i, z, zᵢ, zₙ)
        push!(boxes, Box(l_new, u_new, priority))
    end
    return boxes
end

function compute_r(zᵢ::Vector{Float64}, zᵤ::Vector{Float64})
	return maximum(zᵤ - zᵢ)
end

function compute_ϵ(p::Int, r::Float64)
	return 1 / (2*p*(r-1))
end

function compute_w(z::Vector{Float64}, zᵢ::Vector{Float64})
	return 1 ./ max.(1, z - zᵢ)
end

function holtzman(model::Model, f::Vector{AffExpr}, x::Array{VariableRef}, problem::AbstractMOProblem, w::Vector{Float64}, zᵢ::Vector{Float64}, ϵ::Float64)
	@variable(model, t_max)
	t_constraints = @constraint(model, [j = 1:problem.p], t_max ≥ w[j] * (f[j] - zᵢ[j]))	
	@objective(model, Min, t_max + ϵ * sum(w .* (f - zᵢ)))
	optimize!(model)
	status = termination_status(model)
	obj = objective_value(model)
	f = value.(f)
	x = value.(x)
    delete(model, t_max)
    unregister(model, :t_max)
    delete.(model, t_constraints)
	return status, obj, f, x
end

function join_boxes(A::Box, B::Box, i::Int, z::Vector{Float64}, zᵢ::Vector{Float64}, zₙ::Vector{Float64})
	lᵃ, uᵃ, lᵇ, uᵇ = A.l, A.u, B.l, B.u
    @assert all(uᵃ .≤ uᵇ) "`join` operation not valid. (uᵃ ≰ uᵇ)"
    lᶜ, uᶜ = min.(lᵃ, lᵇ), uᵇ
    ẑ = max.(z, lᶜ)
	return Box(lᶜ, uᶜ, reduced_scaled_priority(lᶜ, uᶜ, i, ẑ, zᵢ, zₙ))
end

function select_next_box(L::Vector{Vector{Box}}, k::Int, p::Int)
    if any(!isempty(Lᵢ) for Lᵢ in L)
		k = k % p + 1
		while isempty(L[k])
			k = k % p + 1
		end
		i = argmax(B.priority for B in L[k])
	end
	return i, k
end

function update!(L::Vector{Vector{Box}}, z::Vector{Float64}, zᵢ::Vector{Float64}, zₙ::Vector{Float64})
	T = [Box[] for _ in 1:length(L)]
	for j in 1:length(L)
		for B in L[j]
			if all(z .< B.u)
				@info "Partitioning $(B.l), $(B.u)..."
				for (i, Bᵢ) in enumerate(p_partition(B, z, zᵢ, zₙ))
					if !isempty(Bᵢ)
						@info "Pushing $(Bᵢ.l), $(Bᵢ.u) into T$i..."
						push!(T[i], Bᵢ)
					end
				end
			else
				push!(T[j], B)
			end
		end
	end
	L .= T
	for k in 1:length(L)
		i = 1; N = length(L[k])
		while i < N
			index_to_remove = []
			for j = i:N
				if i ≠ j
					if all(L[k][i].u .≤ L[k][j].u)
						L[k][i] = join_boxes(L[k][i], L[k][j], k, z, zᵢ, zₙ)
						push!(index_to_remove, j)
					elseif all(L[k][i].u .≥ L[k][j].u)
						L[k][i] = join_boxes(L[k][j], L[k][i], k, z, zᵢ, zₙ)
						push!(index_to_remove, j)
					end
				end
			end
			i += 1
			N -= length(index_to_remove)
			deleteat!(L[k], index_to_remove)
		end
	end
end

function dominguez(problem::AbstractMOProblem)
	Z = Vector{Float64}[]
	Xₑ = Array{Float64}[]
	L = [Box[] for _ in 1:problem.p]
	k = 0
    model, x, f = build_model(problem)
	zᵢ, zₙ = ideal(model, f, problem), nadir(model, f, problem)
    @info "Ideal point: $zᵢ"
    @info "Nadir point: $zₙ"
	r = compute_r(zᵢ, zₙ)
	ϵ = compute_ϵ(problem.p, r)
	@info "r: $r, ϵ: $ϵ"
	push!(L[1], Box(zᵢ, zₙ, 0.))
	iter = 1
	while any(!isempty(Lᵢ) for Lᵢ in L)
		@info "Iteration #$iter"
		i, k = select_next_box(L, k, problem.p)
		B = L[k][i]
		@info "Chosen Box: $(B.l), $(B.u). Direction: $i"
		w = compute_w(B.u, zᵢ)
		@info "w: $w"
		status, obj, z, xₑ = holtzman(model, f, x, problem, w, zᵢ, ϵ)
		@info "Objective: $obj, f(x): $z"
		if (obj < 1) && all(zᵢ .< B.u)
			push!(Z, z); push!(Xₑ, xₑ)
			update!(L, z, zᵢ, zₙ)
		else
			popat!(L[k], i)
		end
		iter += 1
	end
	return Z, Xₑ
end