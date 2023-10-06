struct Rectangle
    l::Vector{Float64}
    u::Vector{Float64}
    dim::Int
    function Rectangle(l::Vector{Float64}, u::Vector{Float64})
        @assert length(l) == length(u) "Dimension mismatch!"
        @assert u ⪴ l "u $(u) should be in the upper right of l $(l)!"
        new(l, u, length(l))
    end
end

function volume(r::Rectangle)
    return prod(r.u - r.l)
end

function volume(r::Rectangle, l::Vector{Float64})
    return prod(r.u - l)
end

function ⊆(Rᵢ::Rectangle, Rⱼ::Rectangle)
	@assert Rᵢ.dim == Rⱼ.dim "Dimension mismatch"
	return Rᵢ.l ⪴ Rⱼ.l && Rᵢ.u ⪳ Rⱼ.u
end

function Base.isempty(r::Rectangle)
	return volume(r) == 0
end

function split_rectangle(r::Rectangle, axis::Int, f::Float64)
	l_new = [i != axis ? r.l[i] : f for i = 1:r.dim]
	u_new = [i != axis ? r.u[i] : f for i = 1:r.dim]
	return Rectangle(r.l, u_new), Rectangle(l_new, r.u)
end

function update_list(L::Vector{Rectangle}, f::Vector{Float64})
	L̄, L = L, Rectangle[]
	for R_i in L̄
		l_i, u_i = R_i.l, R_i.u
		T = [R_i]
		for j = 1:length(f)
			if l_i[j] < f[j] < u_i[j]
				T̄ = Rectangle[]
				for R_t in T
					# @info "Splitting rectangle ($(R_t.l), $(R_t.u))..."
					push!(T̄, split_rectangle(R_t, j, f[j])...)
				end
				T = T̄
			end
		end
		push!(L, T...)
	end
	return L
end

function remove_rect!(L::Vector{Rectangle}, R::Rectangle)
	index_to_remove = []
	for (t, R_t) in enumerate(L)
		if R_t ⊆ R
			@info "Removing ($(R_t.l), $(R_t.u))..."
			push!(index_to_remove, t)
		end
	end
	deleteat!(L, index_to_remove)
end

function kirlik(problem::AbstractMOProblem, k::Int, δ::Int = 1)
    L_history = Int[]
    model, x, f = build_model(problem)
	yI, yN = round.(ideal(model, f, problem)), round.(nadir(model, f, problem))
	@info "Ideal point: $yI, Nadir point: $yN"
	L, Y_N, X_E = Rectangle[], Vector{Float64}[], Array{Float64}[]
	push!(L, Rectangle(project(yI, k), project(yN, k)))
	n_iter = 1
    n_infeasible = 0
	while !isempty(L)
		@info "Begining iteration #$(n_iter)"
		@info "Size of L: $(length(L))"
        push!(L_history, length(L))
		R = L[argmax(volume(R_i, project(yI, k)) for R_i in L)]
		@info "Rectangle chosen: ($(R.l), $(R.u))"
		@info "Solving P..."
        ε = insert!(copy(R.u), k, 0.0)
        ε_constraints = @constraint(model, [j = 1:problem.p; j ≠ k], f[j] <= ε[j] - δ)
        @objective(model, Min, f[k])
        optimize!(model)
        status = termination_status(model)
        @info "P is $status."
		if status == OPTIMAL
            y = round.(value.(f))
			@info "Found solution: $y"
			@info "Solving Q with yₖ=$(y[k])..."
            y_k_constraint = @constraint(model, f[k] == y[k])
            @objective(model, Min, sum(f))
            optimize!(model)
            status = termination_status(model)
			@info "Q is $status."
			if status == OPTIMAL
                y_star = round.(value.(f))
				@info "Found solution: $y_star"
				y_star_proj = project(y_star, k)
				if y_star in Y_N # solution already exists!
					@info "Solution already exists!"
				else # new solution!
					@info "New solution!"
                    x_star = value.(x)
					push!(Y_N, y_star); push!(X_E, x_star)
					@info "Updating L..."
					L = update_list(L, y_star_proj)
				end
				remove_rect!(L, Rectangle(project(y_star, k), R.u))
			else
				@info "Q is infeasible!"
				remove_rect!(L, Rectangle(project(yI, k), R.u))
			end
            delete(model, y_k_constraint)
		else # infeasible!
            n_infeasible += 1
			@info "P is infeasible!"
			remove_rect!(L, Rectangle(project(yI, k), R.u))
		end
        delete.(model, ε_constraints)
		n_iter += 1
	end
    Y_N, X_E, L_history, n_iter, n_infeasible
end
