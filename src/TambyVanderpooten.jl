function select_search_zone(
    U::Dict{Vector{Float64}, Vector{Set{Vector{Float64}}}}, 
    yI::Vector{Float64}, 
    yN::Vector{Float64}
)
    p = length(yI)
    upper_bounds = collect(keys(U))
    hvs = Dict((k, u) => u[k] == yN[k] ? -Inf : prod(project(u, k) - project(yI, k)) 
                for k in 1:p 
                for u in upper_bounds)
    k_star, u_star = argmax(hvs)
    return k_star, u_star
end

function get_child(u::Vector{Float64}, y_star::Vector{Float64}, k::Int)
    @assert length(u) == length(y_star)
    return vcat(u[begin:k-1], y_star[k], u[k+1:end])
end

function update_search_region(
    U::Dict{Vector{Float64}, Vector{Set{Vector{Float64}}}}, 
    y_star::Vector{Float64}, 
    yN::Vector{Float64}
)
    p = length(y_star)
    bounds_to_remove = Vector{Float64}[]
    bounds_to_add = Dict{Vector{Float64}, Vector{Set{Vector{Float64}}}}()
    for u in keys(U)
        if y_star ≺ u
            push!(bounds_to_remove, u)
            for l in 1:p
                u_l = get_child(u, y_star, l)
                N = [k == l ? Set([y_star]) : Set([y for y in U[u][k] if y[l] < y_star[l]]) for k in 1:p]
                if all(!isempty(N[k]) for k in 1:p if k != l && u_l[k] != yN[k])
                    bounds_to_add[u_l] = N
                end
            end
        else
            for k in 1:p
                if (y_star[k] == u[k]) && (project(y_star, k) ≺ project(u, k))
                    push!(U[u][k], y_star)
                end
            end
        end
    end
    for u in bounds_to_remove
        delete!(U, u)
    end
    merge!(U, bounds_to_add)
end

function apply_reduction_rule(
    U::Dict{Vector{Float64}, Vector{Set{Vector{Float64}}}}, 
    V::Vector{Dict{Float64, Vector{Tuple{Vector{Float64}, Vector{Float64}}}}}, 
    yI::Vector{Float64}, 
    y::Vector{Float64}
)
    p = length(y)
    bounds_to_remove = Vector{Float64}[]
    for u_i in keys(U)
        for k in 1:p
            if u_i[k] == yI[k]
                push!(bounds_to_remove, u_i)
                @info "$u_i has been added to the list of bounds to remove"
            else
                if haskey(V[k], u_i[k])
                    if any(project(u_i, k) ⪳ project(u_j, k) for (u_j, _) in V[k][u_i[k]])
                        push!(bounds_to_remove, u_i)
                        @info "$u_i has been added to the list of bounds to remove"
                    end
                end
            end
        end
    end
    for u in bounds_to_remove
        delete!(U, u)
    end
end

function tamby(problem::AbstractMOProblem, δ::Float64 = 1.)
    U_history = Int[]
    model, x, f = build_model(problem)
    yI = round.(ideal(model, f, problem))
    yN = round.(nadir(model, f, problem))
    yN = fill(maximum(yN), problem.p)
    @info "yI = $yI, yN = $yN"
    U = Dict{Vector{Float64}, Vector{Set{Vector{Float64}}}}()
    V = [Dict{Float64, Vector{Tuple{Vector{Float64}, Vector{Float64}}}}() for _ in 1:problem.p]
    U[yN] = [Set{Vector{Float16}}() for _ in 1:problem.p]
    solutions = Dict{Vector{Float64}, Array{Float64}}()
    n_iter = 1
    n_infeasible = 0
    while !isempty(U)
        push!(U_history, length(U))
        @info "Begining iteration #$(n_iter)"
        k, u = select_search_zone(U, yI, yN)
        @info "Search zone selected k, u = $k, $u"
        ε_constraints = @constraint(model, [j = 1:problem.p; j ≠ k], f[j] <= u[j] - δ)
        @objective(model, Min, f[k])
        if u != yN
            x0 = solutions[first(U[u][k])]
            set_start_value.(x, x0)
            @info "Starting solution set."
        end
        optimize!(model)
        status = termination_status(model)
        @info status
        if status == OPTIMAL
            y = round.(value.(f))
            @info "First stage solution y_k = $y"
            y_k_constraint = @constraint(model, f[k] == y[k])
            y_constraints = @constraint(model, [j = 1:problem.p; j ≠ k], f[j] <= y[j])
            delete.(model, ε_constraints)
            @objective(model, Min, sum(f))
            optimize!(model)
            y_star = round.(value.(f))
            x_star = value.(x)
            @info "Second stage solution y* = $y_star"
            if !haskey(V[k], y_star[k])
                V[k][y_star[k]] = [(u, y_star)]
            else
                push!(V[k][y_star[k]], (u, y_star))
            end
            @info "Added to the list of solved problems u, y* = $u, $y_star"
            if y_star ∉ U[u][k]
                @info "Found new nondominated points: y = $y_star"
                solutions[y_star] = x_star
                @info "Updating the search zone with y* = $y_star"
                update_search_region(U, y_star, yN)
            end
            delete.(model, y_constraints)
            delete(model, y_k_constraint)
            @info "Applying the reduction rule..."
            apply_reduction_rule(U, V, yI, y_star)
        else
            delete.(model, ε_constraints)
            @info "Model infeasible!"
            n_infeasible += 1
        end
        @info "Size of Y_N is $(length(keys(solutions)))"
        n_iter += 1
    end
    return keys(solutions), values(solutions), U_history, n_iter, n_infeasible
end