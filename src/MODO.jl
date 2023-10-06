using JuMP, HiGHS

function ≺(yi::Vector{Float64}, yj::Vector{Float64})
    @assert length(yi) == length(yj) "Dimension mismatch: $(length(yi)), $(length(yj))"
    return all(yi .< yj)
end

function ⪯(yi::Vector{Float64}, yj::Vector{Float64})
    @assert length(yi) == length(yj) "Dimension mismatch: $(length(yi)), $(length(yj))"
    return all(yi .<= yj)  && (yi ≠ yj)
end

function ⪳(yi::Vector{Float64}, yj::Vector{Float64})
    @assert length(yi) == length(yj) "Dimension mismatch: $(length(yi)), $(length(yj))"
    return all(yi .<= yj)
end

function ≻(yi::Vector{Float64}, yj::Vector{Float64})
    @assert length(yi) == length(yj) "Dimension mismatch: $(length(yi)), $(length(yj))"
    return all(yi .> yj)
end

function ⪰(yi::Vector{Float64}, yj::Vector{Float64})
    @assert length(yi) == length(yj) "Dimension mismatch: $(length(yi)), $(length(yj))"
    return all(yi .>= yj)  && (yi ≠ yj)
end

function ⪴(yi::Vector{Float64}, yj::Vector{Float64})
    @assert length(yi) == length(yj) "Dimension mismatch: $(length(yi)), $(length(yj))"
    return all(yi .>= yj)
end

function project(x::Vector{Float64}, axis::Int)
	return x[begin:end .!= axis]
end

abstract type AbstractMOProblem end

struct MOKnapsackProblem <: AbstractMOProblem
    p::Int
    n::Int
    W::Float64
    v::Matrix{Float64}
    w::Vector{Float64}
end

MOKnapsackProblem(filename::String) = open(filename) do f
    p = parse(Float64, readline(f))
    n = parse(Float64, readline(f))
    W = parse(Float64, readline(f))
    v = reduce(vcat, transpose(parse.(Float64, split(readline(f)))) for i in 1:p)
    w = parse.(Float64, split(readline(f)))
    MOKnapsackProblem(p, n, W, v, w)
end

struct MOAssignmentProblem <: AbstractMOProblem
    p::Int
    n::Int
    c::Array{Float64, 3}
end

MOAssignmentProblem(filename::String) = open(filename) do f
    p = parse(Int, readline(f))
    n = parse(Int, readline(f))
    c = permutedims(
        reshape(
            reduce(vcat, transpose.(parse.(Float64, split(readline(f))) for i in 1:p*n)), (n, p, n)
            ), 
            [2, 1, 3]
            )
            reduce(vcat, transpose.(parse.(Float64, split(readline(f))) for i in 1:p*n))
            MOAssignmentProblem(p, n, c)
        end

struct MOIntegerLinearProblem <: AbstractMOProblem
  p::Int
  n::Int
  m::Int
  c::Matrix{Float64}
  a::Matrix{Float64}
  b::Vector{Float64}
end

MOIntegerLinearProblem(filename::String) = open(filename) do f
    p = parse(Int, readline(f))
    n = parse(Int, readline(f))
    m = parse(Int, readline(f))
    c = reduce(vcat, transpose(parse.(Float64, split(readline(f)))) for i in 1:p)
    a = reduce(vcat, transpose(parse.(Float64, split(readline(f)))) for i in 1:m)
    b = parse.(Float64, split(readline(f)))
    MOIntegerLinearProblem(p, n, m, c, a, b)
end


function build_model(problem::MOKnapsackProblem)
	model = Model(optimizer_with_attributes(HiGHS.Optimizer))
    set_optimizer_attribute(model, MOI.Silent(), true)
	@variable(model, x[1:problem.n], Bin)
	@constraint(model, sum(problem.w .* x) ≤ problem.W)
	@expression(model, f[j = 1:problem.p], sum(-problem.v[j,:] .* x))
	return model, x, f
end

function build_model(problem::MOAssignmentProblem)
	model = Model(optimizer_with_attributes(HiGHS.Optimizer))
    set_optimizer_attribute(model, MOI.Silent(), true)
	@variable(model, x[1:problem.n, 1:problem.n], Bin)
	@constraint(model, [r = 1:problem.n], sum(x[r,:]) == 1)
	@constraint(model, [l = 1:problem.n], sum(x[:,l]) == 1)
	@expression(model, f[j = 1:problem.p], sum(problem.c[j,:,:] .* x))
	return model, x, f
end

function build_model(problem::MOIntegerLinearProblem)
	model = Model(optimizer_with_attributes(HiGHS.Optimizer))
    set_optimizer_attribute(model, MOI.Silent(), true)
	@variable(model, x[1:problem.n], Bin)
	@constraint(model, [r = 1:problem.m], sum(problem.a[r,:] .* x) ≤ problem.b[r])
	@expression(model, f[j = 1:problem.p], sum(-problem.c[j,:] .* x))
	return model, x, f
end

function solve_P(model::Model, problem::AbstractMOProblem, f::Vector{AffExpr}, k::Int, ε::Vector{Float64}, δ::Int = 1)
	ε̂ = insert!(copy(ε), k, 0)
    ε_constraints = @constraint(model, [j = 1:problem.p; j ≠ k], f[j] <= ε̂[j] - δ)
    set_objective_function(model, f[k])
	optimize!(model)
	status = termination_status(model)
	y = status == OPTIMAL ? value.(f) : nothing
    delete.(model, ε_constraints)
    return status, y
end

function solve_Q(model::Model, problem::AbstractMOProblem, f::Vector{AffExpr}, x::Vector{VariableRef}, 
    k::Int, y_k::Float64, ε::Vector{Float64}, δ::Int = 1)
	ε̂ = insert!(copy(ε), k, 0)
	ε_constraints = @constraint(model, [j = 1:problem.p; j ≠ k], f[j] <= ε̂[j] - δ)
    y_k_constraint = @constraint(model, f[k] == y_k)
	set_objective_function(model, f[k])
	optimize!(model)
	status = termination_status(model)
	y_star = status == OPTIMAL ? value.(f) : nothing
	x_star = status == OPTIMAL ? value.(x) : nothing
    delete.(model, ε_constraints)
    delete(model, y_k_constraint)
    return status, y_star, x_star
end

function solve_two_stage(
    problem::AbstractMOProblem, 
    k::Int, 
    ε::Vector{Float64}, 
    δ::Int = 1
)
    @assert k ≤ problem.p "Problem has $(problem.p) objectives."
    model, x, f = build_model(problem)
    status, y = solve_P(model, problem, f, k, ε, δ)
    if status == OPTIMAL
        status, y_star, x_star = solve_Q(model, problem, f, x, y[k], k, ε, δ)
        if status == OPTIMAL
            return status, y_star, x_star
        else
            return status, nothing, nothing
        end    
    else
        return status, nothing, nothing
    end
end

function solve(model::Model, f::AffExpr, problem::AbstractMOProblem, k::Int)
    @assert k ≤ problem.p "Problem has $(problem.p) objectives."
	set_objective_function(model, f[k])
	optimize!(model)
    return objective_value(model), value.(x)
end

function solve_weighted_sum(model::Model, f::Vector{AffExpr}, problem::AbstractMOProblem, w::Vector{Float64})
    @assert length(w) == problem.p "Problem has $(problem.p) objectives while w has length $(length(w))."
    set_objective_function(model, w'f)
	optimize!(model)
    return value.(f), value.(x), objective_value(model)
end

function ideal(model::Model, f::Vector{AffExpr}, problem::AbstractMOProblem)
	ideal = zeros(length(f))
	for j in 1:length(f)
		set_objective_function(model, f[j])
        set_objective_sense(model, MOI.MIN_SENSE)
		optimize!(model)
		ideal[j] = objective_value(model)
	end
	return ideal
end

function nadir(model::Model, f::Vector{AffExpr}, problem::AbstractMOProblem)
	nadir = zeros(length(f))
	for j in 1:length(f)
		set_objective_function(model, f[j])
        set_objective_sense(model, MOI.MAX_SENSE)
		optimize!(model)
		nadir[j] = objective_value(model)
	end
	return nadir
end

include("KirlikSayin.jl")
include("DominguezRios.jl")
include("TambyVanderpooten.jl")