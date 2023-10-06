function coverage_error(R::AbstractMatrix{Float64}, N::AbstractMatrix{Float64})
	r = vec(maximum(N, dims=1) - minimum(N, dims=1))
    return maximum(
		minimum(
			[maximum(abs.(y .- x) ./ r) for x in eachrow(N), y in eachrow(R)], dims=2
		)
	)
end

coverage_error(R::Vector{Vector{Float64}}, N::Vector{Vector{Float64}}) = 
	coverage_error(transpose(reduce(hcat, R)), transpose(reduce(hcat, N)))

function ε₊(R::AbstractMatrix{Float64}, N::AbstractMatrix{Float64})
	r = vec(maximum(N, dims=1) - minimum(N, dims=1))
    return maximum(
		minimum(
			[maximum((y .- x) ./ r) for x in eachrow(N), y in eachrow(R)], dims=2
		)
	)
end

ε₊(R::Vector{Vector{Float64}}, N::Vector{Vector{Float64}}) = 
	ε₊(transpose(reduce(hcat, R)), transpose(reduce(hcat, N)))

function range_ratio(R::AbstractMatrix{Float64}, N::AbstractMatrix{Float64})
    rᴿ = (maximum(R, dims=1) - minimum(R, dims=1))
    rᴺ = (maximum(N, dims=1) - minimum(N, dims=1))
    return maximum(rᴿ ./ rᴺ)
end

range_ratio(R::Vector{Vector{Float64}}, N::Vector{Vector{Float64}}) = 
	range_ratio(transpose(reduce(hcat, R)), transpose(reduce(hcat, N)))

function onvgr(R::AbstractMatrix{Float64}, N::AbstractMatrix{Float64})
	size(R, 1) / size(N, 1)
end

onvgr(R::Vector{Vector{Float64}}, N::Vector{Vector{Float64}}) = length(R) / length(N)