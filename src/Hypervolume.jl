using Libdl

lib = Libdl.dlopen(joinpath(@__DIR__, "hv.dylib"))
sym = Libdl.dlsym(lib, :fpli_hv)

function hypervolume(data, reference)
    n, d = size(data)
    @assert d == length(reference) "Dimension mismatch"
    hv = ccall(sym, Cdouble, (Ref{Cdouble}, Cint, Cint, Ref{Cdouble}), reduce(vcat, eachrow(data)), d, n, reference)
    return hv
end

hypervolume(data::AbstractMatrix{Float64}) = hypervolume(data, maximum(data, dims=1)')
hypervolume(data::Vector{Vector{Float64}}) = hypervolume(transpose(reduce(hcat, data)))
hypervolume(data::Vector{Vector{Float64}}, reference::AbstractVector{Float64}) = hypervolume(transpose(reduce(hcat, data)), reference)