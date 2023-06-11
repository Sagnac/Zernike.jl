struct WavefrontError
    i::Vector{NamedTuple{(:j, :n, :m, :a), Tuple{Int, Int, Int, Float64}}}
    n_max::Int
    a::Vector{Float64}
    Z::Vector{Polynomial}
end

struct WavefrontOutput
    a::Vector{NamedTuple{(:j, :n, :m, :a), Tuple{Int, Int, Int, Float64}}}
    v::Vector{Float64}
    metrics::NamedTuple{(:PV, :RMS, :Strehl), Tuple{Float64, Float64, Float64}}
    fig::Makie.Figure
end

function (ΔW::WavefrontError)(ρ, θ)::Float64
    (; a, Z) = ΔW
    ∑(aᵢ * Z[i](ρ, θ) for (i, aᵢ) ∈ pairs(a))
end

# fitting function (construction function 1)
function Wf(ρ::Vector, θ::Vector, OPD::Vector, n_max::Int)
    if !allequal(length.((ρ, θ, OPD)))
        error("Vectors must be of equal length.\n")
    end
    j_max = get_j(n_max, n_max)
    Zᵢ = Polynomial[Z(j, Model()) for j = 0:j_max]
    # linear least squares
    A = reduce(hcat, Zᵢ[i].(ρ, θ) for i = 1:j_max+1)
    # Zernike expansion coefficients
    v::Vector{Float64} = A \ OPD
    return v, Zᵢ, n_max
end

# filtering function (construction function 2)
function Ψ(v, Zᵢ, n_max; precision)
    if precision ≠ "full" && !isa(precision, Int)
        precision = 3
    end
    a = @NamedTuple{j::Int, n::Int, m::Int, a::Float64}[]
    av = Float64[]
    Zₐ = Polynomial[]
    clipped = !isassigned(Zᵢ, 1)
    # store the non-trivial coefficients
    for (i, aᵢ) in pairs(v)
        aᵢ = precision == "full" ? aᵢ : round(aᵢ; digits = precision)
        if !iszero(aᵢ)
            if clipped
                Zᵢ[i] = Z(i-1, Model())
            end
            push!(a, (; Zᵢ[i].inds..., a = aᵢ))
            push!(av, aᵢ)
            push!(Zₐ, Zᵢ[i])
        end
    end
    # create the fitted polynomial
    ΔW = WavefrontError(a, n_max, av, Zₐ)
    return ΔW, a
end

# synthesis function
function Λ(ΔW, a, v, n_max; scale::Int)
    scale = scale ∈ 1:100 ? scale : ceil(Int, 100 / √ length(a))
    ρ, θ = polar(n_max, n_max; scale)
    # construct the estimated wavefront error
    ΔWp = ΔW.(ρ', θ)
    W_LaTeX = format_strings(a)
    titles = (plot = W_LaTeX, window = "Estimated wavefront error")
    fig = ZPlot(ρ, θ, ΔWp; titles...)
    WavefrontOutput(a, v, metrics(v, ΔWp), fig)
end

function metrics(v, ΔWp)
    # Peak-to-valley wavefront error
    PV = maximum(ΔWp) - minimum(ΔWp)
    # RMS wavefront error
    # where σ² is the variance (second central moment about the mean)
    # the mean is the first a00 piston term
    σ = ∑(v[i]^2 for i = 2:lastindex(v)) |> sqrt
    Strehl_ratio = exp(-(2π * σ)^2)
    return (PV = PV, RMS = σ, Strehl = Strehl_ratio)
end

# main interface function
function W(ρ::Vector, θ::Vector, OPD::Vector, n_max::Int; precision = 3, scale = 101)
    v, Zᵢ = Wf(ρ, θ, OPD, n_max)
    ΔW, a = Ψ(v, Zᵢ, n_max; precision)
    Λ(ΔW, a, v, n_max; scale)
end

# overload show to clean up the output
function show(io::IO, W::T) where {T <: WavefrontError}
    strip3 = map(-, displaysize(io), (3, 0))
    println(io, T, "(n_max = ", W.n_max, ")")
    println(io, "   ∑aᵢZᵢ(ρ, θ):")
    show(IOContext(io, :limit => true, :displaysize => strip3), "text/plain", W.i)
    print(io, "\n   --> ΔW(ρ, θ)")
end

function show(io::IO, W::T) where {T <: WavefrontOutput}
    strip3 = map(-, displaysize(io), (3, 0))
    println(io, T, "(a, v, metrics, fig)")
    println(io, "Summary:")
    show(IOContext(io, :limit => true, :displaysize => strip3), "text/plain", W.a)
    println(io)
    show(IOContext(io, :compact => true), "text/plain", W.metrics)
end

# extend getindex to allow indexing the output
getindex(W::T, i) where {T <: WavefrontOutput} = getfield(W, fieldnames(T)[i])

# hook into iterate to allow non-property destructuring of the output
iterate(W::WavefrontOutput, i = 1) = (i > 4 ? nothing : (W[i], i + 1))

# methods
function W(ρ::Vector, θ::Vector, OPD::Vector, n_max::Int, ::Model; precision = 3)
    return Ψ(Wf(ρ, θ, OPD, n_max)...; precision)[1]
end

W(; r, t, OPD, n_max, options...) = W(r, t, OPD, n_max; options...)

function W(x::Vector, y::Vector, OPD::Vector; n_max::Int, options...)
    ρ = hypot.(x, y)
    θ = atan.(y, x)
    W(ρ, θ, OPD, n_max; options...)
end

function W(data::Matrix, n_max::Int; options...)
    W((data[:, i] for i = 1:3)..., n_max; options...)
end

W(data; n_max, options...) = W(data, n_max; options...)
