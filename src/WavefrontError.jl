struct WavefrontError
    i::Vector{NamedTuple{(:j, :n, :m, :a), Tuple{Int, Int, Int, Float64}}}
    n_max::Int
    fit_to::Vector{Tuple{Int, Int}}
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

# Type aliases
const FloatVec = AbstractVector{<:AbstractFloat}
const FloatMat = AbstractMatrix{<:AbstractFloat}

function fit(ρ::FloatVec, θ::FloatVec, OPD::FloatVec, Zᵢ::Vector{Polynomial})
    if !allequal(length.((ρ, θ, OPD)))
        error("Vectors must be of equal length.\n")
    end
    # linear least squares
    A = stack(Z.(ρ, θ) for Z ∈ Zᵢ)
    # Zernike expansion coefficients
    v = A \ OPD
    return v, Zᵢ
end

# fitting function (construction function 1)
function Wf(ρ::FloatVec, θ::FloatVec, OPD::FloatVec, n_max::Int)
    j_max = get_j(n_max, n_max)
    Zᵢ = Polynomial[Z(j, Model()) for j = 0:j_max]
    fit(ρ, θ, OPD, Zᵢ)
end

function Wf(ρ::FloatVec, θ::FloatVec, OPD::FloatVec,
            orders::Vector{Tuple{Int, Int}})
    Zᵢ = Polynomial[Zf(mn...) for mn ∈ orders]
    fit(ρ, θ, OPD, Zᵢ)
end

# filtering function (construction function 2)
function Ψ(v, Zᵢ, n_max, orders = Tuple{Int, Int}[]; precision)
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
    # pad the output coefficient vector if needed
    if !isempty(orders)
        v = standardize(v, orders)
    end
    # create the fitted polynomial
    ΔW = WavefrontError(a, n_max, orders, av, Zₐ)
    return ΔW, a, v
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

function metrics(v::FloatVec, ΔWp::FloatMat)
    # Peak-to-valley wavefront error
    PV = maximum(ΔWp) - minimum(ΔWp)
    # RMS wavefront error
    # where σ² is the variance (second central moment about the mean)
    # the mean is the first a00 piston term
    σ = sqrt(v' * v - v[1]^2)
    Strehl_ratio = exp(-(2π * σ)^2)
    return (PV = PV, RMS = σ, Strehl = Strehl_ratio)
end

# main interface function
function W(ρ::FloatVec, θ::FloatVec, OPD::FloatVec, n_max::Int;
           precision = 3, scale = 101)
    v, Zᵢ = Wf(ρ, θ, OPD, n_max)
    ΔW, a = Ψ(v, Zᵢ, n_max; precision)
    Λ(ΔW, a, v, n_max; scale)
end

function W(ρ::FloatVec, θ::FloatVec, OPD::FloatVec, orders::Vector{Tuple{Int, Int}};
           precision = 3, scale = 101)
    n_max = maximum(mn -> mn[2], orders; init = 0)
    v, Zᵢ = Wf(ρ, θ, OPD, orders)
    ΔW, a, v = Ψ(v, Zᵢ, n_max, orders; precision)
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

# pads a subset Zernike expansion coefficient vector to standard length
function standardize(v_sub::FloatVec, orders::Vector{Tuple{Int, Int}})
    j = [get_j(mn...) for mn in orders]
    n_max = get_n(maximum(j))
    j_max = get_j(n_max, n_max)
    v_padded = zeros(eltype(v_sub), j_max + 1)
    v_padded[j.+1] .= v_sub
    return v_padded
end

# methods
function W(ρ::FloatVec, θ::FloatVec, OPD::FloatVec, n_max::Int, ::Model;
           precision = 3)
    v, Zᵢ = Wf(ρ, θ, OPD, n_max)
    return Ψ(v, Zᵢ, n_max; precision)[1]
end

function W(ρ::FloatVec, θ::FloatVec, OPD::FloatVec,
           orders::Vector{Tuple{Int, Int}}, ::Model; precision = 3)
    n_max = maximum(mn -> mn[2], orders; init = 0)
    v, Zᵢ = Wf(ρ, θ, OPD, orders)
    return Ψ(v, Zᵢ, n_max, orders; precision)[1]
end

W(; r, t, OPD, fit_to, options...) = W(r, t, OPD, fit_to; options...)

function W(x::FloatVec, y::FloatVec, OPD::FloatVec; fit_to, options...)
    ρ = hypot.(x, y)
    θ = atan.(y, x)
    W(ρ, θ, OPD, fit_to; options...)
end

# [ρ θ OPD] input format
# W(split_data(data)..., fit_to, etc.)
function split_data(data::FloatMat)
    (data[:, i] for i = 1:3)
end

# extract pupil coordinates
function coords(OPD::FloatMat)
    u, v = size(OPD)
    ρ = repeat(range(0.0, 1.0, v); inner = u)
    θ = repeat(range(0.0, 2π, u); outer = v)
    return ρ, θ
end

# assumes dim(θ) x dim(ρ) matrix polar mapping
# W(OPD', fit_to, etc.) for dim(ρ) x dim(θ) matrix
function W(OPD::FloatMat, fit_to; options...)
    W(coords(OPD)..., vec(OPD), fit_to; options...)
end

function W(OPD::FloatMat, fit_to, ::Model; precision = 3)
    W(coords(OPD)..., vec(OPD), fit_to, Model(); precision)
end
