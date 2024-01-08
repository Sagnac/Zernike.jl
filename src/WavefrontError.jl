const Recap = Vector{NamedTuple{(:j, :n, :m, :a), Tuple{Int, Int, Int, Float64}}}

const precision = 3
const max_precision = 17
const wavefront_finesse = 101

struct WavefrontError <: Phase
    recap::Recap
    v::Vector{Float64}
    n_max::Int
    fit_to::Vector{Tuple{Int, Int}}
    a::Vector{Float64}
    Z::Vector{Polynomial}
    precision::Int
end

struct WavefrontOutput
    recap::Recap
    v::Vector{Float64}
    metrics::NamedTuple{(:pv, :rms, :strehl), Tuple{Float64, Float64, Float64}}
    fig::Makie.Figure
    axis::Axis3
    plot::Surface{Tuple{T, T, T}} where T <: Matrix{Float32}
end

function WavefrontError(a::FloatVec)
    fit_to = []
    v = a
    recap = similar(a, NamedTuple)
    Zᵢ = similar(a, Polynomial)
    local n
    for (i, a) ∈ pairs(a)
        j = i - 1
        m, n = get_mn(j)
        recap[i] = (; j, n, m, a)
        Zᵢ[i] = Z(j, Model)
    end
    n_max = n
    return WavefrontError(recap, v, n_max, fit_to, a, Zᵢ, max_precision)
end

function WavefrontError(orders::Vector{Tuple{Int, Int}}, a::FloatVec)
    length(a) != length(orders) && throw(ArgumentError("Lengths must be equal."))
              allunique(orders) || throw(ArgumentError("Orders must be unique."))
    fit_to = []
    recap = similar(a, NamedTuple)
    Zᵢ = similar(a, Polynomial)
    local n
    for (idx, mn) ∈ pairs(orders)
        m, n = mn
        aᵢ = a[idx]
        j = get_j(m, n)
        recap[idx] = (; j, n, m, a = aᵢ)
        Zᵢ[idx] = Z(m, n, Model)
    end
    n_max = n
    v = standardize(a, orders)
    return WavefrontError(recap, v, n_max, fit_to, a, Zᵢ, max_precision)
end

function WavefrontError(orders::Vector{NamedTuple{(:m, :n), Tuple{Int, Int}}},
                        a::FloatVec)
    WavefrontError([(mn...,) for mn ∈ orders], a)
end

function WavefrontError(orders::Vector{Int}, a::FloatVec)
    WavefrontError([get_mn(j) for j ∈ orders], a)
end

function (ΔW::WavefrontError)(ρ, θ)
    (; a, Z) = ΔW
    ∑(aᵢ * Z[i](ρ, θ) for (i, aᵢ) ∈ pairs(a))
end

function fit(ρ::FloatVec, θ::FloatVec, OPD::FloatVec, Zᵢ::Vector{Polynomial})
    if !allequal(length.((ρ, θ, OPD)))
        error("Vectors must be of equal length.\n")
    end
    # linear least squares
    A = stack(Z.(ρ, θ) for Z ∈ Zᵢ)
    # Zernike expansion coefficients
    A \ OPD
end

function reconstruct(ρ::FloatVec, θ::FloatVec, OPD::FloatVec, n_max::Int)
    j_max = get_j(n_max, n_max)
    Zᵢ = Polynomial[Z(j, Model) for j = 0:j_max]
    v = fit(ρ, θ, OPD, Zᵢ)
    return v, Zᵢ
end

function reconstruct(ρ::FloatVec, θ::FloatVec, OPD::FloatVec,
                     orders::Vector{Tuple{Int, Int}})
    Zᵢ = Polynomial[construct(mn...) for mn ∈ orders]
    v = fit(ρ, θ, OPD, Zᵢ)
    return v, Zᵢ
end

# filtering / sifting function
function Ψ(v, Zᵢ, n_max, orders = Tuple{Int, Int}[]; precision::Int)
    recap = @NamedTuple{j::Int, n::Int, m::Int, a::Float64}[]
    a = Float64[]
    Zₐ = Polynomial[]
    # store the non-trivial coefficients
    for (i, aᵢ) in pairs(v)
        aᵢ = round(aᵢ; digits = precision)
        if !iszero(aᵢ)
            if !isassigned(Zᵢ, i)
                Zᵢ[i] = Z(i-1, Model)
            end
            push!(recap, (; Zᵢ[i].inds..., a = aᵢ))
            push!(a, aᵢ)
            push!(Zₐ, Zᵢ[i])
        end
    end
    # pad the output coefficient vector if needed
    if !isempty(orders)
        v = standardize(v, orders)
    end
    # create the fitted polynomial
    ΔW = WavefrontError(recap, v, n_max, orders, a, Zₐ, precision)
    return ΔW
end

# synthesis function
function Λ(ΔW; finesse::Int)
    (; recap, v, n_max) = ΔW
    finesse::Int = finesse ∈ 1:100 ? finesse : cld(100, √ length(recap))
    ρ, θ = polar(n_max, n_max; finesse)
    # construct the estimated wavefront error
    w = ΔW.(ρ', θ)
    W_LaTeX = format_strings(recap)
    titles = (plot_title = W_LaTeX, window_title = "Estimated wavefront error")
    fig, axis, plot = zplot(ρ, θ, w; titles...)
    WavefrontOutput(recap, v, metrics(v, w), fig, axis, plot)
end

function metrics(ΔW::WavefrontError)
    ρ, θ = polar()
    w = ΔW.(ρ', θ)
    metrics(ΔW.v, w)
end

function metrics(v::FloatVec, w::FloatMat)
    # Peak-to-valley wavefront error
    min_max = extrema(w)
    pv = min_max[2] - min_max[1]
    # RMS wavefront error
    # where σ² is the variance (second central moment about the mean)
    # the mean is the first a00 piston term
    v2 = @view v[begin+1:end]
    σ² = v2' * v2
    strehl_ratio = exp(-4π^2 * σ²)
    return (; pv, rms = √σ², strehl = strehl_ratio)
end

# main interface function
function W(ρ::FloatVec, θ::FloatVec, OPD::FloatVec, n_max::Int;
           precision::Int = precision, finesse::Int = wavefront_finesse)
    v, Zᵢ = reconstruct(ρ, θ, OPD, n_max)
    ΔW = Ψ(v, Zᵢ, n_max; precision)
    Λ(ΔW; finesse)
end

function W(ρ::FloatVec, θ::FloatVec, OPD::FloatVec, orders::Vector{Tuple{Int, Int}};
           precision::Int = precision, finesse::Int = wavefront_finesse)
    n_max = maximum(mn -> mn[2], orders; init = 0)
    v, Zᵢ = reconstruct(ρ, θ, OPD, orders)
    ΔW = Ψ(v, Zᵢ, n_max, orders; precision)
    Λ(ΔW; finesse)
end

# overload show to clean up the output
show(io::IO, W::T) where {T <: WavefrontError} = print(io, T, ": n_max = ", W.n_max)

function show(io::IO, m::MIME"text/plain", W::WavefrontError)
    show(io, W)
    haskey(io, :typeinfo) ? (return) : println(io)
    strip3 = map(-, displaysize(io), (3, 0))
    println(io, "    ∑aᵢZᵢ(ρ, θ):")
    show(IOContext(io, :limit => true, :displaysize => strip3), m, W.recap)
    print(io, "\n    --> ΔW(ρ, θ)")
end

function show(io::IO, W::T) where {T <: WavefrontOutput}
    print(io, T, "(")
    join(io, fieldnames(T), ", ")
    print(io, ")")
end

function show(io::IO, m::MIME"text/plain", W::WavefrontOutput)
    show(io, W)
    haskey(io, :typeinfo) && return
    strip3 = map(-, displaysize(io), (3, 0))
    println(io, "\nSummary:")
    show(IOContext(io, :limit => true, :displaysize => strip3), m, W.recap)
    println(io)
    show(IOContext(io, :compact => true), m, W.metrics)
    display(W.fig)
    return
end

# extend getindex to allow indexing the output
getindex(W::T, i = 1) where {T <: WavefrontOutput} = getfield(W, fieldnames(T)[i])

getindex(W::WavefrontError) = W.v

# hook into iterate to allow non-property destructuring of the output
iterate(W::WavefrontOutput, i = 1) = (i > 6 ? nothing : (W[i], i + 1))

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
function W(ρ::FloatVec, θ::FloatVec, OPD::FloatVec, n_max::Int, ::Type{Model};
           precision::Int = precision)
    v, Zᵢ = reconstruct(ρ, θ, OPD, n_max)
    return Ψ(v, Zᵢ, n_max; precision)
end

function W(ρ::FloatVec, θ::FloatVec, OPD::FloatVec,
           orders::Vector{Tuple{Int, Int}}, ::Type{Model};
           precision::Int = precision)
    n_max = maximum(mn -> mn[2], orders; init = 0)
    v, Zᵢ = reconstruct(ρ, θ, OPD, orders)
    return Ψ(v, Zᵢ, n_max, orders; precision)
end

W(; r, t, OPD, fit_to, options...) = W(r, t, OPD, fit_to; options...)

function W(x::FloatVec, y::FloatVec, OPD::FloatVec; fit_to, options...)
    ρ = hypot.(x, y)
    θ = atan.(y, x)
    W(ρ, θ, OPD, fit_to; options...)
end

# [ρ θ OPD] input format
# W(split_data(data)..., fit_to, etc.)
split_data(data::FloatMat) = (data[:, i] for i = 1:3)

# extract pupil coordinates
function coords(OPD::FloatMat)
    u, v = size(OPD)
    ρ = repeat(range(0.0, 1.0, v); inner = u)
    θ = repeat(range(0.0, 2π, u); outer = v)
    return ρ, θ
end

function coords(ρ::FloatVec, θ::FloatVec)
    ρ2 = ones(length(θ)) * ρ' |> vec
    θ2 = θ * ones(length(ρ))' |> vec
    return ρ2, θ2
end

# assumes dim(θ) x dim(ρ) matrix polar mapping
# W(OPD', fit_to, etc.) for dim(ρ) x dim(θ) matrix
function W(OPD::FloatMat, fit_to; options...)
    W(coords(OPD)..., vec(OPD), fit_to; options...)
end

function W(OPD::FloatMat, fit_to, ::Type{Model}; precision::Int = precision)
    W(coords(OPD)..., vec(OPD), fit_to, Model; precision)
end

function W(ρ::FloatVec, θ::FloatVec, OPD::FloatMat, fit_to; options...)
    W(coords(ρ, θ)..., vec(OPD), fit_to; options...)
end

function W(ρ::FloatVec, θ::FloatVec, OPD::FloatMat, fit_to,
           ::Type{Model}; precision::Int = precision)
    W(coords(ρ, θ)..., vec(OPD), fit_to, Model; precision)
end

# reverse transform;
# for input vectors corresponding to ordered triples over the exit pupil
# under the assumption of equally spaced uniform and regular sampling
# returns uniquely valued coordinate vectors and a phase matrix pinned by them
# such that OPD = ΔW.(ρ', θ)
function map_phase(ρ::FloatVec, θ::FloatVec, OPD::FloatVec)
    ρ = union(ρ)
    θ = union(θ)
    return ρ, θ, reshape(OPD, length.((θ, ρ)))
end
