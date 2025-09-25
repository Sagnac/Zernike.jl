const Recap = Vector{NamedTuple{(:j, :n, :m, :a), Tuple{Int, Int, Int, Float64}}}

const precision = 3
const max_precision = 17
# out-of-bounds sentinel used for dynamically setting default values:
const wavefront_finesse = 101

struct Wavefront{T <: Polynomial} <: Phase
    recap::Recap
    v::Vector{Float64}
    n_max::Int
    fit_to::Vector{Tuple{Int, Int}}
    a::Vector{Float64}
    Z::Vector{Polynomial}
    precision::Int
    ssr::Float64
end

struct WavefrontOutput
    recap::Recap
    v::Vector{Float64}
    ssr::Float64
    metrics::NamedTuple{(:pv, :rms, :strehl), NTuple{3, Float64}}
    W::Wavefront{Polynomial}
    fig::Makie.Figure
    axis::Axis3
    plot::SurfacePlot
end

function Wavefront(recap, v, n_max, fit_to, a, Z, precision, ssr)
    Wavefront{Polynomial}(recap, v, n_max, fit_to, a, Z, precision, ssr)
end

function Wavefront(recap, v, n_max, fit_to, a, Z, precision)
    Wavefront(recap, v, n_max, fit_to, a, Z, precision, 0.0)
end

function Wavefront(a::FloatVec; precision = max_precision)
    fit_to = []
    v = a
    recap = similar(a, NamedTuple)
    Zᵢ = similar(a, Polynomial)
    local n
    for (i, a) ∈ pairs(a)
        j = i - 1
        m, n = get_mn(j)
        recap[i] = (; j, n, m, a)
        Zᵢ[i] = Z(j)
    end
    n_max = n
    return Wavefront(recap, v, n_max, fit_to, a, Zᵢ, precision)
end

function Wavefront(orders::Vector{Tuple{Int, Int}}, a::FloatVec;
                        precision = max_precision)
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
        Zᵢ[idx] = Z(m, n)
    end
    n_max = n
    v = standardize(a, orders)
    return Wavefront(recap, v, n_max, fit_to, a, Zᵢ, precision)
end

function Wavefront(orders::Vector{NamedTuple{(:m, :n), Tuple{Int, Int}}},
                        a::FloatVec)
    Wavefront([(mn...,) for mn ∈ orders], a)
end

function Wavefront(orders::Vector{Int}, a::FloatVec)
    Wavefront([get_mn(j) for j ∈ orders], a)
end

function (ΔW::Wavefront)(ρ::Real, θ::Real = 0)
    (; a, Z) = ΔW
    ∑(aᵢ * Zᵢ(ρ, θ) for (aᵢ, Zᵢ) ∈ zip(a, Z); init = 0.0)
end

(ΔW::Wavefront)(xy::Complex) = ΔW(polar(xy)...)

(W::WavefrontOutput)(t...) = W.W(t...)

function fit(ρ::FloatVec, θ::FloatVec, OPD::FloatVec, Zᵢ::Vector{Polynomial})
    if !allequal(length.((ρ, θ, OPD)))
        error("Vectors must be of equal length.\n")
    end
    # linear least squares
    A = stack(Z.(ρ, θ) for Z ∈ Zᵢ)
    # Zernike expansion coefficients
    v = A \ OPD
    # residual vector
    e = OPD - A * v
    # sum of the squared residuals
    ssr = e'e
    return v, ssr
end

function reconstruct(ρ::FloatVec, θ::FloatVec, OPD::FloatVec, n_max::Int)
    j_max = get_j(n_max)
    Zᵢ = Polynomial[Z(j) for j = 0:j_max]
    v, ssr = fit(ρ, θ, OPD, Zᵢ)
    return v, ssr, Zᵢ
end

function reconstruct(ρ::FloatVec, θ::FloatVec, OPD::FloatVec,
                     orders::Vector{Tuple{Int, Int}})
    Zᵢ = Polynomial[Z(mn...) for mn ∈ orders]
    v, ssr = fit(ρ, θ, OPD, Zᵢ)
    return v, ssr, Zᵢ
end

function reconstruct(OPD::FloatMat, fit_to)
    reconstruct(coords(OPD)..., vec(OPD), fit_to)
end

function reconstruct(ρ::FloatVec, θ::FloatVec, OPD::FloatMat, fit_to)
    reconstruct(coords(ρ, θ)..., vec(OPD), fit_to)
end

# filtering / sifting function
function Ψ(v, Zᵢ, n_max, orders = Tuple{Int, Int}[], ssr = 0.0; precision::Int)
    recap = @NamedTuple{j::Int, n::Int, m::Int, a::Float64}[]
    a = Float64[]
    Zₐ = Polynomial[]
    converted = (-1, -1) ∈ orders
    effective_precision = converted ? max_precision : precision
    # store the non-trivial coefficients
    for (i, aᵢ) in pairs(v)
        aᵢ = round(aᵢ; digits = effective_precision)
        if !iszero(aᵢ)
            if !isassigned(Zᵢ, i)
                Zᵢ[i] = Z(i-1)
            end
            push!(recap, (; Zᵢ[i].inds..., a = aᵢ))
            push!(a, aᵢ)
            push!(Zₐ, Zᵢ[i])
        end
    end
    # pad the output coefficient vector if needed
    if !isempty(orders) && !converted
        v = standardize(v, orders)
    end
    isempty(recap) && push!(recap, (; piston.inds..., a = 0.0))
    # create the fitted polynomial
    return Wavefront(recap, v, n_max, orders, a, Zₐ, precision, ssr)
end

# synthesis function
function Λ(ΔW::Wavefront; finesse::Int)
    (; recap, v, n_max, ssr) = ΔW
    finesse = finesse ∈ 1:100 ? finesse : ceil(Int, 100 / √ length(recap))
    ρ, θ = polar(n_max, n_max; finesse)
    # construct the estimated wavefront error
    w = ΔW.(ρ', θ)
    W_LaTeX = format_strings(ΔW)
    titles = (plot_title = W_LaTeX, window_title = "Estimated wavefront error")
    fig, axis, plot = zplot(ρ, θ, w; titles...)
    WavefrontOutput(recap, v, ssr, metrics(v, w), ΔW, fig, axis, plot)
end

function metrics(ΔW::Wavefront)
    ρ, θ = polar(d_max)
    w = ΔW.(ρ, θ)
    metrics(ΔW.v, w)
end

function metrics(v::FloatVec, w::FloatMat)
    # Peak-to-valley wavefront error
    min_max = extrema(w)
    pv = min_max[2] - min_max[1]
    # RMS wavefront error
    # where σ² is the variance (second central moment about the mean)
    # and the mean is the first a00 piston term
    a = @view v[begin+1:end]
    σ² = a'a
    strehl_ratio = exp(-4π^2 * σ²)
    return (; pv, rms = √σ², strehl = strehl_ratio)
end

# main interface function
function wavefront(ρ::FloatVec, θ::FloatVec, OPD::FloatVec, n_max::Int;
                   precision::Int = precision, finesse::Int = wavefront_finesse)
    ΔW = W(ρ, θ, OPD, n_max; precision)
    Λ(ΔW; finesse)
end

function wavefront(ρ::FloatVec, θ::FloatVec, OPD::FloatVec,
                   orders::Vector{Tuple{Int, Int}};
                   precision::Int = precision, finesse::Int = wavefront_finesse)
    ΔW = W(ρ, θ, OPD, orders; precision)
    Λ(ΔW; finesse)
end

# overload show to clean up the output
show(io::IO, W::T) where {T <: Wavefront} = print(io, T, ": n_max = ", W.n_max)

function show(io::IO, m::MIME"text/plain", W::Wavefront)
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

getindex(W::Wavefront) = W.v

getindex(W::Wavefront, j) = W.v[j.+1]

getindex(W::Wavefront, m::Int, n::Int) = W[get_j(m, n)]

getindex(W::Wavefront, mn::NTuple{2}{Int}) = W[mn[1], mn[2]]

getindex(W::Wavefront, orders::Vector{NTuple{2}{Int}}) = W[get_j.(orders)]

firstindex(W::Wavefront) = 0

lastindex(W::Wavefront) = lastindex(W.v) - 1

function setindex!(W::Wavefront, x, j)
    for (i, t) in pairs(W.recap)
        if t.j == j
            W.recap[i] = merge(t, (; a = x))
            W.a[i] = x
        end
    end
    W.v[j.+1] = x
end

setindex!(W::Wavefront, x, m::Int, n::Int) = setindex!(W, x, get_j(m, n))

setindex!(W::Wavefront, x, mn::NTuple{2}{Int}) = (W[get_j(mn)] = x)

setindex!(W::Wavefront, x, orders::Vector{NTuple{2}{Int}}) = (W[get_j.(orders)] = x)

# reduces precision
function reduce_wave(W::Wavefront, precision::Int)
    @domain(precision < W.precision, precision)
    recap, a, Zᵢ = map(copy, (W.recap, W.a, W.Z))
    fit_to = getfield(W, :fit_to)
    for i ∈ (reverse ∘ eachindex)(a)
        aᵢ = round(a[i]; digits = precision)
        if iszero(aᵢ)
            deleteat!(recap, i)
            deleteat!(a, i)
            deleteat!(Zᵢ, i)
        else
            (; j, n, m) = recap[i]
            recap[i] = (; j, n, m, a = aᵢ)
            a[i] = aᵢ
        end
    end
    return Wavefront(recap, W.v, W.n_max, fit_to, a, Zᵢ, precision, W.ssr)
end

function sieve(v::Vector{Float64}, threshold::Float64)
    map(a -> abs(a) < threshold ? 0.0 : a, v)
end

# pads a subset Zernike expansion coefficient vector to standard length
function standardize(v_sub::FloatVec, j::AbstractVector{Int})
    n_max = get_n(maximum(j))
    j_max = get_j(n_max)
    v_padded = zeros(eltype(v_sub), j_max + 1)
    v_padded[j.+1] .= v_sub
    return v_padded
end

standardize(v_sub::FloatVec) = standardize(v_sub, 0:length(v_sub)-1)

function standardize(v_sub::FloatVec, orders::Vector{Tuple{Int, Int}})
    standardize(v_sub, get_j.(orders))
end

standardize(W::Wavefront) = standardize(W.a, [(i.m, i.n) for i ∈ W.recap])

# methods
function W(ρ::FloatVec, θ::FloatVec, OPD::FloatVec, n_max::Int;
           precision::Int = precision)
    v, ssr, Zᵢ = reconstruct(ρ, θ, OPD, n_max)
    return Ψ(v, Zᵢ, n_max, [], ssr; precision)
end

function W(ρ::FloatVec, θ::FloatVec, OPD::FloatVec,
           orders::Vector{Tuple{Int, Int}};
           precision::Int = precision)
    n_max = maximum(mn -> mn[2], orders; init = 0)
    v, ssr, Zᵢ = reconstruct(ρ, θ, OPD, orders)
    return Ψ(v, Zᵢ, n_max, orders, ssr; precision)
end

wavefront(; r, t, OPD, fit_to, options...) = wavefront(r, t, OPD, fit_to; options...)

function wavefront(x::FloatVec, y::FloatVec, OPD::FloatVec; fit_to, options...)
    ρ, θ = polar(x, y)
    wavefront(ρ, θ, OPD, fit_to; options...)
end

function W(x::FloatVec, y::FloatVec, OPD::FloatVec;
           fit_to, precision::Int = precision)
    ρ, θ = polar(x, y)
    W(ρ, θ, OPD, fit_to; precision)
end

# extract pupil coordinates
function coords(ρ::FloatVec, θ::FloatVec)
    ρ2 = ρ ⊗ ones(length(θ))
    θ2 = ones(length(ρ)) ⊗ θ
    return ρ2, θ2
end

coords(OPD::FloatMat) = coords(polar(size(OPD))...)

# assumes dim(θ) × dim(ρ) matrix polar mapping
# W(OPD', fit_to, etc.) for dim(ρ) × dim(θ) matrix
function wavefront(OPD::FloatMat, fit_to; options...)
    wavefront(coords(OPD)..., vec(OPD), fit_to; options...)
end

function W(OPD::FloatMat, fit_to; precision::Int = precision)
    W(coords(OPD)..., vec(OPD), fit_to; precision)
end

function wavefront(ρ::FloatVec, θ::FloatVec, OPD::FloatMat, fit_to; options...)
    wavefront(coords(ρ, θ)..., vec(OPD), fit_to; options...)
end

function W(ρ::FloatVec, θ::FloatVec, OPD::FloatMat, fit_to;
           precision::Int = precision)
    W(coords(ρ, θ)..., vec(OPD), fit_to; precision)
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
