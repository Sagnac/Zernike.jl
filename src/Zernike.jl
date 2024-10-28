# Zernike.jl
# Version 5.2.2
# 2024-10-23
# https://github.com/Sagnac/Zernike.jl

# Generates Zernike polynomials, models wavefront errors, and plots them using Makie.

module Zernike

export zernike, wavefront, transform, Z, W, P, WavefrontError, get_j, get_mn,
       Noll, Fringe, noll_to_j, j_to_noll, fringe_to_j, j_to_fringe, standardize,
       Standard, Observable, plotconfig, zplot, reduce_wave

const public_names = "public \
    radial_coefficients, wavefront_coefficients, transform_coefficients, \
    metrics, scale, J, Superposition, Product, sieve, format_strings, valid_fringes"

VERSION >= v"1.11.0-DEV.469" && eval(Meta.parse(public_names))

using GLMakie
import .Makie: latexstring, LaTeXString, FigureAxisPlot
import Base: show, getindex, iterate, getproperty, setproperty!, propertynames

const ∑ = sum
const ∏ = prod

⊗(x, y) = kron(x, y)

const ϵ_max = 2^10
const ϵ_fit = 21

const finesse = 100

# Type aliases
const FloatVec = AbstractVector{<:AbstractFloat}
const FloatMat = AbstractMatrix{<:AbstractFloat}

const SurfacePlot = Surface{Tuple{Matrix{Float64}, Matrix{Float64}, Matrix{Float32}}}

abstract type Phase end

struct RadialPolynomial
    λ::Vector{Float64}
    γ::Vector{Float64}
    ν::Vector{Int}
end

struct Harmonic
    m::Int
end

struct Polynomial <: Phase
    inds::NamedTuple{(:j, :n, :m), NTuple{3, Int}}
    N::Float64
    R::RadialPolynomial
    M::Harmonic
end

struct Output
    Z::Polynomial
    fig::Makie.Figure
    axis::Axis3
    plot::SurfacePlot
    coeffs::Vector{Float64}
    latex::LaTeXString
    unicode::String
    inds::String
    high_order::Bool
end

include("Domain.jl")
include("IndexConversions.jl")
include("Polar.jl")
include("RadialCoefficients.jl")
include("FormatStrings.jl")
include("WavefrontError.jl")
include("Arithmetic.jl")
include("ZernikePlot.jl")
include("ScaleAperture.jl")
include("TransformAperture.jl")
include("Docstrings.jl")

function (R::RadialPolynomial)(ρ::Real)
    (; λ) = R
    @domain(0 ≤ ρ ≤ 1, ρ)
    evalpoly(ρ, λ)
end

function (M::Harmonic)(θ::Real)
    (; m) = M
    m < 0 ? -sin(m * θ) : cos(m * θ)
end

function (Z::Polynomial)(ρ::Real, θ::Real = 0)
    (; N, R, M) = Z
    N * R(ρ) * M(θ)
end

(Z::Output)(ρ, θ) = Z.Z(ρ, θ)

# radial order
get_n(j::Int) = ceil(Int, (-3 + √(9 + 8j)) / 2)

# azimuthal frequency
get_m(j::Int, n::Int) = 2j - (n + 2)n

# ISO / ANSI / OSA standard single mode-ordering index
function get_j(m::Int, n::Int)
    μ = abs(m)
    @domain_check_mn
    return ((n + 2)n + m) ÷ 2
end

get_j((m, n)) = get_j(m, n)

get_j(n_max::Int) = get_j(n_max, n_max)

function get_mn(j::Int)
    @domain_check_j
    n = get_n(j)
    m = get_m(j, n)
    return m, n
end

function validate_length(v::Vector)
    len = length(v)
    j_max = len - 1
    n_max = get_n(j_max)
    if j_max ≠ get_j(n_max)
        error(
            """
            Invalid number of coefficients.
            Coefficients for Zernike polynomials up to m = n_max, n = n_max required.
            """
        )
    end
    return len, n_max
end

function radicand(m::Int, n::Int)
    # Kronecker delta δ_{m0}
    δ(m) = m == 0
    (2n + 2) ÷ (1 + δ(m))
end

# This is the naive approach which implements the original explicit formula
# for computing the polynomial coefficients.
function canonical(μ::T, n::T, k::T) where T <: Int
    γ = Vector{Float64}(undef, k + 1)
    for s = 0:k
        t1::Float64 = factorial(n - s)
        t2::Float64 = factorial(s)
        t3::Float64 = factorial(k - s)
        t4::Float64 = factorial(k + μ - s)
        γ[s+1] = (-1)^s * t1 / (t2 * t3 * t4)
    end
    return γ
end

# binds the indices and produces a specific polynomial function
function Z(m::Int, n::Int)
    μ = abs(m)
    # validate
    @domain_check_mn
    # upper bound for the sum (number of terms - 1 [indexing from zero])
    k = (n - μ) ÷ 2
    # OSA single mode index
    j = get_j(m, n)
    # normalization constant following the orthogonality relation
    N = √radicand(m, n)
    # power (exponent)
    ν = Int[n - 2s for s = 0:k]
    # polynomial coefficients
    if μ == n
        λ = [zeros(n); 1.0]
        γ = [1.0]
    elseif n < 21 && Int ≡ Int64
        γ = canonical(μ, n, k)
        λ = zeros(n+1)
        λ[ν.+1] = γ
    else
        λ = Φ(μ, n)[end]
        γ = Float64[λ[νᵢ+1] for νᵢ in ν]
    end
    inds = (j = j, n = n, m = m)
    # radial component
    R = RadialPolynomial(λ, γ, ν)
    # azimuthal / meridional component
    M = Harmonic(m)
    # Zernike polynomial
    return Polynomial(inds, N, R, M)
end

# main interface function
function zernike(m::Int, n::Int; finesse = finesse)
    Z_ = Z(m, n)
    (; γ) = Z_.R
    Z_mn, Z_LaTeX, Z_Unicode = format_strings(Z_)
    inds = chop(string(Z_.inds); head = 1, tail = 1)
    window_title = "Zernike Polynomial: $inds"
    high_order = n ≥ 48
    titles = (; plot_title = Z_.inds.j < 153 ? Z_LaTeX : Z_mn, window_title)
    fig, axis, plot = zplot(Z_; m, n, finesse, high_order, titles...)
    Output(Z_, fig, axis, plot, γ, Z_LaTeX, Z_Unicode, inds, high_order)
end

# overload show to clean up the output
show(io::IO, Z::T) where {T <: Polynomial} = print(io, T, Z.inds, " --> Z(ρ, θ)")

show(io::IO, output::Output) = print(io, output.inds)

function show(io::IO, ::MIME"text/plain", output::Output)
    show(io, output)
    haskey(io, :typeinfo) ? (return) : println(io)
    if output.high_order
        println(io)
        @info "Coefficients are stored in the coeffs field \
               of the current output." output.coeffs
    else
        print(io, "Z = ", output.unicode)
    end
    display(output.fig)
    return
end

getindex(Z::Polynomial) = Z.R.λ

getindex(Z::Polynomial, i::Int) = Z.R.λ[i+1]

# methods
zernike(; m, n, finesse = finesse) = zernike(m, n; finesse)

zernike(j::Int; finesse = finesse) = zernike(get_mn(j)...; finesse)

Z(j::Int) = Z(get_mn(j)...)

const piston = Z(0, 0)

const valid_fringes = fringe_to_j.(1:37)

# API namespace
radial_coefficients(x...) = Φ(x...)[end]

wavefront_coefficents(x...) = reconstruct(x...)[1]

function transform_coefficients(
    v::Vector{Float64},
    ε::Float64,
    δ::ComplexF64         =     0.0im,
    ϕ::Float64            =       0.0,
    ω::NTuple{2, Float64} = (1.0, 0.0)
)
    return S(v, ε, δ, ϕ, ω)[1]
end

end
