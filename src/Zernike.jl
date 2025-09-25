# Zernike.jl
# Version 5.6.0
# 2025-09-04
# https://github.com/Sagnac/Zernike.jl

# Generates Zernike polynomials, models wavefront errors, and plots them using Makie.

module Zernike

export zernike, wavefront, transform, Z, W, Y, Wavefront,
       get_j, get_m, get_n, get_mn, Noll, Fringe, Standard,
       noll_to_j, j_to_noll, fringe_to_j, j_to_fringe, standardize,
       Observable, plotconfig, zplot, Screen, reduce_wave

const public_names = "public \
    radial_coefficients, wavefront_coefficients, transform_coefficients, \
    metrics, scale, S, Superposition, Product, \
    sieve, format_strings, print_strings, valid_fringes, \
    derivatives, Derivative, Gradient, resize!, reset!"

VERSION >= v"1.11.0-DEV.469" && eval(Meta.parse(public_names))

using GLMakie
import .Makie: latexstring, LaTeXString, FigureAxisPlot
import Base: show, getindex, setindex!, firstindex, lastindex, getproperty, complex

const ∑ = sum
const ∏ = prod

⊗(x, y) = kron(x, y)

# magic numbers, limited by performance
const d_max = 1024
const d_fit = 21

const finesse = 100

# Type aliases
const FloatVec = AbstractVector{<:AbstractFloat}
const FloatMat = AbstractMatrix{<:AbstractFloat}

const SurfacePlot = Surface{Tuple{Matrix{Float64}, Matrix{Float64}, Matrix{Float32}}}

abstract type Phase end
abstract type AbstractPolynomial <: Phase end

struct RadialPolynomial
    λ::Vector{Float64}
    γ::Vector{Float64}
    ν::Vector{Int}
end

struct Harmonic
    m::Int
end

struct Polynomial <: AbstractPolynomial
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
include("Wavefront.jl")
include("IndexConversions.jl")
include("Polar.jl")
include("RadialCoefficients.jl")
include("FormatStrings.jl")
include("Arithmetic.jl")
include("ZernikePlot.jl")
include("ScaleAperture.jl")
include("TransformAperture.jl")
include("Derivatives.jl")
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

function (Z::AbstractPolynomial)(ρ::Real, θ::Real = 0)
    (; N, R, M) = Z
    N * R(ρ) * M(θ)
end

(Z::AbstractPolynomial)(xy::Complex) = Z(polar(xy)...)

(Z::Output)(t...) = Z.Z(t...)

# radial order
get_n(j::Int) = ceil(Int, (-3 + √(9 + 8j)) / 2)

# azimuthal frequency
get_m(j::Int, n::Int) = 2j - (n + 2)n

get_m(j::Int) = get_mn(j)[1]

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

N²(m::Int, n::Int) = div(2n + 2, 1 + iszero(m))

# parity sign
psgn(k::Int) = 1 | -(k & 1) # (-1) ^ k

# This is the naive approach which implements the original explicit formula
# for computing the polynomial coefficients.
function canonical(μ::T, n::T, k::T) where T <: Int
    γ = Vector{Float64}(undef, k + 1)
    for s = 0:k
        t1::Float64 = factorial(n - s)
        t2::Float64 = factorial(s)
        t3::Float64 = factorial(k - s)
        t4::Float64 = factorial(k + μ - s)
        γ[s+1] = psgn(s) * t1 / (t2 * t3 * t4)
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
    N = √N²(m, n)
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

function (Z::AbstractPolynomial)(s::Complex)
    ρ, θ = polar(s)
    ρ > 1.0 && return 0.0im
    (; N, R, M) = Z
    (; m) = M
    N * R(ρ) * ei(m * θ)
end

getindex(Z::AbstractPolynomial) = Z.R.λ

getindex(Z::AbstractPolynomial, i) = Z.R.λ[i.+1]

firstindex(Z::AbstractPolynomial) = 0

lastindex(Z::AbstractPolynomial) = lastindex(Z.R.λ) - 1

# methods
zernike(; m, n, finesse = finesse) = zernike(m, n; finesse)

zernike(j::Int; finesse = finesse) = zernike(get_mn(j)...; finesse)

Z(j::Int) = Z(get_mn(j)...)

const piston = Z(0, 0)

const valid_fringes = fringe_to_j.(1:37)

# API namespace
function radial_coefficients(m::Int, n::Int, T::Type{<:Number} = Float64)
    Φ(abs(m), n, T)[end]
end

wavefront_coefficients(x...) = reconstruct(x...)[1:2]

function transform_coefficients(
    v::Vector{Float64},
    ε::Float64,
    δ::ComplexF64         =     0.0im,
    ϕ::Float64            =       0.0,
    ω::NTuple{2, Float64} = (1.0, 0.0)
)
    return Γ(v, ε, δ, ϕ, ω)[1]
end

end
