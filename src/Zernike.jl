# Zernike.jl
# Base Version 6.2.0
# 2025-10-22
# https://github.com/Sagnac/Zernike.jl/tree/base

# Generates Zernike polynomials, & models wavefront errors.

module Zernike

export Z, W, Y, Wavefront, RadialPolynomial,
       get_j, get_m, get_n, get_mn, Noll, Fringe, Standard,
       noll_to_j, j_to_noll, fringe_to_j, j_to_fringe, standardize,
       reduce_wave, mnv, derivatives, ∇, OTF, MTF, PSF

const public_names = "public \
    radial_coefficients, wavefront_coefficients, transform_coefficients, \
    metrics, S, Superposition, Product, \
    sieve, format_strings, print_strings, valid_fringes, \
    grad, lap, Derivative, Gradient, Laplacian, Aberration"

VERSION >= v"1.11.0-DEV.469" && eval(Meta.parse(public_names))

using LaTeXStrings: @L_str, latexstring, LaTeXString
import Base: show, getindex, setindex!, firstindex, lastindex,
             getproperty, propertynames, complex

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

include("Domain.jl")
include("IndexConversions.jl")
include("Wavefront.jl")
include("Coordinates.jl")
include("RadialCoefficients.jl")
include("FormatStrings.jl")
include("Arithmetic.jl")
include("ScaleAperture.jl")
include("TransformAperture.jl")
include("Derivatives.jl")
include("MMT.jl")
include("MTF.jl")
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
    @domain_check(m, n)
    # upper bound for the sum (number of terms - 1 [indexing from zero])
    k = (n - μ) ÷ 2
    # OSA single mode index
    j = get_j(m, n)
    # normalization constant following the orthogonality relation
    N = √N²(m, n)
    # power (exponent)
    ν = collect(n:-2:μ)
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

# overload show to clean up the output
show(io::IO, Z::T) where {T <: Polynomial} = print(io, T, ": ", Z.inds, " ↦ Z(ρ, θ)")

function complex(Z::AbstractPolynomial)
    (; N, R, M) = Z
    (; m) = M
    function(s::Complex)
        ρ, θ = polar(s)
        ρ > 1.0 && return 0.0im
        N * R(ρ) * ei(m * θ)
    end
end

function getproperty(Z::AbstractPolynomial, name::Symbol)
    name in (:j, :m, :n) ? getfield(getfield(Z, :inds), name) : getfield(Z, name)
end

propertynames(::T) where T <: AbstractPolynomial = fieldnames(T)..., :j, :m, :n

getindex(Z::AbstractPolynomial) = Z.R.λ

getindex(Z::AbstractPolynomial, i) = Z.R.λ[i.+1]

firstindex(Z::AbstractPolynomial) = 0

lastindex(Z::AbstractPolynomial) = lastindex(Z.R.λ) - 1

Z(j::Int) = Z(get_mn(j)...)

const piston = Z(0, 0)

const valid_fringes = fringe_to_j.(1:37)

# API namespace
function radial_coefficients(m::Int, n::Int, T::Type{<:Number} = Float64)
    Φ(m, n, T)[end]
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
