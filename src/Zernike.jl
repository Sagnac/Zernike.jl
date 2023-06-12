# Zernike.jl
# Version 1.2.3
# 2023-6-5
# https://github.com/Sagnac/Zernike.jl

# Generates Zernike polynomials, models wavefront errors, and plots them in Makie.

module Zernike

export Z, W, S, Model

using GLMakie
import GLMakie: Makie.latexstring, Makie.LaTeXString
import Base: show, getindex, iterate

const ∑ = sum

struct Model end

struct RadialPolynomial
    γ::Vector{Float64}
    ν::Vector{Int}
    k::Int
end

struct Sinusoid
    m::Int
end

struct Polynomial
    inds::NamedTuple{(:j, :n, :m), Tuple{Int64, Int64, Int64}}
    N::Float64
    R::RadialPolynomial
    M::Sinusoid
end

struct Output
    fig::Makie.Figure
    coeffs::Vector{Float64}
    latex::LaTeXString
end

include("ZernikePlot.jl")
include("WavefrontError.jl")
include("RadialCoefficients.jl")
include("FormatStrings.jl")
# include("ScaleAperture.jl")
include("TransformAperture.jl")
include("Docstrings.jl")

function (R::RadialPolynomial)(ρ)::Float64
    (; γ, ν, k) = R
    ∑(γ[i] * ρ ^ ν[i] for i = 1:k+1)
end

function (M::Sinusoid)(θ)::Float64
    (; m) = M
    m < 0 ? -sin(m * θ) : cos(m * θ)
end

function (Z::Polynomial)(ρ, θ)::Float64
    (; N, R, M) = Z
    N * R(ρ) * M(θ)
end

function polar(m, n; scale::Int = 100)
    m = abs(m)
    scale = clamp(scale, 1, 100)
    ϵ₁ = scale * (ceil(Int, π * n) + 1)
    ϵ₁ = clamp(ϵ₁, ϵ₁, 1000)
    ϵ₂ = scale * (2m + 1)
    ϵ₂ = clamp(ϵ₂, ϵ₂, 1000)
    ρ = range(0, 1, ϵ₁)
    θ = range(0, 2π, ϵ₂)
    return ρ, θ
end

# radial order
get_n(j)::Int = ceil((-3 + √(9 + 8j)) / 2)
# azimuthal frequency
get_m(j, n)::Int = 2j - (n + 2)n
# ISO / ANSI / OSA standard single mode-ordering index
get_j(m, n)::Int = ((n + 2)n + m) ÷ 2

function get_mn(j)
    if j < 0
        error("j must be ≥ 0\n")
    end
    n = get_n(j)
    m = get_m(j, n)
    return m, n
end

function radicand(m, n)
    # Kronecker delta δ_{m0}
    δ(m) = m == 0
    (2n + 2) ÷ (1 + δ(m))
end

#=
# This is the naive approach which implements the original explicit formula
# for computing the polynomial coefficients.
function fact(t::Float64)
    prod(2.0:t)
end

function λ(μ, n, s, k)
    t::Vector{Float64} = [
        n - s;
        s;
        k - s;
        k + μ - s
    ]
    τ = t .|> fact
    (-1)^s * τ[1] ÷ prod(τ[2:4])
end
# =#

# computation and construction function
# binds the indices and produces a specific polynomial function
function Zf(m::Int, n::Int)
    μ = abs(m)
    # validate
    if n < 0 || μ > n || isodd(n - μ)
        error("Bounds:\nn ≥ 0\n|m| ≤ n\nn - |m| even\n")
    end
    # upper bound for the sum (number of terms -1 [indexing from zero])
    k = (n - μ) ÷ 2
    j = get_j(m, n)
    # normalization constant following the orthogonality relation
    N = √radicand(m, n)
    # power (exponent)
    ν = Int[n - 2s for s = 0:k]
    # polynomial coefficients
    λ = Φ(μ, n)
    γ = Float64[λ[νᵢ+1] for νᵢ in ν]
    # γ = Float64[λ(μ, n, s, k) for s = 0:k]
    inds = (j = j, n = n, m = m)
    # radial polynomial
    R = RadialPolynomial(γ, ν, k)
    # azimuthal / meridional component
    M = Sinusoid(m)
    # Zernike polynomial
    Z = Polynomial(inds, N, R, M)
    return Z
end

# main interface function
function Z(m::Int, n::Int; scale::Int = 100)
    Z = Zf(m, n)
    (; γ) = Z.R
    ρ, θ = polar(m, n; scale)
    Zp = Z.(ρ', θ)
    Zmn, Z_LaTeX, Z_Unicode = format_strings(Z)
    indices = replace(Z.inds |> string, '(':')' => "")
    window_title = "Zernike Polynomial: $indices"
    println(indices)
    if n > 54
        println()
        @info "Coefficients are stored in the coeffs field of the current output."
    else
        print("Z = ", Z_Unicode)
    end
    titles = (plot = Z.inds.j < 153 ? Z_LaTeX : Zmn, window = window_title)
    fig = ZPlot(ρ, θ, Zp; titles...)
    Output(fig, γ, Z_LaTeX)
end

# overload show to clean up the output
show(io::IO, Z::T) where {T <: Polynomial} = print(io, T, Z.inds, " --> Z(ρ, θ)")
show(::IO, ::Output) = nothing

# extend getindex to allow indexing the output
getindex(Z::T, i) where {T <: Output} = getfield(Z, fieldnames(T)[i])

# hook into iterate to allow non-property destructuring of the output
iterate(Z::Output, i = 1) = (i > 3 ? nothing : (Z[i], i + 1))

# methods
function Z(m::Int, n::Int, ::Model)
    Zf(m, n)
end

function Z(; m, n, scale::Int = 100)
    Z(m, n; scale)
end

function Z(j::Int; scale::Int = 100)
    Z(get_mn(j)...; scale)
end

function Z(j::Int, ::Model)
    Zf(get_mn(j)...)
end

end
