# Zernike.jl
# Version 0.6.0
# 2023-5-8
# https://github.com/Sagnac/Zernike

# Generates Zernike polynomials, models wavefront errors, and plots them in Makie.

module Zernike

export Z, W, S, Coeffs, Latex, Fit

import Base: @locals
import UnPack: @unpack

const ∑ = sum

struct Coeffs end
struct Latex end
struct Fit end

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

include("ZernikePlot.jl")
include("WavefrontError.jl")
include("RadialCoefficients.jl")
include("FormatStrings.jl")
include("ScaleAperture.jl")

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
get_j(n, m)::Int = ((n + 2)n + m) ÷ 2

function radicand(m, n)
    # Kronecker delta δ_{m0}
    δ(m) = m == 0
    (2n + 2) ÷ (1 + δ(m))
end

#=
# This is used in computing the polynomial coefficients using the original formula.
function fact(t)
    bound = Int ≡ Int32 ? 12 : 20
    [tᵢ > bound ? (factorial ∘ big)(tᵢ) : factorial(tᵢ) for tᵢ ∈ t]
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

    j = get_j(n, m)

    N² = radicand(m, n)

    # normalization constant following the orthogonality relation
    N = √N²

    #=
    # This is the naive approach which implements the original explicit formula.
    function λ(s)
        t::Vector{Int} = [
            n - s;
            s;
            k - s;
            k + μ - s
        ]
        τ = t |> fact
        if j > 278
            τ = convert(Vector{BigInt}, τ)
        end
        (-1)^s * τ[1] ÷ prod(τ[2:4])
    end

    γ = Float64[λ(s) for s = 0:k]
    # =#

    # power (exponent)
    ν = Int[n - 2s for s = 0:k]

    # polynomial coefficients
    λ = Φ(n, μ)

    γ = Float64[λ[νᵢ+1] for νᵢ in ν]

    Z_vars = @locals

    inds = (j = j, n = n, m = m)

    # radial polynomial
    R = RadialPolynomial(γ, ν, k)

    # azimuthal / meridional component
    M = Sinusoid(m)

    # Zernike polynomial
    Z = Polynomial(inds, N, R, M)

    return Z, Z_vars

end

# synthesis function
function Ψ(m::Int, n::Int; scale::Int = 100)

    Z, Z_vars = Zf(m, n)

    @unpack γ = Z_vars

    ρ, θ = polar(m, n; scale)

    Zp = Z.(ρ', θ)

    Zmn, Z_Unicode, Z_LaTeX = format_strings(Z_vars)

    indices = replace(Z.inds |> string, '(':')' => "")
    window_title = "Zernike Polynomial: $indices"
    println(indices)

    if n > 54
        @info "Pass Coeffs() as the third argument to return the coefficients."
    else
        println("Z = ", Z_Unicode)
    end

    titles = (plot = Z.inds.j < 153 ? Z_LaTeX : Zmn, window = window_title)

    fig = ZPlot(ρ, θ, Zp; titles...)

    return fig, γ, Z_LaTeX

end

# main interface function
function Z(m::Int, n::Int; scale::Int = 100)
    fig, = Ψ(m, n; scale)
    return fig
end

# methods

function Z(m::Int, n::Int, ::Coeffs; scale::Int = 100)
    fig, γ = Ψ(m, n; scale)
    return fig, γ
end

function Z(m::Int, n::Int, ::Latex; scale::Int = 100)
    fig, _, Z_LaTeX = Ψ(m, n; scale)
    return fig, Z_LaTeX
end

function Z(m::Int, n::Int, ::Coeffs, ::Latex; scale::Int = 100)
    return Ψ(m, n; scale)
end

function Z(m::Int, n::Int, ::Fit)
    return Zf(m, n)[1]
end

function Z(opts...; m, n, options...)
    if !isempty(opts) && opts[1] isa Fit
        Zf(m, n)[1]
    else
        Z(m, n, opts...; options...)
    end
end

function Z(j::Int, opts...; options...)
    if j < 0
        error("j must be ≥ 0\n")
    end
    n = get_n(j)
    m = get_m(j, n)
    if !isempty(opts) && opts[1] isa Fit
        Zf(m, n)[1]
    else
        Z(m, n, opts...; options...)
    end
end

end
