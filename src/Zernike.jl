# Zernike.jl
# Version 0.6.0
# 2023-5-8
# https://github.com/Sagnac/Zernike

# Generates Zernike polynomials, models wavefront errors, and plots them in Makie.

module Zernike

export Z, W, Coeffs, Latex, Fit

import Base: @locals
import UnPack: @unpack

const ∑ = sum

struct Coeffs end
struct Latex end
struct Fit end

struct Polynomial
    N::Float64
    R::Function
    M::Function
end

include("ZernikePlot.jl")
include("WavefrontError.jl")
include("RadialCoefficients.jl")
include("FormatStrings.jl")

function (Z::Polynomial)(ρ, θ)::Float64
    (; N, R, M) = Z
    N * R(ρ) * M(θ)
end

function polar(m, n; scale::Integer = 100)
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
get_n(j)::Integer = ceil((-3 + √(9 + 8j)) / 2)
# azimuthal frequency
get_m(j, n)::Integer = 2j - (n + 2)n
# ANSI / OSA index
get_j(n, m)::Integer = ((n + 2)n + m) ÷ 2

#=
# This is used in computing the polynomial coefficients using the original formula.
function fact(t)
    bound = Int ≡ Int32 ? 12 : 20
    [tᵢ > bound ? (factorial ∘ big)(tᵢ) : factorial(tᵢ) for tᵢ ∈ t]
end
# =#

# computation and construction function
# binds the indices and produces a specific polynomial function
function Zf(m::Integer, n::Integer)

    μ::Integer = abs(m)

    # validate
    if n < 0 || μ > n || isodd(n - μ)
        error("Bounds:\nn ≥ 0\n|m| ≤ n\nn - |m| even\n")
    end

    # upper bound for the sum (number of terms -1 [indexing from zero])
    k::Integer = (n - μ) ÷ 2
    # ISO / ANSI / OSA standard single mode-ordering index
    j = get_j(n, m)
    # Kronecker delta δ_{m0}
    δ(m) = m == 0
    # radicand
    N²::Int = (2n + 2) ÷ (1 + δ(m))
    # normalization constant following the orthogonality relation
    N = √N²

    #=
    # This is the naive approach which implements the original explicit formula.
    function λ(s)
        t::Vector{Integer} = [
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

    # radial polynomial
    R(ρ)::Float64 = ∑(γ[i] * ρ ^ ν[i] for i = 1:k+1)

    # azimuthal / meridional component
    M(θ) = m < 0 ? -sin(m * θ) : cos(m * θ)

    Z = Polynomial(N, R, M)

    return Z, (j = j, n = n, m = m), Z_vars

end

# synthesis function
function Ψ(m::Integer, n::Integer; scale::Integer = 100)

    Z, (; j, n, m), Z_vars = Zf(m, n)

    @unpack γ = Z_vars

    ρ, θ = polar(m, n; scale)

    Zp = Z.(ρ', θ)

    Zmn, Z_Unicode, Z_LaTeX = format_strings(Z_vars)

    indices = "j = $j, n = $n, m = $m"
    window_title = "Zernike Polynomial: $indices"
    println(indices)

    if n > 54
        @info "Pass Coeffs() as the third argument to return the coefficients."
    else
        println("Z = ", Z_Unicode)
    end

    titles = (plot = j < 153 ? Z_LaTeX : Zmn, window = window_title)

    fig = ZPlot(ρ, θ, Zp; titles...)

    return fig, γ, Z_LaTeX

end

# main interface function
function Z(m::Integer, n::Integer; options...)
    fig, = Ψ(m, n; options...)
    return fig
end

# methods

function Z(m::Integer, n::Integer, ::Coeffs; options...)
    fig, γ = Ψ(m, n; options...)
    return fig, γ
end

function Z(m::Integer, n::Integer, ::Latex; options...)
    fig, _, Z_LaTeX = Ψ(m, n; options...)
    return fig, Z_LaTeX
end

function Z(m::Integer, n::Integer, ::Coeffs, ::Latex; options...)
    return Ψ(m, n; options...)
end

function Z(m::Integer, n::Integer, ::Fit)
    Z, j_n_m = Zf(m, n)
    return Z, j_n_m
end

Z(opts...; m, n, options...) = Z(m, n, opts...; options...)

function Z(j::Integer, options...)
    if j < 0
        error("j must be ≥ 0\n")
    end
    n = get_n(j)
    m = get_m(j, n)
    Z(m, n, options...)
end

end
