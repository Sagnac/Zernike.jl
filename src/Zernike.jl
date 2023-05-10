# Zernike.jl
# Version 0.6.0
# 2023-5-8
# https://github.com/Sagnac/Zernike

# Generates Zernike polynomials, models wavefront errors, and plots them in Makie.

module Zernike

export Z, W, Coeffs, Latex, Fit

using Printf
import Base: @locals
import UnPack: @unpack

const ∑ = sum

include("ZernikePlot.jl")
include("WavefrontError.jl")
include("RadialCoefficients.jl")

struct Coeffs end
struct Latex end
struct Fit end

function polar(m, n)
    ϵ₁ = 100 * (ceil(Int, π * n) + 1)
    ϵ₁ = clamp(ϵ₁, ϵ₁, 1000)
    ϵ₂ = 100 * (2m + 1)
    ϵ₂ = clamp(ϵ₂, ϵ₂, 1000)
    ρ = range(0, 1, ϵ₁)
    θ = range(0, 2π, ϵ₂)
    return ρ, θ
end

get_j(n, m)::Integer = ((n + 2)n + m) ÷ 2

#=
# This is used in computing the polynomial coefficients using the original formula.
function fact(t)
    bound = Int ≡ Int32 ? 12 : 20
    [tᵢ > bound ? (factorial ∘ big)(tᵢ) : factorial(tᵢ) for tᵢ ∈ t]
end
# =#

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

    Z(ρ,θ)::Float64 = N * R(ρ) * M(θ)

    return Z, (j = j, n = n, m = m), Z_vars

end

# synthesis function
function Ψ(m::Integer, n::Integer)

    Z, (; j, n, m), Z_vars = Zf(m, n)

    @unpack γ = Z_vars

    ρ, θ = polar(m, n)

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

function format_strings(Z_vars::Dict)

    @unpack m, n, j, μ, k, N², N, λ, γ, ν = Z_vars

    γₛ = [@sprintf "%d" γₛ for γₛ ∈ abs.(γ)]

    UNICODE = ones(String, 3)
    LaTeX = ones(String, 3)

    # prefactor
    if !isinteger(N)
        UNICODE[1] = "√($N²)"
        LaTeX[1] = "\\sqrt{$N²}"
    elseif N ≠ 1
        UNICODE[1] = LaTeX[1] = @sprintf "%d" N
    end

    superscripts = ['\u2070'; '\u00B9'; '\u00B2'; '\u00B3'; '\u2074':'\u2079']

    ω = ones(String, length(ν))

    for (i, v) in pairs(ν)
        v < 2 && break
        for j in v |> digits |> reverse
            ω[i] *= superscripts[j+1]
        end
    end

    # polynomial terms
    ζ(i) = ν[i] == 1 ? "" : "^{$(ν[i])}"

    for i = 1:k
        UNICODE[2] *= string(γₛ[i], "ρ", ω[i], γ[i+1] > 0 ? " + " : " \u2212 ")
        LaTeX[2] *= string(γₛ[i], "\\rho", ζ(i), γ[i+1] > 0 ? " + " : " - ")
    end

    if ν[end] == 0
        UNICODE[2] *= γₛ[end]
        LaTeX[2] *= γₛ[end]
    elseif γ[end] == 1
        UNICODE[2] *= string("ρ", ω[end])
        LaTeX[2] *= string("\\rho", ν |> lastindex |> ζ)
    else
        UNICODE[2] *= string(γₛ[end], "ρ", ω[end])
        LaTeX[2] *= string(γₛ[end], "\\rho", ν |> lastindex |> ζ)
    end

    # angular term
    υ = μ == 1 ? "" : μ

    if m < 0
        UNICODE[3] = "sin($(υ)θ)"
        LaTeX[3] = "\\sin($(υ)\\theta)"
    elseif m > 0
        UNICODE[3] = "cos($(υ)θ)"
        LaTeX[3] = "\\cos($(υ)\\theta)"
    end

    parentheses = k ≠ 0 ? ("(", ")") : ""

    Zmn = "Z_{$n}^{$m}"

    Z_Unicode = join(UNICODE, parentheses...)
    Z_LaTeX = Zmn * " = " * join(LaTeX, parentheses...)

    return Zmn, Z_Unicode, Z_LaTeX

end

# main interface function
function Z(m::Integer, n::Integer)
    fig, = Ψ(m, n)
    return fig
end

# methods

function Z(m::Integer, n::Integer, ::Coeffs)
    fig, γ = Ψ(m, n)
    return fig, γ
end

function Z(m::Integer, n::Integer, ::Latex)
    fig, _, Z_LaTeX = Ψ(m, n)
    return fig, Z_LaTeX
end

function Z(m::Integer, n::Integer, ::Coeffs, ::Latex)
    return Ψ(m, n)
end

function Z(m::Integer, n::Integer, ::Fit)
    Z, j_n_m = Zf(m, n)
    return Z, j_n_m
end

Z(options...; m, n) = Z(m, n, options...)

function Z(j::Integer, options...)
    if j < 0
        error("j must be ≥ 0\n")
    end
    # radial order
    n::Integer = ceil((-3 + √(9 + 8j)) / 2)
    # azimuthal frequency
    m::Integer = 2j - (n + 2)n
    Z(m, n, options...)
end

end
