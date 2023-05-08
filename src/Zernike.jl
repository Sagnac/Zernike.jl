# Zernike.jl
# Version 0.5.0
# 2023-4-29
# https://github.com/Sagnac/Zernike

# Generates Zernike polynomials, fits wavefront errors, and plots them in Makie

# Julia v1.8.0

module Zernike

export Z, W

using GLMakie # v0.7.3
using Printf

const ρ = range(0, 1, 100)
const θ = range(0, 2π, 100)
const ∑ = sum

include("ZernikePlot.jl")
include("WavefrontError.jl")
include("RadialCoefficients.jl")

get_j(n, m)::Integer = ((n + 2)n + m) ÷ 2

#=
# This is used in computing the polynomial coefficients using the original formula.
function fact(t)
    bound = Int ≡ Int32 ? 12 : 20
    [tᵢ > bound ? (factorial ∘ big)(tᵢ) : factorial(tᵢ) for tᵢ ∈ t]
end
# =#

function Z(m::Integer, n::Integer; fit = false, coeffs = false, latex = false)

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

    # radial polynomial
    R(ρ)::Float64 = ∑(γ[i] * ρ ^ ν[i] for i = 1:k+1)

    # azimuthal / meridional component
    M(θ) = m < 0 ? -sin(m * θ) : cos(m * θ)

    Z(ρ,θ) = N * R(ρ) * M(θ)

    fit && return Z, (n = n, m = m)

    Zp = Z.(ρ', θ)

    indices = "j = $j, n = $n, m = $m"

    Zmn = "Z_{$n}^{$m}"

    window_title = "Zernike Polynomial: $indices"

    println(indices)

    if !latex && coeffs || n > 54
        fig = ZPlot(Zp; plot = Zmn, window = window_title)
        return fig, γ
    end

    # string formatting

    # there's probably an easier and more straightforward way to do this
    # but this approach is intuitive enough
    # and doesn't require importing another package

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

    Z_Unicode = join(UNICODE, parentheses...)
    Z_LaTeX = Zmn * " = " * join(LaTeX, parentheses...)

    # print and display

    n < 55 && println("Z = ", Z_Unicode)

    titles = (plot = j < 153 ? Z_LaTeX : Zmn, window = window_title)
    fig = ZPlot(Zp; titles...)

    if coeffs && latex
        return fig, γ, Z_LaTeX
    elseif coeffs
        return fig, γ
    elseif latex
        return fig, Z_LaTeX
    else
        return fig
    end

end

# methods

Z(; m, n, kwargs...) = Z(m, n; kwargs...)

function Z(j::Integer; kwargs...)
    if j < 0
        error("j must be ≥ 0\n")
    end
    # radial order
    n::Integer = ceil((-3 + √(9 + 8j)) / 2)
    # azimuthal frequency
    m::Integer = 2j - (n + 2)n
    Z(m, n; kwargs...)
end

end
