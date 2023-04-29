# Zernike.jl
# Version 0.5.0
# 2023-4-29
# https://github.com/Sagnac/Zernike

# Generates Zernike polynomials, fits wavefront errors, and plots them in Makie

# Julia v1.8.0

module Zernike

export Z, W

using GLMakie # v0.7.3

const ρ = range(0, 1, 100)
const θ = range(0, 2π, 100)
const ∑ = sum

include("ZernikePlot.jl")
include("WavefrontError.jl")

get_j(n, m)::Integer = ((n + 2)n + m) ÷ 2

function fact(t)
    bound = Int ≡ Int32 ? 12 : 20
    [tᵢ > bound ? (factorial ∘ big)(tᵢ) : factorial(tᵢ) for tᵢ ∈ t]
end

function Z(m::Integer, n::Integer; mode = "plot")

    μ::Integer = abs(m)

    # validate
    if mode != "plot" && mode != "fit"
        @error "Invalid mode.\nValid modes are:\n\"plot\"\n\"fit\""
        return
    end

    if n < 0 || μ > n || isodd(n - μ)
        @error "Bounds:\nn ≥ 0\n|m| ≤ n\nn - |m| even"
        return
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

    # coefficient
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

    γ = Integer[λ(s) for s = 0:k]

    # power (exponent)
    ν = Integer[n - 2s for s = 0:k]

    # radial polynomial
    R(ρ)::Float64 = ∑(γ[i] * ρ ^ ν[i] for i = 1:k+1)

    # azimuthal / meridional component
    M(θ) = m < 0 ? -sin(m * θ) : cos(m * θ)

    Z(ρ,θ) = N * R(ρ) * M(θ)

    mode == "fit" && return (Z = Z, n = n, m = m)

    Zp = Z.(ρ', θ)

    # string formatting

    # there's probably an easier and more straightforward way to do this
    # but this approach is intuitive enough
    # and doesn't require importing another package

    unicode = ones(String, 3)
    latex = ones(String, 3)

    # prefactor
    if !isinteger(N)
        unicode[1] = "√($N²)"
        latex[1] = "\\sqrt{$N²}"
    elseif N ≠ 1
        unicode[1] = latex[1] = string(Int(N))
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
        unicode[2] *= string(abs(γ[i]), "ρ", ω[i], γ[i+1] > 0 ? " + " : " \u2212 ")
        latex[2] *= string(abs(γ[i]), "\\rho", ζ(i), γ[i+1] > 0 ? " + " : " - ")
    end

    if ν[end] == 0
        unicode[2] *= string(abs(γ[end]))
        latex[2] *= string(abs(γ[end]))
    elseif γ[end] == 1
        unicode[2] *= string("ρ", ω[end])
        latex[2] *= string("\\rho", ν |> lastindex |> ζ)
    else
        unicode[2] *= string(abs(γ[end]), "ρ", ω[end])
        latex[2] *= string(abs(γ[end]), "\\rho", ν |> lastindex |> ζ)
    end

    # angular term
    υ = μ == 1 ? "" : μ

    if m < 0
        unicode[3] = "sin($(υ)θ)"
        latex[3] = "\\sin($(υ)\\theta)"
    elseif m > 0
        unicode[3] = "cos($(υ)θ)"
        latex[3] = "\\cos($(υ)\\theta)"
    end

    indices = "j = $j, n = $n, m = $m"

    parentheses = k ≠ 0 ? ("(", ")") : ""

    Z_Unicode = join(unicode, parentheses...)
    Z_LaTeX = "Z_{$n}^{$m} = " * join(latex, parentheses...)

    # print and display

    println(indices)
    println("Z = ", Z_Unicode)

    titles = (plot = Z_LaTeX, window = "Zernike Polynomial: $indices")
    ZPlot(Zp; titles...)

    return

end

# methods

Z(; m, n, mode = "plot") = Z(m, n; mode)

function Z(j::Integer; mode = "plot")
    if j < 0
        @error "j must be ≥ 0"
        return
    end
    # radial order
    n::Integer = ceil((-3 + √(9 + 8j)) / 2)
    # azimuthal frequency
    m::Integer = 2j - (n + 2)n
    Z(m, n; mode)
end

end
