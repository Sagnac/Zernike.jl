using Printf

# there's probably an easier and more straightforward way to do this
# but this approach is intuitive enough
# and doesn't require importing another package

const minus = '\u2212'
const superscripts = ['\u2070'; '\u00B9'; '\u00B2'; '\u00B3'; '\u2074':'\u2079']

function format_strings(Z::T) where T <: AbstractPolynomial
    (; j, n) = Z.inds
    N = Z.N
    (; γ, ν) = Z.R
    (; m) = Z.M
    (iszero(γ) || iszero(N)) && return L"\partial Z_{0}^{0}", L"0", "0"
    μ = abs(m)
    k = length(γ) - 1
    γ_s = [@sprintf("%d", abs(γᵢ)) for γᵢ ∈ γ]
    UNICODE = ones(String, 3)
    LaTeX = ones(String, 3)
    # prefactor
    isint = isinteger(N)
    if !isint
        sgn = sign(N) < 0 ? minus : ""
        N² = round(Int, N ^ 2)
        UNICODE[1] = "$(sgn)√($N²)"
        LaTeX[1] = "$(sgn)\\sqrt{$N²}"
    elseif N ≠ 1
        UNICODE[1] = LaTeX[1] = @sprintf "%d" N
    end
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
        UNICODE[2] *= string(γ_s[i], "ρ", ω[i], γ[i+1] > 0 ? " + " : " $minus ")
        LaTeX[2] *= string(γ_s[i], "\\rho", ζ(i), γ[i+1] > 0 ? " + " : " - ")
    end
    if ν[end] == 0
        UNICODE[2] *= γ_s[end]
        LaTeX[2] *= γ_s[end]
    elseif γ[end] == 1
        UNICODE[2] *= string("ρ", ω[end])
        LaTeX[2] *= string("\\rho", ν |> lastindex |> ζ)
    else
        UNICODE[2] *= string(γ_s[end], "ρ", ω[end])
        LaTeX[2] *= string(γ_s[end], "\\rho", ν |> lastindex |> ζ)
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
    Z_mn = "Z_{$n}^{$(Z.inds.m)}"
    if T <: PartialDerivative
        (; order) = Z
        var = T.parameters[1] == Harmonic ? "θ" : "ρ"
        Z_mn = "\\frac{\\partial^$order $Z_mn}{\\partial $var^$order}"
    end
    Z_Unicode = join(UNICODE, parentheses...)
    Z_LaTeX = latexstring(Z_mn, " = ", join(LaTeX, parentheses...))
    if !iszero(j) && iszero(Z.R.λ[2:end])
        unicode = isint ? Int(N * γ[1]) : "$(Int(γ[1]))" * UNICODE[1]
        latex = replace(unicode, r"√\((\d+)\)" => s"\\sqrt{\1}")
        Z_LaTeX = latexstring(Z_mn, " = ", latex, LaTeX[3])
        Z_Unicode = string(unicode, UNICODE[3])
    end
    return latexstring(Z_mn), Z_LaTeX, Z_Unicode
end

format_strings(m::Int, n::Int) = reverse(format_strings(Z(m, n))[2:3])

function format_strings(recap::Vector)
    W_LaTeX = "ΔW ≈ "
    function ζ(i, sub_index = 0)
        aᵢ = recap[i][:a]
        t = @sprintf("%.3f", sub_index ≠ 1 ? abs(aᵢ) : aᵢ)
        t, "Z_{$(recap[i][:n])}^{$(recap[i][:m])}"
    end
    function η(index, recap)
        for i = 1:length(recap)-1
            W_LaTeX *= string(ζ(index[i], i)..., recap[i+1][:a] > 0 ? " + " : " - ")
        end
    end
    if length(recap) > 8
        sorted_indices = sortperm(recap; by = recapᵢ -> abs(recapᵢ[:a]), rev = true)
        filter!(i -> recap[i][:j] ∉ 0:2, sorted_indices)
        subset_a = getindex(recap, sorted_indices[1:4])
        η(sorted_indices, subset_a)
        W_LaTeX *= string(ζ(sorted_indices[4])...) * "..."
    else
        η(keys(recap), recap)
        W_LaTeX *= string(ζ(lastindex(recap))...)
    end
    return latexstring(W_LaTeX)
end

(Z::AbstractPolynomial)(::Type{String}) = format_strings(Z)[3]
