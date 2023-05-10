using Printf

# there's probably an easier and more straightforward way to do this
# but this approach is intuitive enough
# and doesn't require importing another package

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

function format_strings(a::Vector)

    W_LaTeX = "ΔW ≈ "

    function ζ(i, sub_index = 0)
        aᵢ = a[i][:a]
        t = @sprintf "%.3f" (sub_index ≠ 1 ? abs(aᵢ) : aᵢ)
        t, "Z_{$(a[i][:n])}^{$(a[i][:m])}"
    end

    function η(index, a)
        for i = 1:length(a)-1
            W_LaTeX *= string(ζ(index[i], i)..., a[i+1][:a] > 0 ? " + " : " - ")
        end
    end

    if length(a) > 8
        sorted_indices = sortperm(abs.([aᵢ[:a] for aᵢ ∈ a]); rev = true)
        sorted_indices = [i for i ∈ sorted_indices if a[i][:j] ∉ 0:2]
        subset_a = getindex(a, sorted_indices[1:4])
        η(sorted_indices, subset_a)
        W_LaTeX *= string(ζ(sorted_indices[4])...) * "..."
    else
        η(keys(a), a)
        W_LaTeX *= string(ζ(lastindex(a))...)
    end

    return W_LaTeX

end
