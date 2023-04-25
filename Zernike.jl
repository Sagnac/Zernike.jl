# Zernike.jl
# Version 0.1.2
# 2023-4-25
# https://github.com/Sagnac/Zernike

# Generates Zernike Polynomials and plots them in Makie

# Julia v1.8.0

using GLMakie # v0.7.3

function Z(m::Integer, n::Integer)

    μ::Integer = abs(m)

    if n < 0 || μ > n || isodd(n - μ)
        @error "Bounds:\nn ≥ 0\n|m| ≤ n\nn - |m| even"
        return
    end

    # upper bound for the sum (number of terms -1 [indexing from zero])
    k::Integer = (n - μ) / 2
    # ISO / ANSI / OSA standard single mode-ordering index
    j::Integer = ((n + 2)n + m) / 2
    # Kronecker delta δ_{m0}
    δ(m) = m == 0
    # radicand
    N²::Int = (2n + 2) / (1 + δ(m))
    # normalization constant following the orthogonality relation
    N = √N²

    function fact(t)
        bound = Int ≡ Int32 ? 12 : 20
        [i > bound ? (factorial ∘ big)(i) : factorial(i) for i in t]
    end

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
        (-1)^s * τ[1] / prod(τ[2:4])
    end

    γ = Integer[λ(s) for s = 0:k]

    # power (exponent)
    ν = Integer[n - 2s for s = 0:k]

    # radial polynomial
    R(ρ)::Float64 = sum(γ[i] * ρ ^ ν[i] for i = 1:k+1)

    # azimuthal / meridional component
    M(θ) = m < 0 ? -sin(m * θ) : cos(m * θ)

    Z(ρ,θ) = N * R(ρ) * M(θ)

    ρ = range(0, 1, 100)
    θ = range(0, 2π, 100)
    ρᵪ = [ρⱼ * cos(θᵢ) for θᵢ ∈ θ, ρⱼ ∈ ρ]
    ρᵧ = [ρⱼ * sin(θᵢ) for θᵢ ∈ θ, ρⱼ ∈ ρ]
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
    elseif N != 1
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

    parentheses = k != 0 ? ("(", ")") : ""

    Z_Unicode = join(unicode, parentheses...)
    Z_LaTeX = "Z_{$n}^{$m} = " * join(latex, parentheses...)

    # plot

    resolution = (1400, 1000)
    textsize = 35

    axis3attributes = (
        title = L"%$Z_LaTeX",
        titlesize = textsize,
        xlabel = L"\rho_x",
        ylabel = L"\rho_y",
        zlabel = L"Z",
        xlabelsize = textsize,
        ylabelsize = textsize,
        zlabelsize = textsize,
        azimuth = 0.275π,
        protrusions = 80,
    )

    fig = Figure(; resolution)
    axis3 = Axis3(fig[1,1]; axis3attributes...)
    surface!(axis3, ρᵪ, ρᵧ, Zp; colormap = :oslo)

    # hacky way to produce a top-down heatmap-style view without generating
    # another plot with a different set of data
    # accomplished by adding a toggle which changes the perspective on demand

    zproperties = (
        :zlabelvisible,
        :zgridvisible,
        :zticksvisible,
        :zticklabelsvisible,
        :zspinesvisible,
    )

    label = Label(fig, "Pupil view"; textsize = 0.76textsize)
    pupil = Toggle(fig)
    fig[1,2] = grid!([label pupil]; tellheight = false, valign = :bottom)

    on(pupil.active) do active
        for property in zproperties
            getfield(axis3, property)[] = !active
        end
        axis3.azimuth = active ? -π/2 : 0.275π
        axis3.elevation = active ? π/2 : π/8
        axis3.ylabeloffset = active ? 90 : 40
        axis3.xlabeloffset = active ? 50 : 40
    end

    # print and display

    println(indices)
    println("Z = ", Z_Unicode)

    set_window_config!(title = "Zernike Polynomial: $indices", focus_on_show = true)
    display(fig)

    return

end

# methods

Z(; m, n) = Z(m, n)

function Z(j::Integer)
    # radial order
    n::Integer = ceil((-3 + sqrt(9 + 8j)) / 2)
    # azimuthal frequency
    m::Integer = 2j - (n + 2)n
    Z(m, n)
end

return
