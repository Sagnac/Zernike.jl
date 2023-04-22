# Zernike.jl
# Version 0.1.0
# 2023-4-22
# https://github.com/Sagnac/Zernike

# Generates Zernike Polynomials and plots them in Makie

# Julia v1.8.0

using GLMakie # v0.7.3

function Z(m::Integer, n::Integer)

    μ = abs(m)

    if n < 0 || μ > n || isodd(n - μ)
        @error "Bounds:\nn ≥ 0\n|m| ≤ n\nn - |m| even"
        return
    end

    # number of terms -1 (zero-based initial index)
    k::Int = (n - μ) / 2
    # ISO / ANSI / OSA standard single mode-ordering index
    j::Int = (n * (n + 2) + m) / 2
    # Kronecker delta δₘ₀
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
            (n+μ)/2 - s
        ]
        τ = t |> fact
        if j > 278
            τ = convert(Vector{BigInt}, τ)
        end
        (-1)^s * τ[1] / prod(τ[2:4])
    end

    γ = Integer[λ(s) for s = 0:k]

    # power
    σ = Integer[n - 2s for s = 0:k]

    # radial polynomial
    R(ρ)::Float64 = sum(γ[i] * ρ ^ σ[i] for i = 1:k+1)

    # azimuthal / meridional component
    M(θ) = m < 0 ? -sin(m * θ) : cos(m * θ)

    Z(ρ,θ) = N * R(ρ) * M(θ)

    ρ = range(0, 1, 100)
    θ = range(0, 2π, 100)
    ρᵪ = [ρᵢ * cos(θᵢ) for θᵢ ∈ θ, ρᵢ ∈ ρ]
    ρᵧ = [ρⱼ * sin(θⱼ) for θⱼ ∈ θ, ρⱼ ∈ ρ]
    Zp = Z.(ρ',θ)

    # string formatting

    # there's probably an easier and more straightforward way to do this
    # but this approach is intuitive enough
    # and doesn't require importing another package

    # legend:

    # for Unicode formatting (for printing):
    # Γ = final polynomial
    # η = prefactor
    # Φ = polynomial terms
    # φ = angular term

    # for LaTeX formatting (for the plot title):
    # Λ = final polynomial
    # ϵ = prefactor
    # Ξ = polynomial terms
    # ψ = angular term

    if isinteger(N)
        η = ϵ = N == 1 ? "" : string(Int(N))
    else
        η = "√($N²)"
        ϵ = "\\sqrt{$N²}"
    end

    if m == 0
        φ = ""
        ψ = ""
    elseif m == -1
        φ = "sin(θ)"
        ψ = "\\sin(\\theta)"
    elseif m == 1
        φ = "cos(θ)"
        ψ = "\\cos(\\theta)"
    elseif m < -1
        φ = "sin($(μ)θ)"
        ψ = "\\sin($(μ)\\theta)"
    else
        φ = "cos($(m)θ)"
        ψ = "\\cos($(m)\\theta)"
    end

    if k != 0
        η *= "("
        ϵ *= "("
        φ = ")" * φ
        ψ = ")" * ψ
    end

    superscripts = ['\u2070'; '\u00B9'; '\u00B2'; '\u00B3'; '\u2074':'\u2079']

    ω = fill("", length(σ))
    
    for i in eachindex(σ), j in σ[i] |> digits |> reverse
        ω[i] *= superscripts[j+1]
    end

    for (i, v) in pairs(ω)
        if v == string(superscripts[2])
            ω[i] = ""
        end
    end

    Φ = ""
    Ξ = ""

    for i = 1:k
        Φ *= string(abs(γ[i]), "ρ", ω[i], γ[i+1] > 0 ? " + " : " \u2212 ")
        Ξ *= string(abs(γ[i]), "\\rho", "^{$(σ[i])}", γ[i+1] > 0 ? " + " : " - ")
    end

    if σ[end] == 0
        Φ *= string(abs(γ[end]))
        Ξ *= string(abs(γ[end]))
    elseif γ[end] == 1
        Φ *= string("ρ", ω[end])
        Ξ *= string("\\rho", σ[end] == 1 ? "" : "^{$(σ[end])}")
    else
        Φ *= string(abs(γ[end]), "ρ", ω[end])
        Ξ *= string(abs(γ[end]), "\\rho", σ[end] == 1 ? "" : "^{$(σ[end])}")
    end

    Γ = η * Φ * φ
    Λ = ϵ * Ξ * ψ

    # plot

    resolution = (1400, 1000)
    textsize = 35

    axis3attributes = (
        title = L"Z = %$Λ",
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
    # this toggle allows you to interactively change the perspective

    zproperties = (
        :zlabelvisible,
        :zgridvisible,
        :zticksvisible,
        :zticklabelsvisible,
        :zspinesvisible,
    )

    Label(fig[1,2], "Pupil view";
          tellheight = false, valign = :bottom, textsize = 0.76textsize)
    pupil = Toggle(fig[1,3]; tellheight = false, valign = :bottom)

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

    println("j = $j, m = $m, n = $n")
    println("Z = ", Γ)

    set_window_config!(title = "Zernike Polynomial", focus_on_show = true)
    display(fig)

    return

end

# methods

Z(; m, n) = Z(m, n)

function Z(j::Integer)
    # radial order
    n::Integer = ceil((-3 + sqrt(9 + 8j)) / 2)
    # azimuthal frequency
    m::Integer = 2j - n * (n + 2)
    Z(m, n)
end

return
