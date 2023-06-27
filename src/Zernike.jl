# Zernike.jl
# Version 2.0.0
# 2023-6-21
# https://github.com/Sagnac/Zernike.jl

# Generates Zernike polynomials, models wavefront errors, and plots them in Makie.

module Zernike

export Z, W, P, Model, noll_to_j, standardize, fringe_to_j, standardize!

using GLMakie
import .Makie: latexstring, LaTeXString
import Base: show, getindex, iterate

const ∑ = sum

struct Model end

struct RadialPolynomial
    γ::Vector{Float64}
    ν::Vector{Int}
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
include("ScaleAperture.jl")
include("TransformAperture.jl")
include("Docstrings.jl")

function (R::RadialPolynomial)(ρ)
    (; γ, ν) = R
    γ⋅ρ.^ν
end

function (M::Sinusoid)(θ)
    (; m) = M
    m < 0 ? -sin(m * θ) : cos(m * θ)
end

function (Z::Polynomial)(ρ, θ)
    (; N, R, M) = Z
    N * R(ρ) * M(θ)
end

function polar(m::Int, n::Int; scale::Int = 100)
    m = abs(m)
    scale = clamp(scale, 1, 100)
    ϵ₁ = scale * (ceil(Int, π * n) + 1)
    ϵ₁ = clamp(ϵ₁, ϵ₁, 1000)
    ϵ₂ = scale * (2m + 1)
    ϵ₂ = clamp(ϵ₂, ϵ₂, 1000)
    ρ = range(0.0, 1.0, ϵ₁)
    θ = range(0.0, 2π, ϵ₂)
    return ρ, θ
end

# radial order
get_n(j::Int)::Int = ceil((-3 + √(9 + 8j)) / 2)
# azimuthal frequency
get_m(j::Int, n::Int) = 2j - (n + 2)n
# ISO / ANSI / OSA standard single mode-ordering index
get_j(m::Int, n::Int) = ((n + 2)n + m) ÷ 2

function get_mn(j::Int)
    if j < 0
        error("j must be ≥ 0\n")
    end
    n = get_n(j)
    m = get_m(j, n)
    return m, n
end

function noll_to_j(noll::Int)
    if noll < 1
        error("Noll index must be ≥ 1\n")
    end
    n::Int = trunc(√(2noll - 1) + 0.5) - 1
    n_mod_2 = isodd(n)
    m::Int = 2((2noll - (n + 1)n + 1 + n_mod_2) ÷ 4) - n_mod_2
    m = flipsign(m, iseven(noll) ? 1 : -1)
    get_j(m, n)
end

function standardize!(noll::Vector)
    validate_length(noll)
    invpermute!(noll, [noll_to_j(i) + 1 for i in eachindex(noll)])
end

function fringe_to_j(fringe::Int)
    fringe ∉ 1:37 && error("Invalid Fringe index. fringe ∈ 1:37\n")
    if fringe == 37
        return get_j(0, 12)
    end
    d = trunc(√(fringe - 1)) + 1
    d2 = d^2 - fringe
    d2_mod_2 = isodd(d2)
    m::Int = (d2 + d2_mod_2) ÷ 2
    m = flipsign(m, d2_mod_2 ? -1 : 1)
    n::Int = 2(d - 1) - abs(m)
    get_j(m, n)
end

function standardize(fringe::Vector)
    j = fringe_to_j.(eachindex(fringe))
    n_max = get_n(maximum(j))
    j_max = get_j(n_max, n_max)
    a = zeros(eltype(fringe), j_max + 1)
    # normalize
    N = broadcast(x -> √radicand(x...), get_mn.(j))
    a[j.+1] = fringe ./ N
    return a
end

function validate_length(v::Vector)
    len = length(v)
    j_max = len - 1
    n_max = get_n(j_max)
    if j_max ≠ get_j(n_max, n_max)
        error(
            """
            Invalid number of coefficients.
            Coefficients for Zernike polynomials up to m = n_max, n = n_max required.
            """
        )
    end
    return len, n_max
end

function radicand(m::Int, n::Int)
    # Kronecker delta δ_{m0}
    δ(m) = m == 0
    (2n + 2) ÷ (1 + δ(m))
end

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
    (-1)^s * τ[1] / prod(τ[2:4])
end

# computation and construction function
# binds the indices and produces a specific polynomial function
function Zf(m::Int, n::Int)
    μ = abs(m)
    # validate
    if n < 0 || μ > n || isodd(n - μ)
        error("Bounds:\nn ≥ 0\n|m| ≤ n\nn - |m| even\n")
    end
    # upper bound for the sum (number of terms - 1 [indexing from zero])
    k = (n - μ) ÷ 2
    j = get_j(m, n)
    # normalization constant following the orthogonality relation
    N = √radicand(m, n)
    # power (exponent)
    ν = Int[n - 2s for s = 0:k]
    # polynomial coefficients
    if μ == n
        γ = [1.0]
    else
        λ = Φ(μ, n)[end]
        γ = Float64[λ[νᵢ+1] for νᵢ in ν]
    end
    # γ = Float64[λ(μ, n, s, k) for s = 0:k]
    inds = (j = j, n = n, m = m)
    # radial polynomial
    R = RadialPolynomial(γ, ν)
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
    if n ≥ 48
        high_order = true
        println()
        @info "Coefficients are stored in the coeffs field of the current output." γ
    else
        high_order = false
        print("Z = ", Z_Unicode)
    end
    titles = (plot = Z.inds.j < 153 ? Z_LaTeX : Zmn, window = window_title)
    fig = ZPlot(ρ, θ, Zp; high_order, titles...)
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
Z(m::Int, n::Int, ::Model) = Zf(m, n)

Z(; m, n, scale::Int = 100) = Z(m, n; scale)

Z(j::Int; scale::Int = 100) = Z(get_mn(j)...; scale)

Z(j::Int, ::Model) = Zf(get_mn(j)...)

end
