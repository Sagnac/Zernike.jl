#= This algorithm is based on a mathematical formulation presented in:

https://doi.org/10.1364/JOSAA.24.000569

Linda Lundström and Peter Unsbo, "Transformation of Zernike coefficients: scaled, translated, and rotated wavefronts with circular and elliptical pupils," J. Opt. Soc. Am. A 24, 569-577 (2007)

https://opg.optica.org/josaa/abstract.cfm?uri=josaa-24-3-569

=#

using LinearAlgebra

const b = binomial
const ei = cis # ei(x) = exp(im*x)

function S(v::Vector{T}, ε::T, δ::Complex{T}, ϕ::T, ω::Tuple{T,T}) where T <: Float64
    !(0.0 ≤ ε ≤ 1.0) && error("Bounds: 0.0 ≤ ε ≤ 1.0\n")
    !(0.0 ≤ ε + abs(δ) ≤ 1.0) && error("Bounds: 0.0 ≤ ε + |δ| ≤ 1.0\n")
    !(0.0 < ω[1] ≤ 1.0) && error("Bounds: 0.0 < ω[1] ≤ 1.0\n")
    len, n_max = validate_length(v)
    remap = Dict{Tuple{Int, Int}, Int}()
    order = Tuple{Int, Int, Int}[]
    # normalization factors for complex Zernike polynomials
    N = zeros(Float64, len, len)
    # radial coefficient block-diagonal matrix
    R = zeros(Float64, len, len)
    λ = Φ(n_max, n_max)
    _scale_ = iszero(δ)
    _translate_ = !_scale_
    _rotate_ = !iszero(ϕ)
    _elliptic_ = !isone(ω[1])
    if _scale_
        η_s = zeros(Float64, len, len)
    end
    η_r = _rotate_ ? zeros(ComplexF64, len, len) : I
    i = 0
    for m = -n_max:n_max
        μ = abs(m)
        for n = μ:2:n_max
            i += 1
            k1 = (n - μ) ÷ 2
            push!(remap, (m, n) => i)
            push!(order, (get_j(m, n) + 1, get_j(-m, n) + 1, sign(m)))
            N[i,i] = √(n+1)
            γ = λ[get_i(μ, n)]
            for s = 0:k1
                setindex!(R, γ[n-2s+1], remap[(m, n-2s)], i)
            end
            if _rotate_
                η_r[i,i] = ei(m * ϕ)
            end
            if _scale_
                η_s[i,i] = ε ^ n
            end
        end
    end
    if _translate_ && _elliptic_
        η_s, η_e = translate_ellipse(ε, δ, ω..., remap)
    elseif _translate_
        η_s = translate(ε, δ, remap)
        η_e = I
    elseif _elliptic_
        η_e = elliptical(ω..., remap)
    else
        η_e = I
    end
    # net transformation matrix
    η = η_e * η_r * η_s
    # complex Zernike expansion coefficients
    c = to_complex(v, order)
    # conversion matrix
    C = (R * N) \ (η * R * N)
    # transformed set of complex Zernike expansion coefficients
    c′ = C * c
    # transformed set of standard Zernike expansion coefficients
    v2 = to_real(c′, order)
    return v2, n_max
end

macro init_kernel!()
    quote
        k2 = (n + m) ÷ 2
        k3 = (n - m) ÷ 2
        i = remap[(m, n)]
    end |> esc
end

macro translation_kernel!()
    quote
        n′ = n - p - q
        m′ = m - p + q
        t = b(k2,p) * b(k3,q) * ε^n′ * ρ^(p+q) * ei((p-q)θ)
        η_s[remap[(m′, n′)], i] += t
    end |> esc
end

macro elliptical_kernel!()
    quote
        m′ = m - 2p + 2q
        t = 0.5^n * b(k2,p) * b(k3,q) * (ξ+1)^(n-p-q) * (ξ-1)^(p+q) * ei(2(p-q)φ)
        η_e[remap[(m′, n)], i] += t
    end |> esc
end

function translate(ε::Float64, δ::ComplexF64, remap::Dict)
    ρ, θ = polar(δ)
    len = length(remap)
    n_max = get_n(len - 1)
    η_s = zeros(ComplexF64, len, len)
    for m = -n_max:n_max, n = abs(m):2:n_max
        @init_kernel!
        for p = 0:k2, q = 0:k3
            @translation_kernel!
        end
    end
    return η_s
end

function elliptical(ξ::Float64, φ::Float64, remap::Dict)
    len = length(remap)
    n_max = get_n(len - 1)
    η_e = zeros(ComplexF64, len, len)
    for m = -n_max:n_max, n = abs(m):2:n_max
        @init_kernel!
        for p = 0:k2, q = 0:k3
            @elliptical_kernel!
        end
    end
    return η_e
end

function translate_ellipse(ε::Float64, δ::ComplexF64, ξ::Float64, φ::Float64,
                           remap::Dict)
    ρ, θ = polar(δ)
    len = length(remap)
    n_max = get_n(len - 1)
    η_s = zeros(ComplexF64, len, len)
    η_e = zeros(ComplexF64, len, len)
    for m = -n_max:n_max, n = abs(m):2:n_max
        @init_kernel!
        for p = 0:k2, q = 0:k3
            @translation_kernel!
            @elliptical_kernel!
        end
    end
    return η_s, η_e
end

function to_complex(v::Vector{Float64}, order::Vector{Tuple{Int, Int, Int}})
    c = similar(v, ComplexF64)
    for (i, j) in pairs(order)
        if j[3] < 0
            c[i] = complex(v[j[2]], v[j[1]]) / √2
        elseif j[3] > 0
            c[i] = complex(v[j[1]], -v[j[2]]) / √2
        else
            c[i] = v[j[1]]
        end
    end
    return c
end

function to_real(c::Vector{Complex{Float64}}, order::Vector{Tuple{Int, Int, Int}})
    c2 = c[sortperm(order; by = first)]
    v2 = similar(c2, Float64)
    for (i, j) in pairs(order)
        if j[3] < 0
            v2[j[1]] = (c2[j[1]] - c2[j[2]])  / √2 |> imag
        elseif j[3] > 0
            v2[j[1]] = (c2[j[1]] + c2[j[2]]) / √2 |> real
        else
            v2[j[1]] = c2[j[1]] |> real
        end
    end
    return v2
end

function P(
    v::Vector{T},
    ε::T,
    δ::Complex{T} = 0.0im,
    ϕ::T = 0.0,
    ω::Tuple{T, T} = (1.0, 0.0);
    precision::Int = precision, finesse::Int = wavefront_finesse
) where T <: Float64
    v2, n_max = S(v, ε, δ, ϕ, ω)
    Zᵢ = similar(v, Polynomial)
    ΔW = Ψ(v2, Zᵢ, n_max; precision)
    Λ(ΔW; finesse)
end

function P(v::Vector{T}, ε::T, δ::Complex{T}, ϕ::T, ω::Tuple{T, T}, ::Type{Model};
           precision::Int = precision) where T <: Float64
    v2, n_max = S(v, ε, δ, ϕ, ω)
    Zᵢ = similar(v, Polynomial)
    Ψ(v2, Zᵢ, n_max; precision)
end

function P(v::Vector{T}, ε::T, ::Type{Model};
           precision::Int = precision) where T <: Float64
    P(v, ε, 0.0im, zero(T), (1.0, 0.0), Model; precision)
end

function P(v::Vector{T}, ε::T, δ::Complex{T}, ::Type{Model};
           precision::Int = precision) where T <: Float64
    P(v, ε, δ, zero(T), (1.0, 0.0), Model; precision)
end

function P(v::Vector{T}, ε::T, δ::Complex{T}, ϕ::T, ::Type{Model};
           precision::Int = precision) where T <: Float64
    P(v, ε, δ, ϕ, (1.0, 0.0), Model; precision)
end
