using SparseArrays

const Factor = Union{RadialPolynomial, Harmonic}

struct Derivative{T <: Factor} <: AbstractPolynomial
    order::Int
    inds::NamedTuple{(:j, :n, :m), NTuple{3, Int}}
    N::Float64
    R::RadialPolynomial
    M::Harmonic
end

struct Gradient{T <: Polynomial}
    r::Derivative{RadialPolynomial}
    t::Derivative{Harmonic}
    Gradient{Polynomial}(Z::Polynomial) = new(derivatives(Z)...)
end

Gradient(Z::Polynomial) = Gradient{Polynomial}(Z)

(g::Gradient)(ρ::Real, θ::Real = 0) = [g.r(ρ, θ), g.t(ρ, θ) / ρ]

function (g::Gradient)(xy::Complex)
    ρ, θ = polar(xy)
    ∇Z = g(ρ, θ)
    s, c = sc = sincos(θ)
    ∂x = (c, -s) ⋅ ∇Z
    ∂y = sc ⋅ ∇Z
    return [∂x, ∂y]
end

struct Laplacian{T <: Polynomial}
    r1::Derivative{RadialPolynomial}
    r2::Derivative{RadialPolynomial}
    t::Derivative{Harmonic}
    function Laplacian{T}(Z::T) where T <: Polynomial
        new(derivatives(Z)[1], derivatives(Z, 2)...)
    end
end

Laplacian(Z::T) where T <: Polynomial = Laplacian{T}(Z)

function (l::Laplacian)(ρ::Real, θ::Real = 0)
    l.r1(ρ, θ) / ρ + l.r2(ρ, θ) + l.t(ρ, θ) / ρ ^ 2
end

(l::Laplacian)(xy::Complex) = l(polar(xy)...)

derivatives(λ::Vector, order::Int) = spdiagm(1 => 1.0:length(λ)-1) ^ order * λ

function derivatives(Z::Polynomial, order::Int = 1)
    @domain(order > 0, order)
    (; inds, N, R, M) = Z
    (; λ, γ, ν) = R
    (; m) = M
    λ′ = derivatives(λ, order)
    ν′ = Int[νᵢ - order for νᵢ ∈ ν if νᵢ ≥ order]
    γ′ = Float64[λ′[νᵢ+1] for νᵢ ∈ ν′]
    N′ = N * psgn(div(order * (order + sign(m)), 2)) * abs(float(m)) ^ order
    m *= psgn(order)
    R′ = RadialPolynomial(λ′, γ′, ν′)
    M′ = Harmonic(m)
    ∂ρ = Derivative{RadialPolynomial}(order, inds, N, R′, M)
    ∂θ = Derivative{Harmonic}(order, inds, N′, R, M′)
    return ∂ρ, ∂θ
end

function show(io::IO, ∂::T) where T <: Derivative
    print(io, T, " order: ", ∂.order)
end

show(io::IO, ::T) where T <: Union{Gradient, Laplacian} = print(io, T)

function Wavefront(g::Gradient)
    (; m, n) = g.r.inds
    Wavefront(Gradient, m, n)
end

function Wavefront(::Type{<:Gradient}, m::Int, n::Int)
    cx, cy = grad(m, n)
    order = conjugate_indices(n - 1)
    ax = to_real(cx, order, 1)
    ay = to_real(cy, order, 1)
    ∂x = Wavefront(ax)
    ∂y = Wavefront(ay)
    return ∂x, ∂y
end

function conjugate_indices(n_max::Int)
    order = Vector{NTuple{3}{Int}}(undef, get_j(max(n_max, 0)) + 1)
    for i in eachindex(order)
        m, n = get_mn(i - 1)
        order[i] = (i, get_j(-m, n) + 1, sign(m))
    end
    return order
end

#= The following algorithm computes Cartesian derivatives in terms of Zernike polynomials. It is based on mathematical formulas found in:

https://doi.org/10.1364/JOSAA.31.001604

A. J. E. M. Janssen, "Zernike expansion of derivatives and Laplacians of the Zernike circle polynomials," J. Opt. Soc. Am. A 31, 1604-1613 (2014)

https://opg.optica.org/josaa/abstract.cfm?uri=josaa-31-7-1604

Note:

Complex polynomials correspond to a linear combination of two real, standard
Zernike polynomials. As such, derivatives must be computed for conjugate indices
and appropriately combined before converting the complex coefficients to
the real ones.

=#

function grad(m::Int, n::Int)
    len = get_j(max(n - 1, 0)) + 1
    c = [zeros(ComplexF64, len) for i = 1:4]
    # c[1:2] ≡ ∂/∂x (Z(±m, n)), c[3:4] ≡ ∂/∂y (Z(±m, n))
    μ = abs(m)
    for ci = n:-2:μ
        n′ = ci - 1
        n′ < 0 && break
        for (i, m′) ∈ enumerate(((μ + 1, μ - 1), (-μ + 1, -μ - 1)))
            for (i′, m′′) ∈ enumerate(m′)
                if abs(m′′) ≤ n′
                    i′′ = get_j(m′′, n′) + 1
                    c[i][i′′] = ci
                    c[i+2][i′′] = psgn(i′) * im * ci
                end
            end
        end
    end
    cx, cy = to_complex(c, m)
    return cx, cy, c
end

function to_complex(c::Vector{Vector{ComplexF64}}, m::Int)
    if m < 0
        cx, cy = (-(c[i], c[i+1]) / 2im for i = (1, 3))
    else
        cx, cy = (+(c[i], c[i+1]) / 2   for i = (1, 3))
    end
    return cx, cy
end
