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
    ∂Z_∂x = (c, -s) ⋅ ∇Z
    ∂Z_∂y = sc ⋅ ∇Z
    return ∂Z_∂x, ∂Z_∂y
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
    ∂Z_∂ρ = Derivative{RadialPolynomial}(order, inds, N, R′, M)
    ∂Z_∂θ = Derivative{Harmonic}(order, inds, N′, R, M′)
    return ∂Z_∂ρ, ∂Z_∂θ
end

function show(io::IO, ∂::T) where T <: Derivative
    print(io, T, " order: ", ∂.order)
end

show(io::IO, ::T) where T <: Union{Gradient, Laplacian} = print(io, T)
