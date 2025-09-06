using SparseArrays

const Factor = Union{RadialPolynomial, Harmonic}

struct PartialDerivative{T <: Factor} <: AbstractPolynomial
    order::Int
    inds::NamedTuple{(:j, :n, :m), NTuple{3, Int}}
    N::Float64
    R::RadialPolynomial
    M::Harmonic
end

struct Gradient{T <: Polynomial}
    r::PartialDerivative{RadialPolynomial}
    t::PartialDerivative{Harmonic}
    Gradient{Polynomial}(Z::Polynomial) = new(derivatives(Z)...)
end

Gradient(Z::Polynomial) = Gradient{Polynomial}(Z)

(g::Gradient)(ρ::Real, θ::Real = 0) = [g.r(ρ, θ), g.t(ρ, θ) / ρ]

function derivatives(Z::Polynomial, order::Int = 1)
    @domain(order > 0, order)
    (; inds, N, R, M) = Z
    (; λ, γ, ν) = R
    (; m) = M
    D = spdiagm(1 => 1.0:length(λ)-1)
    λ′ = D ^ order * λ
    ν′ = Int[νᵢ - order for νᵢ ∈ ν if νᵢ ≥ order]
    γ′ = Float64[λ′[νᵢ+1] for νᵢ ∈ ν′]
    N′ = N * (-1) ^ div(order * (order + sign(m)), 2) * abs(float(m)) ^ order
    m *= (-1) ^ order
    R′ = RadialPolynomial(λ′, γ′, ν′)
    M′ = Harmonic(m)
    ∂Z_∂ρ = PartialDerivative{RadialPolynomial}(order, inds, N, R′, M)
    ∂Z_∂θ = PartialDerivative{Harmonic}(order, inds, N′, R, M′)
    return ∂Z_∂ρ, ∂Z_∂θ
end

function show(io::IO, ∂::T) where T <: PartialDerivative
    print(io, T, " order: ", ∂.order)
end

show(io::IO, ::T) where T <: Gradient = print(io, T)
