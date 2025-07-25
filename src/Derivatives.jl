struct PartialDerivative{T <: Union{RadialPolynomial, Harmonic}} <: AbstractPolynomial
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

"""
    Zernike.Gradient(Z::Polynomial)

Returns ∇Z(ρ, θ).
"""
Gradient(Z::Polynomial) = Gradient{Polynomial}(Z)

(g::Gradient)(ρ::Real, θ::Real = 0) = [g.r(ρ, θ), g.t(ρ, θ) / ρ]

p(i, order) = ∏(i-order+1:float(i))

s(m, order)::Int = order % 4 ∈ (0, 1+2(m > 0)) || -1

# TODO: this needs unit tests
"""
    Zernike.derivatives(Z::Polynomial, order::Int = 1)

Computes the nth order partial derivatives of Z(ρ, θ).
"""
function derivatives(Z::Polynomial, order::Int = 1)
    order < 1 && throw(DomainError(order))
    (; inds, N, R, M) = Z
    (; λ, γ, ν) = R
    (; m) = M
    λ′ = zeros(length(λ))
    i = eachindex(λ) .- 1
    j = @. i ≥ order
    @. λ′[j] = λ[j] * p(i[j], order)
    λ′ = shift(λ′, -order)
    ν′ = Int[νᵢ - order for νᵢ ∈ ν if νᵢ ≥ order]
    γ′ = Float64[λ′[νᵢ+1] for νᵢ ∈ ν′]
    N′ = N * s(m, order) * abs(m) ^ order
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
