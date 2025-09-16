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

function Wavefront(g::Gradient)
    (; m, n) = g.r.inds
    index_order_1 = NTuple{3}{Int}[]
    index_order_2 = NTuple{3}{Int}[]
    mn_orders_1 = NTuple{2}{Int}[]
    mn_orders_2 = NTuple{2}{Int}[]
    c1 = Float64[]
    c2 = Float64[]
    m1 = m - 1
    m2 = m + 1
    for s = 0:div(n - abs(m), 2)
        cs = n - 2s
        n′ = cs - 1
        n′ < 0 && break
        if abs(m1) ≤ n′
            set_order!(index_order_1, m1, n′)
            push!(mn_orders_1, (m1, n′))
            push!(c1, cs)
        end
        if abs(m2) ≤ n′
            set_order!(index_order_2, m2, n′)
            push!(mn_orders_2, (m2, n′))
            push!(c2, cs)
        end
    end
    if isempty(c1)
        c1 = [0.0]
        index_order_1 = [(1, 1, 0)]
        mn_orders_1 = [(0, 0)]
    end
    if isempty(c2)
        c2 = [0.0]
        index_order_2 = [(1, 1, 0)]
        mn_orders_2 = [(0, 0)]
    end
    c1 = standardize(c1, mn_orders_1)
    c2 = standardize(c2, mn_orders_2)
    c1 = to_complex(c1, index_order_1)
    c2 = to_complex(c2, index_order_2)
    cx = c1 + c2
    cy = im * (c1 - c2)
    ax = to_real(cx, index_order_1)
    ay = to_real(cy, index_order_2)
    ∂Z_∂x = Wavefront(ax)
    ∂Z_∂y = Wavefront(ay)
    return (; ∂Z_∂x, ∂Z_∂y)
end
