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

function complex(g::Gradient)
    function(xy::Complex)
        abs(xy) > 1.0 && return 0.0im
        complex(g(xy)...)
    end
end

grad(Z::Polynomial; kw...) = Wavefront(Gradient, Z.inds.m, Z.inds.n; kw...)

function Wavefront(g::Gradient; kw...)
    (; m, n) = g.r.inds
    Wavefront(Gradient, m, n; kw...)
end

function Wavefront(::Type{<:Gradient}, m::Int, n::Int; normalize::Bool = true)
    c = grad(m, n)
    cx, cy = to_complex(c, m)
    order = conjugate_indices(n - 1)
    ax = to_real(cx, order, 1)
    ay = to_real(cy, order, 1)
    if normalize
        N = √N²(m, n)
        ax *= N
        ay *= N
    end
    ∂x = Wavefront(ax)
    ∂y = Wavefront(ay)
    return ∂x, ∂y
end

function Wavefront(::Type{<:Gradient}, j::Int; kw...)
    @domain_check_j
    Wavefront(Gradient, get_mn(j)...; kw...)
end

function derivatives(W::Wavefront)
    [mapreduce(*, +, ∂, W.a) for ∂ ∈ eachrow(stack(grad(Z) for Z ∈ W.Z))]
end

function Wavefront(l::Laplacian; kw...)
    (; m, n) = l.r1.inds
    Wavefront(Laplacian, m, n; kw...)
end

function Wavefront(::Type{<:Laplacian}, m::Int, n::Int; normalize::Bool = true)
    a = lap(m, n)
    normalize && (a *= √N²(m, n))
    return Wavefront(a)
end

function Wavefront(::Type{<:Laplacian}, j::Int; kw...)
    @domain_check_j
    Wavefront(Laplacian, get_mn(j)...; kw...)
end

function conjugate_indices(n_max::Int)
    order = Vector{NTuple{3}{Int}}(undef, to_i(n_max))
    for i in eachindex(order)
        m, n = get_mn(i - 1)
        order[i] = (i, to_i(-m, n), sign(m))
    end
    return order
end

to_i(n_max::Int) = get_j(max(n_max, 0)) + 1

#= The following algorithms compute Cartesian derivatives in terms of Zernike polynomials. It is based on mathematical formulas found in:

https://doi.org/10.1364/JOSAA.31.001604

A. J. E. M. Janssen, "Zernike expansion of derivatives and Laplacians of the Zernike circle polynomials," J. Opt. Soc. Am. A 31, 1604-1613 (2014)

https://opg.optica.org/josaa/abstract.cfm?uri=josaa-31-7-1604

Note:

Complex polynomials correspond to a linear combination of two real, standard
Zernike polynomials. As such, derivatives must be computed for conjugate indices
and appropriately combined before converting the complex coefficients to
the real ones. For second order derivatives this is not necessary, hence the
Laplacian can be computed directly.

=#

function grad(m::Int, n::Int)
    μ = abs(m)
    @domain_check_mn
    len = to_i(n - 1)
    c = zeros(ComplexF64, (len, 4))
    # c[:,1:2] ≡ ∂/∂x (Z(±|m|, n)), c[:,3:4] ≡ ∂/∂y (Z(±|m|, n))
    t = (μ + 1, μ - 1)
    t′ = (t, .-reverse(t))
    for s = n:-2:μ
        n′ = s - 1
        n′ < 0 && break
        for (i, m′) ∈ enumerate(t′), (i′, m′′) ∈ enumerate(m′)
            abs(m′′) > n′ && continue
            i′′ = to_i(m′′, n′)
            c[i′′,i] = s
            c[i′′,i+2] = psgn(i′) * im * s
        end
    end
    return c
end

function grad(m::Int, n::Int, ::Type{Vector})
    c = grad(m, n)
    m > 0 ? (return c[:,1], c[:,3]) : (return c[:,2], c[:,4])
end

function to_complex(c::Matrix{ComplexF64}, m::Int)
    c = eachcol(c)
    if m < 0
        cx, cy = (-(c[i], c[i+1]) / 2im for i = (1, 3))
    else
        cx, cy = (+(c[i], c[i+1]) / 2   for i = (1, 3))
    end
    return cx, cy
end

function lap(m::Int, n::Int)
    μ = abs(m)
    @domain_check_mn
    a = zeros(to_i(n - 2))
    for s = μ:2:n-2
        a[to_i(m, s)] = (s + 1) * (n + s + 2) * (n - s)
    end
    return a
end

function W(B_plus::Vector{ComplexF64}, B_minus::Vector{ComplexF64})
    _, n_max = validate_length(B_plus)
    n′_max = n_max + 1
    α = Vector{ComplexF64}(undef, to_i(n′_max))
    for i in eachindex(α)
        m, n = get_mn(i - 1)
        m′ = (m + 1, m - 1)
        n′ = (n - 1, n + 1, n + 2)
        C = [compute_C(m, ni) for ni ∈ (n, n′[3])]
        φ1 = (_get(B_plus, m′[1], n′[1]) + _get(B_minus, m′[2], n′[1])) / 2
        φ2 = (_get(B_plus, m′[1], n′[2]) + _get(B_minus, m′[2], n′[2])) / 2
        α[i] = C[1] * φ1 - C[2] * φ2
    end
    α[1] = 0.0im
    a = to_real(α, conjugate_indices(n′_max), 1)
    return a
end

compute_B(m::Int, n::Int) = compute_B(grad(m, n, Vector)...)

function compute_B(cx::Vector{ComplexF64}, cy::Vector{ComplexF64})
    B_plus, B_minus = (op(cx, im * cy) for op ∈ (+, -))
    return B_plus, B_minus
end

function _get(B::Vector{ComplexF64}, m::Int, n::Int)
    try
        B[to_i(m, n)]
    catch
        0.0im
    end
end

function compute_C(m::Int, n::Int)::ComplexF64
    μ = abs(m)
    μ == n ? inv(μ) : inv(2n)
end

function W(∂x::Vector{Float64}, ∂y::Vector{Float64}; normalize::Bool = true)
    @assert length(∂x) == length(∂y)
    _, n_max = validate_length(∂x)
    order = conjugate_indices(n_max)
    cx, cy = (to_complex(∂, order, 2) for ∂ in (∂x, ∂y))
    a = W(compute_B(cx, cy)...)
    normalize && (a ./= N(0:length(a)-1))
    return a
end

W(∂x::Wavefront, ∂y::Wavefront; normalize::Bool = true) = W(∂x[], ∂y[]; normalize)

function Wavefront(::Type{<:Derivative},
                   ∂x::Vector{Float64}, ∂y::Vector{Float64};
                   normalize::Bool = true)
    Wavefront(W(∂x, ∂y; normalize))
end
