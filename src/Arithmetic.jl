import Base: +, -, *, /, ^, sum, prod, factorial, convert, promote_rule

struct MixedPhase{S, T <: WavefrontError} <: Phase
    b::Float64
    W::Vector{T}
end

Superposition{T} = MixedPhase{:sum, T}

Superposition(b, W) = Superposition{WavefrontError}(b, W)

Superposition(W) = Superposition(1.0, W)

Product{T} = MixedPhase{:product, T}

Product(b, W) = Product{WavefrontError}(b, W)

Product(W) = Product(1.0, W)

function (ΔW::Superposition{<:WavefrontError})(ρ, θ)
    ΔW.b * ∑(W(ρ, θ) for W ∈ ΔW.W; init = 0.0)
end

∑(v::Vector{<:Phase}) = Superposition(v)

(ΔW::Product{<:WavefrontError})(ρ, θ) = ΔW.b * ∏(W(ρ, θ) for W ∈ ΔW.W; init = 1.0)

∏(v::Vector{<:Phase}) = Product(v)

factorial(z::Polynomial) = Product([Z(j) for j = 0:z.inds.j])

φ1::Phase + φ2::Phase = add_subtract(+, φ1, φ2)
φ1::Phase - φ2::Phase = add_subtract(-, φ1, φ2)

α::Real * φ::Phase = multiply(float(α), φ)
φ::Phase * α::Real = α * φ
φ::Phase / α::Real = inv(α) * φ
φ1::Phase * φ2::Phase = multiply(φ1, φ2)
φ::Phase ^ n::Integer = exponentiate(φ, n)

φ::Phase + α::Real = φ + α * piston
α::Real + φ::Phase = φ + α

+(φ::Phase) = φ
-(φ::Phase) = -1 * φ

const converted = [(-1, -1)]

function convert(::Type{<:WavefrontError}, Z::Polynomial)
    a1 = 1.0
    a = [a1]
    (; m, n, j) = Z.inds
    recap = [(; j, n, m, a = a1)]
    v = standardize(a, [(m, n)])
    fit_to = converted
    WavefrontError(recap, v, n, fit_to, a, [Z], precision)
end

function convert(::Type{<:Superposition}, Z::Polynomial)
    Superposition([convert(WavefrontError, Z)])
end

function convert(::Type{<:Superposition}, W::WavefrontError)
    Superposition([W])
end

function convert(::Type{<:WavefrontError}, ΔW::Product)
    OPD = ΔW.(polar()...)
    n_max, precision = params(ΔW)
    W(OPD, n_max; precision)
end

promote_rule(::Type{<:WavefrontError}, ::Type{Polynomial}) = WavefrontError
promote_rule(::Type{<:Superposition}, ::Type{Polynomial}) = Superposition
promote_rule(::Type{<:Superposition}, ::Type{<:WavefrontError}) = Superposition
promote_rule(::Type{Polynomial}, ::Type{<:Product}) = WavefrontError
promote_rule(::Type{<:WavefrontError}, ::Type{<:Product}) = WavefrontError

function add_subtract(f, φ1::Phase, φ2::Phase)
    add_subtract(f, promote(φ1, φ2)...)
end

function add_subtract(f, W1::T, W2::T) where T <: WavefrontError
    v1 = W1.v
    v2 = W2.v
    len1 = length(v1)
    len2 = length(v2)
    if len2 > len1
        v1 = [v1; zeros(len2 - len1)]
    elseif len1 > len2
        v2 = [v2; zeros(len1 - len2)]
    end
    v3 = f(v1, v2)
    Z3 = similar(v3, Polynomial)
    foreach(Z -> setindex!(Z3, Z, Z.inds.j + 1), W1.Z ∪ W2.Z)
    n_max = max(W1.n_max, W2.n_max)
    orders = setdiff(W1.fit_to ∪ W2.fit_to, converted)
    precision = max(W1.precision, W2.precision)
    Ψ(v3, Z3, n_max, orders; precision)
end

function add_subtract(f, Z1::T, Z2::T) where T <: Polynomial
    Z1.inds.j == Z2.inds.j && return 2 * Z1
    orders = [(Z.inds.m, Z.inds.n) for Z ∈ (Z1, Z2)]
    a = [1.0, f(1.0)]
    WavefrontError(orders, a)
end

function add_subtract(f, ΔW1::T, ΔW2::T) where T <: Superposition
    Superposition([ΔW1.b * ΔW1.W; f(ΔW2.b) * ΔW2.W])
end

function add_subtract(f, ΔW1::T, ΔW2::T) where T <: Product
    wavefront_expansion(f, ΔW1, ΔW2)
end

add_subtract(f, ΔW1::MixedPhase, ΔW2::MixedPhase) = wavefront_expansion(f, ΔW1, ΔW2)

function multiply(α::Float64, W::WavefrontError)
    (; v) = W
    Z = similar(v, Polynomial)
    foreach(Zᵢ -> setindex!(Z, Zᵢ, Zᵢ.inds.j + 1), W.Z)
    Ψ(α * v, Z, W.n_max, W.fit_to; W.precision)
end

multiply(α::Float64, Z::Polynomial) = WavefrontError([(Z.inds.m, Z.inds.n)], [α])

multiply(b::Float64, ΔW::T) where T <: MixedPhase = T(b, ΔW.W)

multiply(φ1::Phase, φ2::Phase) = multiply(promote(φ1, φ2)...)

function multiply(φ1::T, φ2::T) where T <: Phase
    φ1::WavefrontError = φ1
    φ2::WavefrontError = φ2
    ρ, θ = polar()
    φ3 = @. φ1(ρ, θ) * φ2(ρ, θ)
    n_max = φ1.n_max + φ2.n_max
    precision = max(φ1.precision, φ2.precision)
    W(φ3, n_max; precision)
end

multiply(ΔW1::MixedPhase, ΔW2::MixedPhase) = wavefront_expansion(*, ΔW1, ΔW2)

function exponentiate(φ::Phase, n::Integer)
    φ::WavefrontError = φ
    ρ, θ = polar()
    φ_n = @. φ(ρ, θ) ^ n
    n_max = φ.n_max * n
    W(φ_n, n_max; φ.precision)
end

function params(ΔW::MixedPhase)
    return (maximum(W -> getfield(W, i), ΔW.W) for i = (:n_max, :precision))
end

function wavefront_expansion(f, φ1::MixedPhase, φ2::MixedPhase)
    ρ, θ = polar()
    φ3 = @. f(φ1(ρ, θ), φ2(ρ, θ))
    n1, p1 = params(φ1)
    n2, p2 = params(φ2)
    n_max = (f == *) ? n1 + n2 : max(n1, n2)
    precision = max(p1, p2)
    W(φ3, n_max; precision)
end

function show(io::IO, W::T) where {T <: MixedPhase}
    println(io, T, "\nW field:")
    show(io, "text/plain", W.W)
end

function show(io::IO, m::MIME"text/plain", W::T) where {T <: MixedPhase}
    if haskey(io, :typeinfo)
        print(summary(W.W))
        return
    else
        show(io, W)
    end
end
