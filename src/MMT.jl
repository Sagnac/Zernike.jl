#= The following algorithm is based on an analysis presented in:

https://arxiv.org/abs/2411.18361

Matthieu Cadiot, Jonathan Jaquette, Jean-Philippe Lessard, Akitoshi Takayasu. Validated matrix multiplication transform for orthogonal polynomials with applications to computer-assisted proofs for PDEs. (2024)

=#

using Jacobi: jacobi, jacobi_zeros

P(x, n, μ) = jacobi(x, n, 0, μ)

nodes(n, μ) = jacobi_zeros(n, 0, μ)

function product_expansion(m_1::Int, a::Vector, m_2::Int, b::Vector)
    μ_1 = abs(m_1)
    μ_2 = abs(m_2)
    n_1 = radial_n_max(μ_1, a)
    n_2 = radial_n_max(μ_2, b)
    m_3 = m_1 + m_2
    m_bar = m_1 * m_2 < 0 ? min(μ_1, μ_2) : 0
    n_3 = n_1 + n_2 + m_bar
    len = n_3 + 1
    μ_3 = abs(m_3)
    x = nodes(len, μ_3)
    MMT_1 = [P(xᵢ, n, μ_1) for xᵢ ∈ x, n ∈ 0:n_1]
    MMT_2 = [P(xᵢ, n, μ_2) for xᵢ ∈ x, n ∈ 0:n_2]
    MMT_3 = [P(xᵢ, n, μ_3) for xᵢ ∈ x, n ∈ 0:n_3]
    iMMT = inv(MMT_3)
    dn_1 = n_3 - n_1
    dn_2 = n_3 - n_2
    dl_1 = len - length(a)
    dl_2 = len - length(b)
    MMT_1 = hcat(MMT_1, zeros(len, dn_1))
    MMT_2 = hcat(MMT_2, zeros(len, dn_2))
    a_padded = vcat(a, zeros(dl_1))
    b_padded = vcat(b, zeros(dl_2))
    f_0 = @. ((x + 1) / 2) ^ m_bar
    f_1 = MMT_1 * a_padded
    f_2 = MMT_2 * b_padded
    c = iMMT * (f_0 .* f_1 .* f_2)
    return m_3, c
end

const threshold = nextfloat(0.0)

function *(w1::Wavefront{T}, w2::Wavefront{T};
              threshold = threshold) where T <: RadialPolynomial
    m_1 = w1.m
    m_2 = w2.m
    a = w1.a
    b = w2.v
    m_3, c = product_expansion(m_1, a, m_2, b)
    c = sieve(c, threshold)
    Wavefront{T}(m_3, c)
end
