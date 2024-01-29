#= This algorithm is based on a formula presented in:

https://doi.org/10.2971/jeos.2007.07012

Janssen, A., & Dirksen, P. (2007). Computing Zernike polynomials of arbitrary degree using the discrete Fourier transform. Journal Of The European Optical Society - Rapid Publications, 2.

https://www.jeos.org/index.php/jeos_rp/article/view/07012

=#

function Π(v::Vector{T}, ε::T) where T <: Float64
    !(0.0 ≤ ε ≤ 1.0) && bounds("Bounds: 0.0 ≤ ε ≤ 1.0\n", ε)
    _, n_max = validate_length(v)
    v2 = similar(v, Float64)
    n = 0
    m = 0
    for i in eachindex(v)
        a = v[i]
        N_mn = √radicand(m, n)
        R0 = Z(n, n).R(ε)
        v2[i] = a * N_mn * R0
        ii = i
        for n′ = n+2:2:n_max
            ii += 2n′
            a = v[ii]
            N = √radicand(m, n′)
            R1 = Z(n, n′).R(ε)
            R2 = Z(n+2, n′).R(ε)
            v2[i] += a * N * (R1 - R2)
        end
        v2[i] /= N_mn
        if m == n
            n += 1
            m = -n
        else
            m += 2
        end
    end
    return v2, n_max
end

function scale(v::Vector{Float64}, ε::Float64;
           precision::Int = precision, finesse::Int = wavefront_finesse)
    v2, n_max = Π(v, ε)
    Zᵢ = similar(v, Polynomial)
    ΔW = Ψ(v2, Zᵢ, n_max; precision)
    Λ(ΔW; finesse)
end

function J(v::Vector{Float64}, ε::Float64; precision::Int = precision)
    v2, n_max = Π(v, ε)
    Zᵢ = similar(v, Polynomial)
    Ψ(v2, Zᵢ, n_max; precision)
end
