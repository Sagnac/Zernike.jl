#= This algorithm is based on a mathematical formulation presented in:

https://doi.org/10.1364/JOSAA.24.000569

Linda Lundström and Peter Unsbo, "Transformation of Zernike coefficients: scaled, translated, and rotated wavefronts with circular and elliptical pupils," J. Opt. Soc. Am. A 24, 569-577 (2007)

https://opg.optica.org/josaa/abstract.cfm?uri=josaa-24-3-569

=#

const b = binomial
const ℯ = cis # ℯ(x) = exp(im*x)

function transform(v::Vector{T}, ε::T, δ::Complex{T}, ϕ::T) where T <: Float64
    !(0.0 ≤ ε + abs(δ) ≤ 1.0) && error("Bounds: 0.0 ≤ ε + abs(δ) ≤ 1.0\n")
    len = length(v)
    n_max = get_n(len - 1)
    remap = Dict{Tuple{Int, Int}, Int}()
    order = Tuple{Int, Int, Int}[]
    # normalization factors for complex Zernike polynomials
    N = zeros(Float64, len, len)
    # radial coefficient block-diagonal matrix
    R = copy(N)
    if iszero(δ)
        # scaling
        ηₛ = copy(N)
    end
    if !iszero(ϕ)
        # rotation
        ηᵣ = zeros(ComplexF64, len, len)
    end
    i = 0
    for m = -n_max:n_max
        μ = abs(m)
        for n = μ:2:n_max
            i += 1
            k1 = (n - μ) ÷ 2
            push!(remap, (m, n) => i)
            push!(order, (get_j(m, n) + 1, get_j(-m, n) + 1, sign(m)))
            N[i,i] = √(n+1)
            λ = Φ(μ, n)
            for s = 0:k1
                setindex!(R, λ[n-2s+1], remap[(m, n-2s)], i)
            end
            if !iszero(ϕ)
                ηᵣ[i,i] = ℯ(m * ϕ)
            end
            if iszero(δ)
                ηₛ[i,i] = ε ^ n
            end
        end
    end
    if !iszero(δ)
        ηₜ = translate(ε, δ, remap)
        if !iszero(ϕ)
            # net transformation matrix
            η = ηᵣ * ηₜ
        else
            η = ηₜ
        end
    elseif !iszero(ϕ)
        η = ηᵣ * ηₛ
    else
        η = ηₛ
    end
    c = to_complex(v, order)
    # conversion matrix
    C = (R * N) \ (η * R * N)
    # transformed set of complex Zernike expansion coefficients
    c′ = C * c
    # transformed set of standard Zernike expansion coefficients
    v2 = to_real(c′, order)
    return v2, n_max
end

function translate(ε::Float64, δ::ComplexF64, remap::Dict)
    ρₜ, θₜ = abs(δ), angle(δ)
    len = length(remap)
    n_max = get_n(len - 1)
    ηₜ = zeros(ComplexF64, len, len)
    i = 0
    for m = -n_max:n_max, n = abs(m):2:n_max
        i += 1
        k2 = (n + m) ÷ 2
        k3 = (n - m) ÷ 2
        for p = 0:k2, q = 0:k3
            n′ = n - p - q
            m′ = m - p + q
            z = b(k2,p) * b(k3,q) * ε^n′ * ρₜ^(p+q) * ℯ((p-q)θₜ)
            setindex!(ηₜ, z, remap[(m′, n′)], i)
        end
    end
    return ηₜ
end

function to_complex(v::Vector{Float64}, order::Vector{Tuple{Int, Int, Int}})
    c = Vector{Complex{Float64}}(undef, length(v)) 
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
    v2 = Vector{Float64}(undef, length(c2))
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

function S(v::Vector{T}, ε::T, δ::Complex{T} = 0.0im, ϕ::T = 0.0;
           precision = 3, scale::Int = 101) where T <: Float64
    v2, n_max = transform(v, ε, δ, ϕ)
    Zᵢ = Vector{Polynomial}(undef, length(v))
    ΔW, b = Ψ(v2, Zᵢ, n_max; precision)
    Λ(ΔW, b, v2, n_max; scale)
end

function S(v::Vector{T}, ε::T, ::Model;
           precision = 3) where T <: Float64
    v2, n_max = transform(v, ε, 0.0im, zero(T))
    Zᵢ = Vector{Polynomial}(undef, length(v))
    Ψ(v2, Zᵢ, n_max; precision)[1]
end

function S(v::Vector{T}, ε::T, δ::Complex{T}, ::Model;
           precision = 3) where T <: Float64
    v2, n_max = transform(v, ε, δ, zero(T))
    Zᵢ = Vector{Polynomial}(undef, length(v))
    Ψ(v2, Zᵢ, n_max; precision)[1]
end

function S(v::Vector{T}, ε::T, δ::Complex{T}, ϕ::T, ::Model;
           precision = 3) where T <: Float64
    v2, n_max = transform(v, ε, δ, ϕ)
    Zᵢ = Vector{Polynomial}(undef, length(v))
    Ψ(v2, Zᵢ, n_max; precision)[1]
end
