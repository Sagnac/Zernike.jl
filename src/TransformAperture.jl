const b = binomial
const ℯ = cis # ℯ(x) = exp(im*x)

function transform(ε::T, δ::Complex{T}, ϕ::T, v::Vector{T}) where T <: Float64
    len = length(v)
    n_max = get_n(len - 1)
    τ = abs(δ)
    σ = angle(δ)
    order = Tuple{Int, Int, Int}[]
    mapping = remap(n_max)
    N = zeros(Float64, len, len)
    R = copy(N)
    ηₛ = copy(N)
    ηᵣ = zeros(ComplexF64, len, len)
    ηₜ = zeros(Complex{Float64}, len, len)
    i = 1
    for m = -n_max:n_max
        μ = abs(m)
        for n = μ:2:n_max
            j = get_j(m, n)
            k1 = (n - μ) ÷ 2
            push!(order, (j + 1, get_j(-m, n) + 1, sign(m)))
            N[i,i] = √(n+1)
            λ = Φ(μ, n)
            for s = 0:k1
                setindex!(R, λ[n-2s+1], mapping[(m, n-2s)], i)
            end
            ηₛ[i,i] = ε ^ n
            ηᵣ[i,i] = ℯ(m * ϕ)
            k2 = (n + m) ÷ 2
            k3 = (n - m) ÷ 2
            for p = 0:k2, q = 0:k3
                n′ = n - p - q
                m′ = m - p + q
                z = b(k2,p) * b(k3,q) * ε^n′ * τ^(p+q) * ℯ((p-q)σ)
                setindex!(ηₜ, z, mapping[(m′, n′)], i)
            end
            i += 1
        end
    end
    c = to_complex(v, order)
    C = (R * N) \ (ηₜ * R * N)
    c′ = C * c
    v2 = to_real(c′, order)
    return v2
end

function remap(n_max)
    mapping = Dict{Tuple{Int, Int}, Int}()
    i = 1
    for m = -n_max:n_max, n = abs(m):2:n_max
        push!(mapping, (m, n) => i)
        i += 1
    end
    return mapping
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
