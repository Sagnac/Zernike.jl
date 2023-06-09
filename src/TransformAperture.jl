function transform(ε::Float64, v::Vector{Float64})
    len = length(v)
    n_max = get_n(len - 1)
    order = Tuple{Int, Int, Int}[]
    remap = Dict{Tuple{Int, Int}, Int}()
    N = zeros(Float64, len, len)
    R = copy(N)
    η = copy(N)
    i = 1
    for m = -n_max:n_max
        μ = abs(m)
        for n = μ:2:n_max
            j = get_j(m, n)
            k = (n - μ) ÷ 2
            push!(remap, (m, n) => i)
            push!(order, (j + 1, get_j(-m, n) + 1, sign(m)))
            N[i,i] = √(n+1)
            λ = Φ(μ, n)
            for s = 0:k
                setindex!(R, λ[n-2s+1], remap[(m, n-2s)], i)
            end
            η[i,i] = ε ^ n
            i += 1
        end
    end
    c = to_complex(v, order)
    C = (R * N) \ (η * R * N)
    c′ = C * c
    v2 = to_real(c′, order)
    return v2
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
