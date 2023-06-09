function transform(v::Vector{Float64})
    j_max = length(v) - 1
    n_max = get_n(j_max)
    order = Tuple{Int, Int, Int}[]
    remap = Dict{Tuple{Int, Int}, Int}()
    N = zeros(Float64, j_max + 1, j_max + 1)
    R = zeros(Float64, j_max + 1, j_max + 1)
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
            i += 1
        end
    end
    return N, R
end

function to_complex(v::Vector{Float64}, order::Vector{Tuple{Int, Int, Int}})
    c = Vector{Complex{Float64}}(undef, length(v)) 
    for (i, j) in pairs(order)
        if j[3] < 0
            c[i] = (v[j[2]] + im * v[j[1]]) / √2
        elseif j[3] > 0
            c[i] = (v[j[1]] - im * v[j[2]]) / √2
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
