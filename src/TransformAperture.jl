function transform(v::Vector{Float64})
    n_max = get_n(length(v) - 1)
    order = Tuple{Int, Int, Int}[]
    for m = -n_max:n_max
        for n = abs(m):2:n_max
            j = get_j(m, n)
            push!(order, (j + 1, get_j(-m, n) + 1, sign(m)))
        end
    end
    c = to_complex(v, order)
    v2 = to_real(c, order)
    return c, v2
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
