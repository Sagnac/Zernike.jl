function transform(v)
    j_max::Int = length(v) - 1
    n_max::Int = get_n(j_max)
    order = []
    for m = -n_max:n_max
        μ = abs(m)
        for n = μ:2:n_max
            j = get_j(m, n)
            if m != 0
                push!(order, (j + 1, sign(m), get_j(-m, n) + 1))
            else
                push!(order, (j + 1, 0))
            end
        end
    end
    c = Vector{Complex{Float64}}(undef, j_max + 1) 
    for (i, j) in pairs(order)
        if j[2] < 0
            c[i] = (v[j[3]] + im * v[j[1]]) / √2
        elseif j[2] > 0
            c[i] = (v[j[1]] - im * v[j[3]]) / √2
        else
            c[i] = v[j[1]]
        end
    end
    v2 = Vector{Float64}(undef, j_max + 1)
    c2 = Vector{Complex{Float64}}(undef, j_max + 1) 
    for (i, j) in pairs(order)
        c2[j[1]] = c[i]
    end
    permute!(c, sortperm(order; by = first))
    for (i, j) in pairs(order)
        if j[2] < 0
            v2[j[1]] = (c[j[1]] - c[j[3]])  / √2 |> imag
        elseif j[2] > 0
            v2[j[1]] = (c[j[1]] + c[j[3]]) / √2 |> real
        else
            v2[j[1]] = c[j[1]] |> real
        end
    end
    return c, order, v2
end
