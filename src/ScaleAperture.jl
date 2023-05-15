#= This algorithm is based on a formula presented by Janssen and Dirksen in:

https://doi.org/10.2971/jeos.2007.07012

Janssen, A., & Dirksen, P. (2007). Computing Zernike polynomials of arbitrary degree using the discrete Fourier transform. Journal Of The European Optical Society - Rapid Publications, 2.

https://www.jeos.org/index.php/jeos_rp/article/view/07012

=#

function S(ε, v)

    j_max::Int = length(v) - 1
    n_max::Int = get_n(j_max)

    b = Vector{Float64}(undef, j_max + 1)

    n = 0
    m = 0

    @inbounds for i in eachindex(v)

        V = Float64[]

        a = v[i]
        Nmn = √radicand(m, n)
        R0 = Z(n, n, Fit()).R(ε)
        push!(V, a * Nmn * R0)

        ii = i

        for n′ = n+2:2:n_max

            ii += 2n′

            a = v[ii]
            N = √radicand(m, n′)
            R1 = Z(n, n′, Fit()).R(ε)
            R2 = Z(n+2, n′, Fit()).R(ε)

            push!(V, a * N * (R1 - R2))

        end

        b[i] = ∑(V) / Nmn

        m += 2

        if m > n
            n += 1
            m = -n
        end

    end

    return b

end
