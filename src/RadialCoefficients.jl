#= This algorithm is based on a recursive relation presented in:

https://doi.org/10.1364/OL.38.002487

B. Honarvar Shakibaei and R. Paramesran, "Recursive formula to compute Zernike radial polynomials," Opt. Lett.  38, 2487-2489 (2013).

https://opg.optica.org/ol/abstract.cfm?uri=ol-38-14-2487

=#

import ShiftedArrays: circshift as shift

get_i(m, n) = ((n + 2)n + 1) ÷ 4 + (m + iseven(n) + 1) ÷ 2

# TODO: write a non-allocating version
function Φ(m_max::Int, n_max::Int, T::Type{<:Number} = Float64)
    μ = abs(m_max)
    @domain_check(m_max, n_max)
    m_max = μ
    i = 0
    λ = Vector{Vector{T}}(undef, get_i(m_max, n_max))
    for n = 0:n_max
        for m = isodd(n):2:ifelse(n ≠ n_max, n, m_max)
            i += 1
            if m == n
                λ[i] = zeros(T, n_max + 1)
                λ[i][n+1] = one(T)
            elseif m == 0
                λ[i] = 2shift(λ[i-n÷2], 1) - λ[i-n]
            else
                Δ = n÷2+1
                λ[i] = shift(λ[i-Δ] + λ[i-Δ+1], 1) - λ[i-n]
            end
            # @assert i == get_i(m, n)
        end
    end
    return λ
end
