#= This algorithm is based on a recursive relation presented in:

https://doi.org/10.1364/OL.38.002487

B. Honarvar Shakibaei and R. Paramesran, "Recursive formula to compute Zernike radial polynomials," Opt. Lett.  38, 2487-2489 (2013).

https://opg.optica.org/ol/abstract.cfm?uri=ol-38-14-2487

=#

import ShiftedArrays: circshift as shift

function get_i(m_max, n_max)
    n_mod_2 = isodd(n_max)
    i_max = ((n_max + 2)n_max + n_mod_2) ÷ 4 + (m_max + !n_mod_2 + 1) ÷ 2
    return i_max
end

function Φ(m_max::Int, n_max::Int)
    if m_max < 0 || n_max < 0 || m_max > n_max || isodd(n_max - m_max)
        error(
            """
            In method: Φ(m_max, n_max)
            Bounds:
            m_max ≥ 0
            n_max ≥ 0
            m_max ≤ n_max
            n_max - m_max even
            """
        )
    end
    if m_max == n_max
        λᵢ = zeros(Float64, n_max + 1)
        λᵢ[end] = 1.0
        return λᵢ
    end
    λ = Vector{Float64}[]
    i = 0
    n_even = true
    for n = 0:n_max
        for m = !n_even:2:ifelse(n ≠ n_max, n, m_max)
            i += 1
            if m == n
                λᵢ = zeros(Float64, n_max + 1)
                λᵢ[n+1] = 1.0
                n_even = !n_even
            elseif m == 0
                λᵢ = 2shift(λ[i-n÷2], 1) - λ[i-n]
            else
                Δ = (n + 1 + n_even) ÷ 2
                λᵢ = shift(λ[i-Δ] + λ[i-Δ+1], 1) - λ[i-n]
            end
            push!(λ, λᵢ)
        end
    end
    return λ[end]
end
