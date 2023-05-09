#= This algorithm is based on a recursive relation presented in:

https://doi.org/10.1364/OL.38.002487

B. Honarvar Shakibaei and R. Paramesran, "Recursive formula to compute Zernike radial polynomials," Opt. Lett.  38, 2487-2489 (2013).

https://opg.optica.org/ol/abstract.cfm?uri=ol-38-14-2487

=#

import ShiftedArrays: circshift as shift

function Φ(n_max::Integer, m_max::Integer)

    n_mod_2 = isodd(n_max)

    i_max = ((n_max + 2)n_max + n_mod_2) ÷ 4 + (m_max + !n_mod_2 + 1) ÷ 2

    λ = [zeros(Float64, n_max + 1) for i = 1:i_max]

    i = 0
    n_even = true

    @inbounds for n = 0:n_max

        for m = !n_even:2:ifelse(n ≠ n_max, n, m_max)

            i += 1

            if m == n
                λ[i][n+1] = 1.0
                n_even = !n_even
            elseif m == 0
                λ[i] = 2shift(λ[i-n÷2], 1) - λ[i-n]
            else
                Δ = (n + 1 + n_even) ÷ 2
                λ[i] = shift(λ[i-Δ] + λ[i-Δ+1], 1) - λ[i-n]
            end

        end

    end

    return λ[end]

end
