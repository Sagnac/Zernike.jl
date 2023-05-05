# Estimates wavefront error by expressing the aberrations as a linear combination
# of weighted Zernike polynomials. The representation is approximate,
# especially for a small set of data.
function W(r::Vector, ϕ::Vector, OPD::Vector, n_max::Integer)

    if !allequal(length.((r, ϕ, OPD)))
        @error "Vectors must be of equal length."
        return
    end

    j_max = get_j(n_max, n_max)

    Zᵢ = [Z(j; mode = "fit") for j = 0:j_max]

    # linear least squares
    A = reduce(hcat, Zᵢ[i][:Z].(r, ϕ) for i = 1:j_max+1)

    # Zernike expansion coefficients
    v = A \ OPD

    a = NamedTuple[]

    # store the non-trivial coefficients
    for (i, aᵢ) in pairs(v)
        aᵢ = round(aᵢ; digits = 3)
        if !iszero(aᵢ)
            push!(a, (j = i - 1, n = Zᵢ[i].n, m = Zᵢ[i].m, a = aᵢ))
        end
    end

    # create the fitted polynomial
    ΔW(ρ, θ) = ∑(v[i] * Zᵢ[i][:Z](ρ, θ) for i = 1:j_max+1)

    # construct the estimated wavefront error
    ΔWp = ΔW.(ρ', θ)

    Z_LaTeX = "ΔW ≈ "

    function ζ(i)
        aᵢ = a[i][:a]
        (i > 1 ? abs(aᵢ) : aᵢ), "Z_{$(a[i][:n])}^{$(a[i][:m])}"
    end

    function η(index, a)
        for i = 1:length(a)-1
            Z_LaTeX *= string(ζ(index[i])..., a[i+1][:a] > 0 ? " + " : " - ")
        end
    end

    if length(a) > 8
        sorted_indices = sortperm(abs.([aᵢ[:a] for aᵢ ∈ a]); rev = true)
        subset_a = getindex(a, sorted_indices[1:4])
        η(sorted_indices, subset_a)
        Z_LaTeX *= string(ζ(sorted_indices[4])...) * "..."
    else
        η(keys(a), a)
        Z_LaTeX *= string(ζ(lastindex(a))...)
    end

    titles = (plot = Z_LaTeX, window = "Estimated wavefront error")

    # Peak-to-valley wavefront error
    PV = maximum(ΔWp) - minimum(ΔWp)

    # RMS wavefront error
    # where σ² is the variance (second central moment about the mean)
    # the mean is the first a00 piston term
    σ = ∑(v[i]^2 for i = 2:lastindex(v)) |> sqrt

    Strehl_ratio = exp(-(2π * σ)^2)

    metrics = (PV = PV, RMS = σ, Strehl = Strehl_ratio)

    fig = ZPlot(ΔWp; titles...)

    return fig, a, v, metrics

end

# methods

W(; r, t, OPD, n_max) = W(r, t, OPD, n_max)

function W(x::Vector, y::Vector, OPD::Vector; n_max::Integer)
    r = hypot.(x, y)
    ϕ = atan.(y, x)
    W(r, ϕ, OPD, n_max)
end

function W(data::Matrix, n_max::Integer)
    r, ϕ, OPD = [data[:, i] for i = 1:3]
    W(r, ϕ, OPD, n_max)
end

W(data; n_max) = W(data, n_max)
