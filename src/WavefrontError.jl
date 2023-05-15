# Estimates wavefront error by expressing the aberrations as a linear combination
# of weighted Zernike polynomials. The representation is approximate,
# especially for a small set of data.

struct WavefrontError
    a::Vector{Float64}
    Z::Vector{Polynomial}
end

function (ΔW::WavefrontError)(ρ, θ)::Float64
    (; a, Z) = ΔW
    ∑(aᵢ * Z[i](ρ, θ) for (i, aᵢ) ∈ pairs(a))
end

# construction function
function Wf(r::Vector, ϕ::Vector, OPD::Vector, n_max::Int; precision = 3)

    if !allequal(length.((r, ϕ, OPD)))
        error("Vectors must be of equal length.\n")
    end

    j_max = get_j(n_max, n_max)

    Zᵢ = Polynomial[Z(j, Fit()) for j = 0:j_max]

    # linear least squares
    A = reduce(hcat, Zᵢ[i].(r, ϕ) for i = 1:j_max+1)

    # Zernike expansion coefficients
    v::Vector{Float64} = A \ OPD

    ΔW, a = Ξ(v, Zᵢ; precision)

    return ΔW, a, v

end

# filtering function
function Ξ(v, Zᵢ; precision)

    a = @NamedTuple{j::Int, n::Int, m::Int, a::Float64}[]

    av = Float64[]

    Zₐ = Polynomial[]

    clipped = !isassigned(Zᵢ, 1)

    # store the non-trivial coefficients
    for (i, aᵢ) in pairs(v)
        aᵢ = ifelse(precision == "full", aᵢ, round(aᵢ; digits = precision))
        if !iszero(aᵢ)
            if clipped
                Zᵢ[i] = Z(i-1, Fit())
            end
            push!(a, (; Zᵢ[i].inds..., a = aᵢ))
            push!(av, aᵢ)
            push!(Zₐ, Zᵢ[i])
        end
    end

    # create the fitted polynomial
    ΔW = WavefrontError(av, Zₐ)

    return ΔW, a

end

# synthesis function
function Λ(ΔW, a, v, n_max; scale)

    ρ, θ = polar(n_max, n_max; scale)

    # construct the estimated wavefront error
    ΔWp = ΔW.(ρ', θ)

    W_LaTeX = format_strings(a)

    titles = (plot = W_LaTeX, window = "Estimated wavefront error")

    fig = ZPlot(ρ, θ, ΔWp; titles...)

    return a, v, metrics(v, ΔWp), fig

end

function metrics(v, ΔWp)

    # Peak-to-valley wavefront error
    PV = maximum(ΔWp) - minimum(ΔWp)

    # RMS wavefront error
    # where σ² is the variance (second central moment about the mean)
    # the mean is the first a00 piston term
    σ = ∑(v[i]^2 for i = 2:lastindex(v)) |> sqrt

    Strehl_ratio = exp(-(2π * σ)^2)

    return (PV = PV, RMS = σ, Strehl = Strehl_ratio)

end

# main interface function
function W(r::T, ϕ::T, OPD::T, n_max::Int; options...) where T <: Vector

    if haskey(options, :precision) &&
       (options[:precision] == "full" || options[:precision] isa Int)
        precision = options[:precision]
    else
        precision = 3
    end

    ΔW, a, v = Wf(r, ϕ, OPD, n_max; precision)

    if haskey(options, :scale) &&
       (options[:scale] isa Int && options[:scale] ∈ 1:100)
       scale = options[:scale]
    else
       scale = ceil(Int, 100 / √ length(a))
    end

    Λ(ΔW, a, v, n_max; scale)

end

# methods

function W(r::Vector, ϕ::Vector, OPD::Vector, n_max::Int, ::Fit; options...)
    return Wf(r, ϕ, OPD, n_max; options...)[1]
end

W(; r, t, OPD, n_max, options...) = W(r, t, OPD, n_max; options...)

function W(x::Vector, y::Vector, OPD::Vector; n_max::Int, options...)
    r = hypot.(x, y)
    ϕ = atan.(y, x)
    W(r, ϕ, OPD, n_max; options...)
end

function W(data::Matrix, n_max::Int; options...)
    r, ϕ, OPD = [data[:, i] for i = 1:3]
    W(r, ϕ, OPD, n_max; options...)
end

W(data; n_max, options...) = W(data, n_max; options...)
