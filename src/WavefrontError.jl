# Estimates wavefront error by expressing the aberrations as a linear combination
# of weighted Zernike polynomials. The representation is approximate,
# especially for a small set of data.
function Wf(r::Vector, ϕ::Vector, OPD::Vector, n_max::Integer; precision = 3)

    if !allequal(length.((r, ϕ, OPD)))
        error("Vectors must be of equal length.\n")
    end

    j_max = get_j(n_max, n_max)

    Zᵢ = Function[]

    inds = @NamedTuple{j::Int, n::Int, m::Int}[]

    for j = 0:j_max
        Zⱼ, j_n_m = Z(j, Fit())
        push!(Zᵢ, Zⱼ)
        push!(inds, j_n_m)
    end

    # linear least squares
    A = reduce(hcat, Zᵢ[i].(r, ϕ) for i = 1:j_max+1)

    # Zernike expansion coefficients
    v::Vector{Float64} = A \ OPD

    a = @NamedTuple{j::Int, n::Int, m::Int, a::Float64}[]

    # store the non-trivial coefficients
    for (i, aᵢ) in pairs(v)
        aᵢ = ifelse(precision == "full", aᵢ, round(aᵢ; digits = precision))
        if !iszero(aᵢ)
            push!(a, (; inds[i]..., a = aᵢ))
        end
    end

    # create the fitted polynomial
    ΔW(ρ, θ)::Float64 = ∑(aᵢ[:a] * Zᵢ[aᵢ.j+1](ρ, θ) for aᵢ ∈ a)

    return a, v, ΔW

end

function metrics(ΔWp, v)

    # Peak-to-valley wavefront error
    PV = maximum(ΔWp) - minimum(ΔWp)

    # RMS wavefront error
    # where σ² is the variance (second central moment about the mean)
    # the mean is the first a00 piston term
    σ = ∑(v[i]^2 for i = 2:lastindex(v)) |> sqrt

    Strehl_ratio = exp(-(2π * σ)^2)

    return (PV = PV, RMS = σ, Strehl = Strehl_ratio)

end

function W(r::T, ϕ::T, OPD::T, n_max::Integer; options...) where T <: Vector

    if haskey(options, :precision) &&
       (options[:precision] == "full" || options[:precision] isa Integer)
        precision = options[:precision]
    else
        precision = 3
    end

    a, v, ΔW = Wf(r, ϕ, OPD, n_max; precision)

    if haskey(options, :scale) &&
       (options[:scale] isa Integer && options[:scale] ∈ 1:100)
       scale = options[:scale]
    else
       scale = ceil(Int, 100 / √ length(a))
    end

    ρ, θ = polar(n_max, n_max; scale)

    # construct the estimated wavefront error
    ΔWp = ΔW.(ρ', θ)

    W_LaTeX = format_strings(a)

    titles = (plot = W_LaTeX, window = "Estimated wavefront error")

    fig = ZPlot(ρ, θ, ΔWp; titles...)
        
    return a, v, metrics(ΔWp, v), fig

end

# methods

W(; r, t, OPD, n_max, options...) = W(r, t, OPD, n_max; options...)

function W(x::Vector, y::Vector, OPD::Vector; n_max::Integer, options...)
    r = hypot.(x, y)
    ϕ = atan.(y, x)
    W(r, ϕ, OPD, n_max; options...)
end

function W(data::Matrix, n_max::Integer; options...)
    r, ϕ, OPD = [data[:, i] for i = 1:3]
    W(r, ϕ, OPD, n_max; options...)
end

W(data; n_max, options...) = W(data, n_max; options...)
