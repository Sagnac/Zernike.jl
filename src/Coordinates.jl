function polar((u, v)::NTuple{2, Int})
    θ = range(0.0, 2π, u)
    ρ = range(0.0, 1.0, v)
    return ρ, θ
end

function polar(x::FloatVec, y::FloatVec)
    ρ = hypot.(x, y)
    θ = atan.(y, x)
    return ρ, θ
end

function polar(d::Int = d_fit)
    ρ, θ = polar((d, d))
    return ρ', θ
end

function polar(m::Int, n::Int; finesse::Int = finesse)
    m = abs(m)
    ρ = range(0.0, 1.0, finesse)
    θ = range(0.0, 2π, finesse)
    return ρ, θ
end

polar(z::Complex) = abs(z), angle(z)

function polar_mat(ρ, θ)
    x = @. ρ' * cos(θ)
    y = @. ρ' * sin(θ)
    return x, y
end

# extract pupil coordinates
function coords(ρ::FloatVec, θ::FloatVec)
    ρ2 = ρ ⊗ ones(length(θ))
    θ2 = ones(length(ρ)) ⊗ θ
    return ρ2, θ2
end

function cartesian_coords(x::FloatVec, y::FloatVec)
    y, x = coords(y, x)
    return x, y
end

coords(OPD::FloatMat) = coords(polar(size(OPD))...)
