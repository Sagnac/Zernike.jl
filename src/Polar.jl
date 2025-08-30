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

function polar(ϵ::Int = ϵ_fit)
    ρ, θ = polar((ϵ, ϵ))
    return ρ', θ
end

function polar(m::Int, n::Int; finesse = finesse)
    m = abs(m)
    finesse = clamp(finesse, 1, 100)
    ϵ_n = finesse * (ceil(Int, π * n) + 1)
    ϵ_n = min(ϵ_n, ϵ_max)
    ϵ_m = finesse * (2m + 1)
    ϵ_m = min(ϵ_m, ϵ_max)
    ρ = range(0.0, 1.0, ϵ_n)
    θ = range(0.0, 2π, ϵ_m)
    return ρ, θ
end

polar(z::Complex) = abs(z), angle(z)

function polar_mat(ρ, θ)
    x = [ρⱼ * cos(θᵢ) for θᵢ ∈ θ, ρⱼ ∈ ρ]
    y = [ρⱼ * sin(θᵢ) for θᵢ ∈ θ, ρⱼ ∈ ρ]
    return x, y
end
