using Test
using Zernike
using Zernike: Z, W, P, radicand, Φ, get_i, canonical, coords, reconstruct,
               validate_length, map_phase, format_strings, get_mn, LaTeXString,
               latexstring, J, metrics, polar

@testset "fringe" begin
    @test_throws "37" fringe_to_j(38)
    fringe = [1, 3, 2, 6, 4, 5, 11, 8, 7, 10, 18, 13, 9, 12, 17]
    @test fringe_to_j.(fringe) == 0:14
    v_fringe = collect(1.0:18.0)
    N_ = (√radicand(m, n) for n = 0:4 for m = -n:2:n)
    compare = map(N_, standardize(v_fringe), v_fringe[fringe]) do N, v, o
                  N * v ≈ o
              end
    @test compare |> all
end

@testset "noll" begin
    @test_throws "index" noll_to_j(-1)
    noll = [1, 3, 2, 5, 4, 6, 9, 7, 8, 10, 15, 13, 11, 12, 14]
    @test noll_to_j.(noll) == 0:14
    v = zeros(15)
    v[noll] = 1:15
    @test !isequal(v, 1:15)
    standardize!(v)
    @test isequal(v, 1:15)
end

function compare_coefficients(j_max)
    γ = Vector{Float64}[]
    zλ = Vector{Vector{Float64}}(undef, j_max + 1)
    Φλ = similar(zλ)
    mn = Tuple{Int, Int}[]
    kv = Int[]
    for j = 0:j_max
        m, n = get_mn(j)
        μ = abs(m)
        k = (n - μ) ÷ 2
        if m ≥ 0
            push!(mn, (m, n))
            push!(kv, k)
            push!(γ, canonical(μ, n, k))
        end
        i = j + 1
        zλ[i] = zernike(m, n, Model).R.λ
        Φλ[i] = Zernike.coefficients(μ, n)[end]
    end
    m_max, n_max = last(mn)
    λ = Φ(m_max, n_max)
    for (i, v) ∈ pairs(λ)
        n = mn[i][2]
        λ[i] = v[range(start = n + 1, step = -2, length = kv[i] + 1)]
    end
    @test get_i(m_max, n_max) == length(γ) == length(λ)
    @test zλ == Φλ
    γ == λ
end

@testset "radial coefficients" begin
    @test_throws "Bounds" Φ(-1, 1)
    @test_throws "Bounds" Φ(4, 0)
    @test_throws "Bounds" Φ(2, 5)
    @test Z(4, 4, Model).R.γ == [1.0]
    @test get_i(101, 101) == length(Φ(101, 101))
    i = 0
    for n = 0:800, m = (n % 2):2:n
        i += 1
    end
    @test [get_i(m, n) for n = 0:800 for m = (n & 1):2:n] == 1:i
    γ = [20, -30, 12, -1]
    @test canonical(0, 6, 3) == γ
    @test Φ(0, 6)[end][7:-2:1] == γ
    @test Z(0, 6, Model).R.γ == γ
    j = 230
    @test compare_coefficients(j)
end

Z62 = Z(2, 6, Model)

@testset "zernike polynomials" begin
    @test_throws "j" Z(-1)
    @test_throws "Bounds" Z(1, -1)
    @test_throws "Bounds" Z(-5, 3)
    @test_throws "Bounds" Z(0, 3)
    polynomial(ρ, θ) = √14 * (15ρ^6 - 20ρ^4 + 6ρ^2)cos(2θ)
    @test Z62(0.3, 0.7) ≈ polynomial(0.3, 0.7)
    (; R) = Z(0, 30, Model)
    @test R(1.0) ≈ sum(R.γ) ≈ sum(R.λ) ≈ 1.0
end

ρ = range(0.0, 1.0, 21)
θ = range(0, 2π, 21)
z = Z62.(ρ', θ)
OPD = z + [10sinc(5r) for i = 1:21, r ∈ ρ]
r, t = coords(OPD)
OPD_vec = vec(OPD)
v = reconstruct(r, t, OPD_vec, 8)[1]

@testset "wavefront error" begin
    @test_throws "number" validate_length(ones(5))
    @test_throws "length" W(zeros(2), zeros(2), zeros(3), 4)
    v_sub = [2.0, 4.0, 6.0, 6.2, 8.0]
    orders = [(0, 2), (0, 4), (0, 6), (2, 6), (0, 8)]
    v_full = standardize(v_sub, orders)
    @test length(v_full) == 45
    @test getindex(v_full, [5, 13, 25, 26, 41]) == v_sub
    ΔW = W(z, 6, Model; precision = 0)
    @test ΔW(0.3, 0.7) ≈ Z62(0.3, 0.7)
    @test map_phase(r, t, OPD_vec) == (ρ, θ, OPD)
    inds_a = getfield(W(OPD, 8, Model; precision = 7), :recap)
    @testset "expansion coefficients" for i in inds_a
        @test i[:a] ≈ v[i[:j] + 1] atol = 1e-7
    end
    Z44 = Z(4, 4, Model)
    ρ_rs = rand(2^8)
    θ_rs = 2π * rand(2^8)
    OPD_matrix = 7.0 * Z44.(ρ_rs', θ_rs)
    OPD_coords = [(i, j) for j in θ_rs, i in ρ_rs]
    unrolled_coords = (vec(getindex.(OPD_coords, i)) for i = 1:2)
    ΔW44_vector_phase = W(unrolled_coords..., vec(OPD_matrix), 4, Model)
    ΔW44_matrix_phase = W(ρ_rs, θ_rs, OPD_matrix, 4, Model)
    @test ΔW44_vector_phase.v == ΔW44_matrix_phase.v
    a = rand(5)
    named_orders = [(m = i[1], n = i[2]) for i in orders]
    j_orders = [4, 12, 24, 25, 40]
    W1 = WavefrontError(a)
    W2 = WavefrontError(orders, a)
    W3 = WavefrontError(named_orders, a)
    W4 = WavefrontError(j_orders, a)
    @test eltype(W1.recap) <: NamedTuple
    @test W1.a == W1.v
    @test isempty(W1.fit_to)
    @test W1.n_max == 2
    @testset "constructor equality" for i in (:recap, :v, :n_max, :fit_to, :a)
        @test getfield(W2, i) == getfield(W3, i) == getfield(W4, i)
    end
    idx_orders = j_orders .+ 1
    @test W2.v[idx_orders] == a
    @test iszero(W2.v[setdiff(1:45, idx_orders)])
    @test isempty(W2.fit_to)
    @test W2.n_max == 8
end

@testset "metrics" begin
    r = range(0.0, 1.0, 2^10)
    w = [0.25r^2 for _ ∈ r, r ∈ r]
    ΔW = W(w, 2, Model)
    (; pv, rms, strehl) = metrics(ΔW)
    @test pv ≈ 1/4 atol = 1e-2
    @test rms ≈ 1/14 atol = 1e-2
    @test strehl ≈ 0.8 atol = 1e-1
end

@testset "format strings" begin
    latex1, latex2, unicode = format_strings(Z62)
    @test typeof(latex1) == typeof(latex2) == LaTeXString
    @test unicode == "√(14)(15ρ⁶ − 20ρ⁴ + 6ρ²)cos(2θ)"
    recap = []
    for i = 1:10
        j = i - 1
        m, n = get_mn(j)
        push!(recap, (; j, n, m, a = -float(i)))
    end
    abbreviated = "ΔW ≈ -10.000Z_{3}^{3} - 9.000Z_{3}^{1} \
                        - 8.000Z_{3}^{-1} - 7.000Z_{3}^{-3}..."
    @test format_strings(recap) == latexstring(abbreviated)
end

ΔW = W(OPD, 8, Model; precision = "full")
ε = 0.75

@testset "scale" begin
    @test_throws "Bounds" J(v, -1.0, Model)
    @test_throws "Bounds" J(v, 2.0, Model)
    @test_throws "Bounds" P(v, -1.0, Model)
    @test_throws "Bounds" P(v, 2.0, Model)
    ΔW_J_s = J(v, ε, Model; precision = "full")
    ΔW_P_s = P(v, ε, Model; precision = "full")
    @test ΔW_J_s.a ≈ ΔW_P_s.a
    @test ΔW(ε, π/2) ≈ ΔW_J_s(1.0, π/2)
    @test ΔW(ε, π/2) ≈ ΔW_P_s(1.0, π/2)
end

δ = 0.15 + 0.15im
ρ_t, θ_t = polar(δ)
θ1 = 0.68π
φ = π - θ1 + θ_t
ρ2 = √(ε^2 + ρ_t^2 - 2*ε*ρ_t*cos(φ))
θ2 = asin(ε * sin(φ) / ρ2) + θ_t

@testset "translate" begin
    @test_throws "Bounds" P(v, -1.0, 0.3 + 0.0im)
    @test_throws "Bounds" P(v, 1.0, 0.1 + 0.0im)
    ΔW_t = P(v, ε, δ, Model; precision = "full")
    @test ΔW_t(1.0, θ1) ≈ ΔW(ρ2, θ2)
    @test ΔW_t(0.0, 0.0) ≈ ΔW(ρ_t, θ_t)
end

ϕ = 0.64
θ1 = 0.21
φ = π - θ1 + θ_t - ϕ
ρ2 = √(ε^2 + ρ_t^2 - 2*ε*ρ_t*cos(φ))
θ2 = asin(ε * sin(φ) / ρ2) + θ_t

@testset "rotate" begin
    ΔW_r1 = P(v, 1.0, 0.0im, ϕ, Model; precision = "full")
    @test ΔW_r1(1.0, 0.0) ≈ ΔW(1.0, ϕ)
    ΔW_r2 = P(v, ε, δ, ϕ, Model; precision = "full")
    @test ΔW_r2(1.0, θ1) ≈ ΔW(ρ2, θ2)
end

@testset "elliptical" begin
    ΔW_e = P([0.0, 0.0, 1.0], 1.0, 0.0im, 0.0, (0.148, 0.0), Model)
    @test getfield(ΔW_e, :recap) == [(j = 2, n = 1, m = 1, a = 0.148)]
end
