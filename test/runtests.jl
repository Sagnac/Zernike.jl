using Test
using Zernike
using Zernike: radicand, Φ, get_i, λ, coords, Wf, validate_length, map_phase,
               format_strings, get_mn, LaTeXString, latexstring, J

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
    @test [λ(0, 6, s, 3) for s = 0:3] == γ
    @test Φ(0, 6)[end][7:-2:1] == γ
    @test Z(0, 6, Model).R.γ == γ
end

Z62 = Z(2, 6, Model)

@testset "zernike polynomials" begin
    @test_throws "j" Z(-1)
    @test_throws "Bounds" Z(1, -1)
    @test_throws "Bounds" Z(-5, 3)
    @test_throws "Bounds" Z(0, 3)
    polynomial(ρ, θ) = √14 * (15ρ^6 - 20ρ^4 + 6ρ^2)cos(2θ)
    @test Z62(0.3, 0.7) ≈ polynomial(0.3, 0.7)
end

ρ = range(0.0, 1.0, 21)
θ = range(0, 2π, 21)
Zp = Z62.(ρ', θ)
OPD = Zp + [10sinc(5r) for i = 1:21, r ∈ ρ]
r, t = coords(OPD)
OPD_vec = vec(OPD)
v = Wf(r, t, OPD_vec, 8)[1]

@testset "wavefront error" begin
    @test_throws "number" validate_length(ones(5))
    @test_throws "length" W(zeros(2), zeros(2), zeros(3), 4)
    v_sub = [2.0, 4.0, 6.0, 6.2, 8.0]
    orders = [(0, 2), (0, 4), (0, 6), (2, 6), (0, 8)]
    v_full = standardize(v_sub, orders)
    @test length(v_full) == 45
    @test getindex(v_full, [5, 13, 25, 26, 41]) == v_sub
    ΔW = W(Zp, 6, Model; precision = 0)
    @test ΔW(0.3, 0.7) ≈ Z62(0.3, 0.7)
    @test map_phase(r, t, OPD_vec) == (ρ, θ, OPD)
    inds_a = getfield(W(OPD, 8, Model; precision = 7), :i)
    @testset "expansion coefficients" for i in inds_a
        @test i[:a] ≈ v[i[:j] + 1] atol = 1e-7
    end
end

@testset "format strings" begin
    latex1, latex2, unicode = format_strings(Z62)
    @test typeof(latex1) == typeof(latex2) == LaTeXString
    @test unicode == "√(14)(15ρ⁶ − 20ρ⁴ + 6ρ²)cos(2θ)"
    a = []
    for i = 1:10
        j = i - 1
        m, n = get_mn(j)
        push!(a, (; j, n, m, a = -float(i)))
    end
    abbreviated = "ΔW ≈ -10.000Z_{3}^{3} - 9.000Z_{3}^{1} \
                        - 8.000Z_{3}^{-1} - 7.000Z_{3}^{-3}..."
    @test format_strings(a) == latexstring(abbreviated)
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
ρ_t, θ_t = abs(δ), angle(δ)
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
    @test getfield(ΔW_e, :i) == [(j = 2, n = 1, m = 1, a = 0.148)]
end
