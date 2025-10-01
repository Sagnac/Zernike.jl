# normalization constant for Fringe coefficients
N(j::AbstractVector{Int}) = [√N²(m, n) for (m, n) ∈ get_mn.(j)]

struct Standard
    v::Vector{Float64}
end

struct Noll
    v::Vector{Float64}
end

struct Fringe
    v::Vector{Float64}
end

Standard(noll::Noll) = Standard(standardize(noll))
Standard(fringe::Fringe) = Standard(standardize(fringe))

Noll(fringe::Fringe) = (Noll ∘ Standard ∘ standardize)(fringe)

Fringe(noll::Noll) = (Fringe ∘ Standard ∘ standardize)(noll)

function Noll(s::Standard)
    sv = s.v
    j = eachindex(sv) .- 1
    noll = j_to_noll.(j)
    v = zeros(maximum(noll))
    v[noll] = sv
    return Noll(v)
end

function Fringe(s::Standard)
    sv = s.v
    j = [i for i ∈ 0:length(sv)-1 if i ∈ valid_fringes]
    fringe = j_to_fringe.(j)
    v = zeros(maximum(fringe))
    v[fringe] .= N(j) .* sv[j.+1]
    return Fringe(v)
end

function noll_to_j(noll::Int)
    @domain(noll ≥ 1, noll)
    n::Int = trunc(√(2noll - 1) + 0.5) - 1
    n_mod_2 = isodd(n)
    m::Int = 2((2noll - (n + 1)n + 1 + n_mod_2) ÷ 4) - n_mod_2
    m = flipsign(m, iseven(noll) ? 1 : -1)
    get_j(m, n)
end

function j_to_noll(j::Int)
    @domain_check_j
    m, n = get_mn(j)
    i = n % 4 ∈ 0:1
    k = m > 0 && i || m < 0 && !i ? 0 : 1
    (n + 1)n ÷ 2 + abs(m) + k
end

function fringe_to_j(fringe::Int)
    @domain(fringe ∈ 1:37, fringe)
    if fringe == 37
        return get_j(0, 12)
    end
    d = trunc(Int, √(fringe - 1)) + 1
    d2 = d^2 - fringe
    m::Int = (d2 + 1) ÷ 2
    m = flipsign(m, isodd(d2) ? -1 : 1)
    n::Int = 2(d - 1) - abs(m)
    get_j(m, n)
end

function j_to_fringe(j::Int)
    @domain_check_j
    mn = get_mn(j)
    mn == (0, 12) && return 37
    m, n = mn
    μ = abs(m)
    fringe = ((μ + n) ÷ 2 + 1)^2 - 2μ + fld(1 - sign(m), 2)
    @domain(fringe ∈ 1:36,
        """
        \nPolynomial does not have a Fringe representation.
        See Zernike.valid_fringes for valid polynomial indices.
        """,
        j
    )
    return fringe
end

function standardize(noll::Noll)
    (; v) = noll
    p = [noll_to_j(i) + 1 for i in eachindex(v)]
    u = zeros(maximum(p))
    u[p] = v
    return u
end

function standardize(fringe::Fringe)
    (; v) = fringe
    j = fringe_to_j.(eachindex(v))
    n_max = get_n(maximum(j))
    j_max = get_j(n_max)
    a = zeros(eltype(v), j_max + 1)
    a[j.+1] = v ./ N(j)
    return a
end

# radial order
get_n(j::Int) = ceil(Int, (-3 + √(9 + 8j)) / 2)

# azimuthal frequency
get_m(j::Int, n::Int) = 2j - (n + 2)n

get_m(j::Int) = get_mn(j)[1]

# ISO / ANSI / OSA standard single mode-ordering index
function get_j(m::Int, n::Int)
    μ = abs(m)
    @domain_check_mn
    return ((n + 2)n + m) ÷ 2
end

get_j((m, n)) = get_j(m, n)

get_j(n_max::Int) = get_j(n_max, n_max)

function get_mn(j::Int)
    @domain_check_j
    n = get_n(j)
    m = get_m(j, n)
    return m, n
end

radial_n_max(μ, a) = μ + 2 * (length(a) - 1)

mnv(v) = hcat(stack(get_mn(j) for j ∈ 0:length(v)-1; dims = 1), Vector{Number}(v))

to_i(m::Int, n::Int) = get_j(m, n) + 1
