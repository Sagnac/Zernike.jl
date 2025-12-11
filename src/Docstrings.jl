"""
    Z(m, n)

Return the Zernike `Polynomial` function `Z(ρ, θ)` corresponding to indices `m` and `n`.

See also: [`W`](@ref), [`Y`](@ref).
"""
Z

"""
    W(ρ, θ, OPD, n_max)

Return the `Wavefront` function `ΔW(ρ, θ)` corresponding to an `n_max` fit.

----

    W(∂x::Vector{Float64}, ∂y::Vector{Float64}; [normalized = true])

Compute the original wavefront error coefficients given the Zernike expansion coefficient vectors of the partial derivatives of the wavefront error.

`normalized` refers to whether the derivatives are associated with unnormalized Zernike polynomials.

----

    W(∂x::Wavefront, ∂y::Wavefront)

Equivalent to the above, but for input wavefronts. Useful if the derivative data was fit into a Zernike basis using `W`. Returns a `Wavefront`.

See also: [`Wavefront`](@ref), [`Z`](@ref), [`Y`](@ref).
"""
W

"""
    Y(v, ε, [δ], [ϕ], [ω])

Return the `Wavefront` function `ΔW(ρ, θ)` corresponding to the input transform parameters.

See also: [`Z`](@ref), [`W`](@ref).
"""
Y

"""
`Zernike.Polynomial`

Callable type: function `Z(ρ, θ)` bound to a given set of Zernike indices `m` and `n`.

The single argument method `Zₘₙ(ρ)` will radially evaluate the polynomial with angle zero.

Fields:

* `inds`: named tuple containing the Zernike polynomial indices;
* `N`: normalization factor;
* `R`: `RadialPolynomial` callable type: function `R(ρ)`;
* `M`: `Harmonic` callable type: function `M(θ)`.

This type can be indexed (zero-based) to return a specific radial coefficient corresponding to the term with exponent `i`. Calling `getindex` without an explicit index will return the full vector of coefficients.

See also: [`Zernike.Wavefront`](@ref).
"""
Polynomial

"""
`Zernike.Wavefront`

Callable type: function `ΔW(ρ, θ)` bound to a given set of Zernike polynomial functions `Zᵢ(ρ, θ)` and their corresponding expansion coefficients `aᵢ`.

Specifically, `ΔW(ρ, θ)` = `∑aᵢZᵢ(ρ, θ)`

The single argument method `ΔW(ρ)` will radially evaluate the polynomials with angle zero.

Fields:

* `recap`: vector of named tuples containing the Zernike polynomial indices and the corresponding expansion coefficients;
* `v`: full vector of the unfiltered full-precision standardized expansion coefficients up to `n_max`;
* `n_max`: maximum radial degree fit to;
* `fit_to`: vector of `(m, n)` tuples specifying the polynomials used for the fit;
* `a`: vector of the Zernike expansion coefficients corresponding to each polynomial present;
* `Z`: the respective Zernike polynomial functions;
* `precision`: the precision with which the `a` & `Z` values were determined;
* `ssr`: the sum of the squared residuals from the fit.

The `fit_to` field is an empty vector if the default full range up to `n_max` (`0:j_max`) was used with no `orders` specified. Note that these orders could differ from the polynomials determined after the fit; they are simply what was passed to the fitting function and may refer to polynomials not present in the reconstruction if after filtering the corresponding coefficients are zero.

This type can be indexed (zero-based) to return a specific Zernike expansion coefficient corresponding to the Zernike polynomial of index `j`. Calling `getindex` without an explicit index will return the full vector of coefficients.

See also: [`Zernike.Polynomial`](@ref).

----

    Wavefront(G::Gradient{Polynomial})

Return the Zernike coefficients of the gradient as a `Wavefront`.

See also: [`grad`](@ref).

----

    Wavefront(L::Laplacian)

Return the Zernike coefficients of the Laplacian as a `Wavefront`.

See also: [`lap`](@ref).

----

    Wavefront{RadialPolynomial}(m::Int, a::Vector)

Construct a radial series in a Wavefront from an azimuthal index and a coefficient vector associated with a sequence of radial polynomials.

----

    Wavefront{RadialPolynomial}(m::Int, n::Vector{Int}, a::Vector)

Construct a Wavefront radial polynomial series from an azimuthal index, a vector of radial indices, and a coefficient vector associated with them; the indices are subject to the usual restrictions for Zernike polynomials.

See also [`*`](@ref).

----

    Wavefront(aberr::Aberration)

Convert a set of Seidel aberrations to Zernike polynomials.
"""
Wavefront

"""
    get_j(m::Int, n::Int)

Return the single mode-ordering index `j` corresponding to azimuthal & radial indices `(m, n)`.

----

    get_j(n_max::Int)

Return the single mode-ordering index `j` corresponding to the maximum radial index `n_max`; equivalent to `get_j(n_max, n_max)`.

See also: [`get_mn`](@ref), [`noll_to_j`](@ref), [`j_to_noll`](@ref), [`fringe_to_j`](@ref), [`j_to_fringe`](@ref), [`standardize`](@ref).
"""
get_j

"""
    get_mn(j::Int)

Return the azimuthal & radial indices `(m, n)` given the single mode-ordering index `j`.

See also: [`get_j`](@ref), [`noll_to_j`](@ref), [`j_to_noll`](@ref), [`fringe_to_j`](@ref), [`j_to_fringe`](@ref), [`standardize`](@ref).
"""
get_mn

"""
    get_m(j::Int)

Return the azimuthal index `m` given the single mode-ordering index `j`.

See also: [`get_n`](@ref), [`get_mn`](@ref), [`get_j`](@ref), [`noll_to_j`](@ref), [`j_to_noll`](@ref), [`fringe_to_j`](@ref), [`j_to_fringe`](@ref), [`standardize`](@ref).
"""
get_m

"""
    get_n(j::Int)

Return the radial index `n` given the single mode-ordering index `j`.

See also: [`get_m`](@ref), [`get_mn`](@ref), [`get_j`](@ref), [`noll_to_j`](@ref), [`j_to_noll`](@ref), [`fringe_to_j`](@ref), [`j_to_fringe`](@ref), [`standardize`](@ref).
"""
get_n

"""
    noll_to_j(noll::Int)

Convert Noll indices to ANSI standard indices.

See also: [`j_to_noll`](@ref), [`fringe_to_j`](@ref), [`j_to_fringe`](@ref), [`get_j`](@ref), [`get_mn`](@ref), [`standardize`](@ref).
"""
noll_to_j

"""
    j_to_noll(j::Int)

Convert ANSI standard indices to Noll indices.

See also: [`noll_to_j`](@ref), [`fringe_to_j`](@ref), [`j_to_fringe`](@ref), [`get_j`](@ref), [`get_mn`](@ref), [`standardize`](@ref).
"""
j_to_noll

"""
    fringe_to_j(fringe::Int)

Convert Fringe indices to ANSI standard indices.

Only indices 1:37 are valid.

See also: [`j_to_fringe`](@ref), [`noll_to_j`](@ref), [`j_to_noll`](@ref), [`standardize`](@ref), [`get_j`](@ref), [`get_mn`](@ref).
"""
fringe_to_j

"""
    j_to_fringe(j::Int)

Convert ANSI standard indices to Fringe indices.

Call `fringe_to_j.(1:37)` to return valid indices.

See also: [`fringe_to_j`](@ref), [`noll_to_j`](@ref), [`j_to_noll`](@ref), [`standardize`](@ref), [`get_j`](@ref), [`get_mn`](@ref).
"""
j_to_fringe

"""
    standardize(noll::Noll)
    standardize(fringe::Fringe)

Format a `Noll` or `Fringe` specified Zernike expansion coefficient vector according to the ANSI standard.

Floating-point coefficient vectors need to be wrapped in the index types (e.g. `standardize(Fringe(v))`).

The `Fringe` method expects unnormalized coefficients; the input coefficients will be re-ordered and normalized in line with the orthonormal standard. As Fringe is a 37 polynomial subset of the full set of Zernike polynomials any coefficients in the standard order missing a counterpart in the input vector will be set to zero.

See also: [`noll_to_j`](@ref), [`j_to_noll`](@ref), [`fringe_to_j`](@ref), [`j_to_fringe`](@ref), [`get_j`](@ref), [`get_mn`](@ref).

----

    standardize(v_sub::FloatVec, [j::AbstractVector{Int}])
    standardize(v_sub::Vector, orders::Vector{Tuple{Int, Int}})
    standardize(W::Wavefront)

Pad a subset Zernike expansion coefficient vector to the full standard length up to `n_max` (`1:j_max+1`).

The tuples in `orders` must be of the form `(m, n)` associated with the respective coefficients at each index in `v_sub`.

`j` is a vector of single-mode ordering indices associated with the coefficients; if this is not supplied the coefficients will be assumed to be in order (`0:j`).

The `Wavefront` method pads the `W.a` coefficient vector.

"""
standardize

"""
    Standard(v::Vector{Float64})

Wraps a standard ANSI / ISO vector of Zernike polynomial single indices for inter-conversions, viz. `Noll(s::Standard)` & `Fringe(s::Standard)`.

See also: [`standardize`](@ref).
"""
Standard

"""
    OTF(ΔW::Wavefront)

Return the Optical Transfer Function as a complex valued matrix for the input [`Wavefront`](@ref).

See also: [`MTF`](@ref).
"""
OTF

"""
    MTF(ΔW::Wavefront)

Return the Modulation Transfer Function as a real valued matrix for the input [`Wavefront`](@ref).

See also: [`OTF`](@ref), [`PSF`](@ref).
"""
MTF

"""
    PSF(ΔW::Wavefront)

Return the incoherent Point Spread Function as a real valued matrix for the input [`Wavefront`](@ref).

The output is normalized by the diffraction limited peak intensity i.e. the maximum is the Strehl ratio.

The keyword argument `s` (default `50.0`) determines the scale, i.e. the Fourier space step size is `0.25 / s`.

See also: See also: [`OTF`](@ref), [`MTF`](@ref).
"""
PSF

"""
    metrics(ΔW::Wavefront)

Compute wavefront error metrics. Returns a named 3-tuple with the peak-to-valley error, RMS wavefront error, and Strehl ratio.
"""
metrics

"""
    radial_coefficients(m::Int, n::Int, T::Type{<:Number} = Float64)

Compute `Zernike` radial polynomial coefficients as type `T` using an algorithm based on `Honarvar & Paramesran's` recursive relation suitable for high orders.

The coefficients in the vector correspond to terms with powers in ascending order for a Zernike polynomial with indices `m` & `n` subject to the following requirements:

    n ≥ 0
    |m| ≤ n
    n ≡ m (mod 2).

See also: [`wavefront_coefficients`](@ref), [`transform_coefficients`](@ref).
"""
radial_coefficients

"""
    wavefront_coefficients(ρ, θ, OPD, n_max)

Returns a 2-tuple with the full vector of Zernike expansion coefficients obtained through the least squares fit and the corresponding sum of the squared residuals.

See also: [`radial_coefficients`](@ref), [`transform_coefficients`](@ref).
"""
wavefront_coefficients

"""
    transform_coefficients(v, ε, δ, ϕ, ω)

Directly compute `Zernike` wavefront error expansion coefficients under pupil transformations. The argument types are the same as in `transform`.

Returns a 2-tuple with the new coefficient vector and order `n_max`.

See also: [`radial_coefficients`](@ref), [`wavefront_coefficients`](@ref).
"""
transform_coefficients

"""
    S(v, ε; precision)

Scale the pupil over a wavefront using an algorithm based on `Janssen & Dirksen's` formula and return a new `Wavefront`.

`v` is the set of Zernike wavefront error expansion coefficients and `ε` is the scaling factor.

See also: [`Y`](@ref).
"""
S

"""
    map_phase(ρ, θ, OPD)

Reverse dimensional coordinate transform with respect to the main wavefront error method. Returns the OPD as a matrix along with the corresponding unique coordinate vectors. Assumes uniform sampling.

"""
map_phase

"""
    reduce_wave(W::Wavefront, precision::Int)

Reduces `Wavefront` precision.
"""
reduce_wave

"""
    sieve(v::Vector{Float64}, threshold::Float64)

Zero out any elements lower than the threshold.

----

    sieve(a::Vector)

Return `Zernike` non-sequential indices and the corresponding wavefront expansion coefficients for any non-zero coefficients in the input vector.
"""
sieve

"""
    format_strings(Z::AbstractPolynomial)

Return a 3-tuple with the index formatted `LaTeX` variable name, the full `LaTeX` string equation, and the `Unicode` string representation of the polynomial. Calling the polynomial with the `String` type is a shortcut which yields just the symbolic `Unicode` string.

```jldoctest
julia> Z(-8, 8)(String)
"√(18)ρ⁸sin(8θ)"
```

----

    format_strings(ΔW::Wavefront)

Return a (possibly truncated) symbolic representation of the wavefront error in a Zernike basis as a `LaTeXString`.

See also: [`print_strings`](@ref).
"""
format_strings

"""
    print_strings([io::IO], j::AbstractVector; [limit = true])

Print the Unicode string representations of select Zernike polynomials; `j` here can be a range / vector of `Zernike` single indices; the output stream defaults to `stdout` and the vertical output is limited by your `displaysize` unless you pass the keyword argument `limit = false`.

----

    print_strings(j_max::Int)

Print the Unicode string representations of the first `j + 1` Zernike polynomials from single index `0` to `j_max`.

See also: [`format_strings`](@ref).
"""
print_strings

"""
    Zernike.Gradient(Z::Polynomial)

Returns the gradient `∇Z(ρ, θ)` in a polar basis.

If called with a complex number `x + iy` this function returns the vector `[∂Z/∂x, ∂Z/∂y]` at the point `(x, y)` in Cartesian coordinates instead.

See also: [`grad`](@ref), [`lap`](@ref), [`derivatives`](@ref), [`Laplacian`](@ref).
"""
Gradient

"""
    Zernike.Laplacian(Z::Polynomial)

Returns the Laplacian `ΔZ(ρ, θ)`.

This function can be evaluated in Cartesian coordinates if passed a complex number instead (i.e. `∇²Z(x, y) ≡ ΔZ(xy::Complex)`).

See also: [`grad`](@ref), [`lap`](@ref), [`derivatives`](@ref), [`Gradient`](@ref).
"""
Laplacian

"""
    derivatives(Z::Polynomial, order::Int = 1)

Computes the nth order partial derivatives of `Z(ρ, θ)` and returns the two-tuple (`∂Z/∂ρ`, `∂Z/∂θ`).

See also: [`grad`](@ref), [`lap`](@ref), [`Gradient`](@ref), [`Laplacian`](@ref).
"""
derivatives

"""
    ∇(Z)
    grad(Z::Polynomial)

Return the gradient of the Zernike polynomial expressed as a `Wavefront` in Zernike polynomial expansion coefficients.

Equivalent to:

    Wavefront(g::Gradient{Polynomial})
    Wavefront(::Type{<:Gradient}, m, n; [normalize = true])
    Wavefront(::Type{<:Gradient}, j; [normalize = true])

with the last two methods allowing unnormalized input. Normalized means the polynomial is expressed as `Z(ρ, θ) = N * R(ρ) * M(θ)` with `N` being the normalization prefactor required so that `π`-normalized integration with respect to the areal measure of `Z²` over the unit disk yields unity .

----

    grad(m, n, ::Type{Matrix{Complex}})

Return a complex matrix encoding the partial derivatives of a complex unnormalized `Zernike` polynomial at indices `m` & `n`.

The first two columns of the matrix refer to the `x` partial derivatives of `Z(|m|, n)` & `Z(-|m|, n)`, respectively. Likewise, the last two columns refer to the `y` partial derivatives of `Z(|m|, n)` & `Z(-|m|, n)`.

----

    grad(m, n, ::Type{Vector{Complex}})

Convenience method which returns the relevant sign-dependent complex gradient from the above mentioned `grad(m, n, Matrix{Complex})` call.

----

    grad(m, n, ::Type{Vector{Real}}; [normalize = true])

Return the gradient of the Zernike polynomial as a tuple of real-valued Zernike coefficient vectors.

----

    ∇(W)
    grad(W::Wavefront)

Computes the gradient of the wavefront.

See also: [`Zernike.lap`](@ref), [`Zernike.Gradient`](@ref), [`Zernike.Laplacian`](@ref), [`derivatives`](@ref).
"""
grad

"""
    lap(Z::Polynomial)

Return the Laplacian of the Zernike polynomial expressed as a `Wavefront` in Zernike polynomial expansion coefficients.

Equivalent to:

    Wavefront(l::Laplacian)
    Wavefront(::Type{<:Laplacian}, m, n; [normalize = true])
    Wavefront(::Type{<:Laplacian}, j; [normalize = true])

with the last two methods allowing unnormalized input. Normalized means the polynomial is expressed as `Z(ρ, θ) = N * R(ρ) * M(θ)` with `N` being the normalization prefactor required so that `π`-normalized integration with respect to the areal measure of `Z²` over the unit disk yields unity .

----

    lap(m, n; [normalize = true])

Return the Laplacian of the Zernike polynomial as a real-valued coefficient vector.

See also: [`grad`](@ref), [`Gradient`](@ref), [`Laplacian`](@ref), [`derivatives`](@ref).
"""
lap

"""
    mnv(v)

`Zernike` tool which takes a vector `v` (e.g. a vector of single-index `j` integers or a real or complex floating point coefficient vector) and produces a (length(v) × 3) `Matrix{Number}` with columns corresponding to the indices `m`, `n`, & the vector `v`.

See also: [`get_n`](@ref), [`get_m`](@ref), [`get_mn`](@ref), [`get_j`](@ref).
"""
mnv

"""
    *(w1::Wavefront{T}, w2::Wavefront{T}; [threshold = 5.0E-324]) where T <: RadialPolynomial
    w1 * w2

Compute the radial product expansion of two `Zernike` radial polynomial sequences.

Returns the coefficients of the new series embedded in a Wavefront. Based on a paper by Cadiot et al. (2024) using Matrix Multiplication Transforms and their inverses to perform quadrature.

The coefficients will be filtered according to `threshold`, with any absolute values less than it being set to zero.

See also: [`Wavefront`](@ref).
"""
*

@doc (@doc grad) ∇
