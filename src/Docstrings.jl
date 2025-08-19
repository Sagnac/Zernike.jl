"""
    zernike(m, n)
    zernike(j)

Plot a Zernike polynomial of azimuthal order `m` and radial degree `n`.

The single index `j` begins at zero and follows the ANSI Z80.28-2004 / ISO 24157:2008 / Optica (OSA) standard.

Returns a `Zernike.Output` type which contains (among other things):

* `Z`: the `Polynomial` function `Z(ρ, θ)`;
* `fig`: the `Makie` figure;
* `coeffs`: vector of radial polynomial coefficients;
* `latex`: `LaTeX` string of the Zernike polynomial;
* `unicode`: `Unicode` string of the Zernike polynomial.

The coefficients belong to terms with exponent `n - 2(i - 1)` where `i` is the vector's index.

The radial polynomial coefficients are computed using a fast and accurate algorithm suitable for high orders; it is based on a recursive relation presented by Honarvar & Paramesran (2013) doi:10.1364/OL.38.002487.

See also: [`Z`](@ref), [`wavefront`](@ref), [`transform`](@ref), [`radial_coefficients`](@ref).

----

# Keyword argument options:

    zernike(m, n; [finesse::Int = 100])

`finesse`: `{1 ≤ finesse ≤ 100}`: multiplicative factor determining the size of the plotted matrix; the total number of elements is capped at 1 million.
"""
zernike

"""
    Z(m, n)

Return the Zernike `Polynomial` function `Z(ρ, θ)` corresponding to indices `m` and `n`.

See also: [`zernike`](@ref), [`W`](@ref), [`P`](@ref).
"""
Z

"""
    wavefront(ρ, θ, OPD, n_max)

Fit wavefront errors up to order `n_max`.

Estimates wavefront error by expressing optical aberrations as a linear combination of weighted Zernike polynomials using a linear least squares method. The accuracy of this type of wavefront reconstruction represented as an expanded series depends upon a sufficiently sampled phase field and a suitable choice of the fitting order `n_max`.

# Main arguments

`ρ`, `θ`, and `OPD` must be floating-point vectors of equal length; at each specific index the values are elements of an ordered triple over the exit pupil.

* `ρ`: normalized radial exit pupil position variable `{0 ≤ ρ ≤ 1}`;
* `θ`: angular exit pupil variable in radians `(mod 2π)`, defined positive counter-clockwise from the horizontal x-axis;
* `OPD`: measured optical path difference in waves;
* `n_max`: maximum radial degree to fit to.

# Return values

Returns seven values contained within a `WavefrontOutput` type, with fields:

* `recap`: vector of named tuples containing the Zernike polynomial indices and the corresponding expansion coefficients rounded according to `precision`;
* `v`: full vector of Zernike wavefront error expansion coefficients;
* `metrics`: named 3-tuple with the peak-to-valley error, RMS wavefront error, and Strehl ratio;
* `W`: the `WavefrontError` function `ΔW(ρ, θ)`;
* `fig`: the plotted `Makie` figure;
* `axis`: the plot axis;
* `plot`: the surface plot object.

See also: [`W`](@ref), [`zernike`](@ref), [`transform`](@ref), [`reconstruct`](@ref).

----

    wavefront(ρ, θ, OPD, orders::Vector{Tuple{Int, Int}})

Fit wavefront errors to specific Zernike polynomials specified in `orders` containing Zernike `(m, n)` tuples.

----

    wavefront(OPD, fit_to; options...)

Fitting method accepting a floating-point matrix of phase data _uniformly_ produced in a polar coordinate system over the pupil.

The matrix is expected to be a polar grid of regularly spaced periodic samples with the first element referring to the value at the origin. The first axis of the matrix (the rows) must correspond to the angular variable `θ` while the second axis (the columns) must correspond to the radial variable `ρ`.

`fit_to` can be either `n_max::Int` or `orders::Vector{Tuple{Int, Int}}`.

----

    wavefront(ρ::Vector, θ::Vector, OPD::Matrix, fit_to; options...)

Fitting method accepting coordinate vectors and a floating-point matrix of corresponding phase data produced in a polar coordinate system over the pupil under the aforementioned dimensional ordering assumption. This method does not assume equally spaced samples.

----

    wavefront(x, y, OPD; fit_to, options...)

Fitting method accepting normalized Cartesian coordinate data.

----

# Keyword argument options:

    wavefront(ρ, θ, OPD, n_max; [precision = 3], [finesse::Int])

* `precision`: number of digits to use after the decimal point in computing the expansion coefficients. Results will be rounded according to this precision and any polynomials with zero-valued coefficients will be ignored when pulling in the Zernike functions while constructing the composite wavefront error; this means lower precision values yield faster results.

* `finesse`: `{1 ≤ finesse ≤ 100}`: multiplicative factor determining the size of the plotted matrix; the total number of elements is capped at 2^20 (~ 1 million).
"""
wavefront

"""
    W(ρ, θ, OPD, n_max)

Return the `WavefrontError` function `ΔW(ρ, θ)` corresponding to an `n_max` fit.

See also: [`wavefront`](@ref), [`Z`](@ref), [`P`](@ref).
"""
W

"""
    transform(v::Vector{T}, ε::T, [δ::Complex{T}], [ϕ::T], [ω::Tuple{T,T}]) where T <: Float64

Compute a new set of Zernike wavefront error expansion coefficients under a given set of transformation factors and plot the result.

Available transformations are scaling, translation, & rotation for circular and elliptical exit pupils. These are essentially coordinate transformations in the pupil plane over the wavefront map.

# Main arguments

* `v`: vector of full Zernike expansion coefficients ordered in accordance with the ANSI / OSA single index standard. This is the `v` vector returned by `wavefront(ρ, θ, OPD, n_max)`;
* `ε`: scaling factor `{0 ≤ ε ≤ 1}`;
* `δ`: translational complex coordinates (displacement of the pupil center in the complex plane);
* `ϕ`: rotation of the pupil in radians `(mod 2π)`, defined positive counter-clockwise from the horizontal x-axis;
* `ω`: elliptical pupil transform parameters; 2-tuple where `ω[1]` is the ratio of the minor radius to the major radius of the ellipse and `ω[2]` is the angle defined positive counter-clockwise from the horizontal coordinate axis of the exit pupil to the minor axis of the ellipse.

The order the transformations are applied is:\\
scaling --> translation --> rotation --> elliptical transform.

See also: [`P`](@ref), [`zernike`](@ref), [`wavefront`](@ref), [`transform_coefficients`](@ref).

----

# Keyword argument options:

    transform(v, ε, [δ], [ϕ], [ω]; [precision = 3], [finesse::Int])

* `precision`: number of digits to use after the decimal point in computing the expansion coefficients. Results will be rounded according to this precision and any polynomials with zero-valued coefficients will be ignored when pulling in the Zernike functions while constructing the composite wavefront error; this means lower precision values yield faster results.

* `finesse`: `{1 ≤ finesse ≤ 100}`: multiplicative factor determining the size of the plotted matrix; the total number of elements is capped at 2^20 (~ 1 million).

# Extended help

`ε` = `r₂/r₁` where `r₂` is the new smaller radius, `r₁` the original

In particular the radial variable corresponding to the rescaled exit pupil is normalized such that:\\
`ρ` = `r/r₂`; `{0 ≤ ρ ≤ 1}`\\
`r`: radial pupil position, `r₂`: max. radius\\
`ΔW₂(ρ₂, θ)` = `ΔW₁(ερ₂, θ)`

For translation the shift must be within the bounds of the scaling applied such that:\\
`0.0 ≤ ε + |δ| ≤ 1.0`.

For elliptical pupils (usually the result of measuring the wavefront off-axis), the major radius is defined such that it equals the radius of the circle and so `ω[1]` is the fraction of the circular pupil covered by the minor radius (this is approximated well by a cosine projection factor for angles up to 40 degrees); `ω[2]` is then the direction of the stretching applied under transformation in converting the ellipse to a circle before fitting the expansion coefficients.

The transformed expansion coefficients are computed using a fast and accurate algorithm suitable for high orders; it is based on a formulation presented by Lundström & Unsbo (2007) doi:10.1364/JOSAA.24.000569.
"""
transform

"""
    P(v, ε, [δ], [ϕ], [ω])

Return the `WavefrontError` function `ΔW(ρ, θ)` corresponding to the input transform parameters.

See also: [`transform`](@ref), [`Z`](@ref), [`W`](@ref).
"""
P

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

See also: [`Zernike.WavefrontError`](@ref).
"""
Polynomial

"""
`Zernike.WavefrontError`

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
* `precision`: the precision with which the `a` & `Z` values were determined.

The `fit_to` field is an empty vector if the default full range up to `n_max` (`0:j_max`) was used with no `orders` specified. Note that these orders could differ from the polynomials determined after the fit; they are simply what was passed to the fitting function and may refer to polynomials not present in the reconstruction if after filtering the corresponding coefficients are zero.

This type can be indexed (zero-based) to return a specific Zernike expansion coefficient corresponding to the Zernike polynomial of index `j`. Calling `getindex` without an explicit index will return the full vector of coefficients.

See also: [`Zernike.Polynomial`](@ref).
"""
WavefrontError

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

    standardize(v_sub::FloatVec, j::AbstractVector{Int})
    standardize(v_sub::Vector, orders::Vector{Tuple{Int, Int}})
    standardize(W::WavefrontError)

Pad a subset Zernike expansion coefficient vector to the full standard length up to `n_max` (`1:j_max+1`).

The tuples in `orders` must be of the form `(m, n)` associated with the respective coefficients at each index in `v_sub`.

`j` is a vector of single-mode ordering indices associated with the coefficients.

The `WavefrontError` method pads the `W.a` coefficient vector.

"""
standardize

"""
    Standard(v::Vector{Float64})

Wraps a standard ANSI / ISO vector of Zernike polynomial single indices for inter-conversions, viz. `Noll(s::Standard)` & `Fringe(s::Standard)`.

See also: [`standardize`](@ref).
"""
Standard

"""
`Zernike` plot settings.

# Fields / Options:

* `size`::**Tuple{Float64, Float64}**: window size (DPI scaled resolution);
* `fontsize`::**Float64**: text size;
* `colormap`::**Symbol**: Default: `:oslo`;
* `focus_on_show`::**Bool**: whether the window is focused on generation (default: `true`).

There are two additional properties which trigger a settings refresh:

* `plotconfig.reset = true` will reset all of the settings to their defaults;
* `plotconfig.resize = true` will reset only the `size` and `fontsize` settings. This is useful if your primary monitor changes or you want to return to the automatically determined values.

See also: [`zplot`](@ref).
"""
plotconfig

"""
    zplot(P; kwargs...)
    zplot(ρ, θ, w; kwargs...)

Plot `Polynomial` and `WavefrontError` input function types `P` as well as quantized wavefront errors; for the latter the arguments must be a collection of discretized `ρ`, `θ`, `w` samples where the polar variable objects refer to either ranges or vectors and the wavefront error `w` is a matrix.

# Keyword arguments:

* `size`::**Tuple{Float64, Float64}**: window size (DPI scaled resolution);
* `fontsize`::**Float64**: text size;
* `colormap`::**Symbol**: Default: `:oslo`;
* `focus_on_show`::**Bool**: whether the window is focused on generation (default: `true`);
* `window_title`::**String**: window title;
* `plot_title`::**Union{String, LaTeXString}**: plot title;
* `m`::**Int**: azimuthal order (used to determine matrix size);
* `n`::**Int**: radial order (used to determine matrix size);
* `finesse`::**Int**: `{1 ≤ finesse ≤ 100}`: (used to determine matrix size);
* `high_order`::**Bool**: whether to apply a logarithmic transform (default: `false`).

Plots can be updated on demand by passing an `Observable` and changing its value.

For example:

```julia
w = Observable(WavefrontError([0.0, -1.0, 1.0]))
zplot(w)

# update
w[] = WavefrontError([0.0, 1.0, 1.0])
```

See also: [`plotconfig`](@ref).
"""
zplot

"""
    metrics(ΔW::WavefrontError)

Compute wavefront error metrics. Returns a named 3-tuple with the peak-to-valley error, RMS wavefront error, and Strehl ratio.
"""
metrics

"""
    radial_coefficients(m_max, n_max)

Compute the `Zernike` radial polynomial coefficients using an algorithm based on `Honarvar & Paramesran's` recursive relation suitable for high orders.

Returns a vector of coefficient vectors for each unique radial polynomial up to Zernike indices `m_max`, `n_max`. The coefficient vector elements correspond to radial degree `n` in ascending order.

`m_max` and `n_max` must be of type `Int` and must meet the following additional requirements:

    m_max ≥ 0
    n_max ≥ 0
    m_max ≤ n_max
    n_max ≡ m_max (mod 2).

See also: [`wavefront_coefficients`](@ref), [`transform_coefficients`](@ref).
"""
radial_coefficients

"""
    wavefront_coefficients(ρ, θ, OPD, n_max)

Returns the full vector of Zernike expansion coefficients obtained through the least squares fit.

See also: [`wavefront`](@ref), [`radial_coefficients`](@ref), [`transform_coefficients`](@ref).
"""
wavefront_coefficients

"""
    transform_coefficients(v, ε, δ, ϕ, ω)

Directly compute `Zernike` wavefront error expansion coefficients under pupil transformations. The argument types are the same as in `transform`.

Returns a 2-tuple with the new coefficient vector and order `n_max`.

See also: [`transform`](@ref), [`radial_coefficients`](@ref), [`wavefront_coefficients`](@ref).
"""
transform_coefficients

"""
    reconstruct(ρ::Vector, θ::Vector, OPD::Vector, fit_to::Union{Int, Vector{Tuple{Int, Int}}})
    reconstruct(ρ::Vector, θ::Vector, OPD::Matrix, fit_to)
    reconstruct(OPD::Matrix, fit_to)

Fit wavefront errors in terms of Zernike polynomials without computing extra results or plotting the wavefront error.

Returns a 2-tuple with the full vector of expansion coefficients and the corresponding Zernike polynomials.

See also: [`wavefront`](@ref).
"""
reconstruct

"""
    scale(v, ε; precision, finesse)

Scale the pupil over a wavefront using an algorithm based on `Janssen & Dirksen's` formula and plot the result.

`v` is the set of Zernike wavefront error expansion coefficients and `ε` is the scaling factor.

See also: [`transform`](@ref), [`J`](@ref).
"""
scale

"""
    J(v, ε; precision)

Scale the pupil over a wavefront using an algorithm based on `Janssen & Dirksen's` formula and return a new `WavefrontError`.

`v` is the set of Zernike wavefront error expansion coefficients and `ε` is the scaling factor.

See also: [`P`](@ref), [`scale`](@ref).
"""
J

"""
    map_phase(ρ, θ, OPD)

Reverse dimensional coordinate transform with respect to the main wavefront error method. Returns the OPD as a matrix along with the corresponding unique coordinate vectors. Assumes uniform sampling.

See also: [`wavefront`](@ref).
"""
map_phase

"""
    reduce_wave(W::WavefrontError, precision::Int)

Reduces `WavefrontError` precision.
"""
reduce_wave

"""
    sieve(v::Vector{Float64}, threshold::Float64)

Zero out any elements lower than the threshold.
"""
sieve

"""
    format_strings(m::Int, n::Int)

Return a 2-tuple with the `Unicode` and `LaTeX` string representations of the Zernike polynomial.

```jldoctest
julia> unicode, latex = Zernike.format_strings(-8, 8);

julia> unicode
"√(18)ρ⁸sin(8θ)"
```
"""
format_strings

"""
    Zernike.Gradient(Z::Polynomial)

Returns `∇Z(ρ, θ)`.
"""
Gradient

"""
    Zernike.derivatives(Z::Polynomial, order::Int = 1)

Computes the nth order partial derivatives of `Z(ρ, θ)` and returns the two-tuple (`∂Z/∂ρ`, `∂Z/∂θ`).
"""
derivatives
