"""
    Z(m, n)
    Z(j)

Plot a Zernike polynomial of azimuthal order `m` and radial degree `n`.

The single index `j` begins at zero and follows the ANSI Z80.28-2004 / ISO 24157:2008 / Optica (OSA) standard.

Returns three values contained within a Zernike.Output type, with fields:

1. `fig`: the Makie figure;
2. `coeffs`: vector of radial polynomial coefficients;
3. `latex`: LaTeX string of the Zernike polynomial.

These can also be accessed through indexing and regular non-property destructuring.

The coefficients belong to terms with exponent `n - 2(i - 1)` where `i` is the vector's index.

The radial polynomial coefficients are computed using a fast and accurate algorithm suitable for high orders; it is based on a recursive relation presented by Honarvar & Paramesran (2013) doi:10.1364/OL.38.002487.

See also [`W`](@ref), [`P`](@ref).

----

# Positional argument options:

    Z(m, n, ::Type{Model})

Return the Zernike polynomial function `Z(ρ, θ)` corresponding to indices `m` and `n`.

# Keyword argument options:

    Z(m, n; [finesse::Int = 100])

`finesse`: `{1 ≤ finesse ≤ 100}`: multiplicative factor determining the size of the plotted matrix; the total number of elements is capped at 1 million.
"""
Z

"""
    W(ρ, θ, OPD, n_max)

Fit wavefront errors up to order n_max.

Estimates wavefront error by expressing optical aberrations as a linear combination of weighted Zernike polynomials using a linear least squares method. The accuracy of this type of wavefront reconstruction represented as an expanded series depends upon a sufficiently sampled phase field and a suitable choice of the fitting order n_max.

# Main arguments

`ρ`, `θ`, and `OPD` must be floating-point vectors of equal length; at each specific index the values are elements of an ordered triple over the exit pupil.

* `ρ`: normalized radial exit pupil position variable `{0 ≤ ρ ≤ 1}`;
* `θ`: angular exit pupil variable in radians `(mod 2π)`, defined positive counter-clockwise from the horizontal x-axis;
* `OPD`: measured optical path difference in waves;
* `n_max`: maximum radial degree to fit to.

# Return values

Returns four values contained within a WavefrontOutput type, with fields:

1. `a`: vector of named tuples containing the Zernike polynomial indices and the corresponding expansion coefficients rounded according to `precision`;
2. `v`: full vector of Zernike wavefront error expansion coefficients;
3. `metrics`: named 3-tuple with the peak-to-valley error, RMS wavefront error, and Strehl ratio;
4. `fig`: the plotted Makie figure.

These can also be accessed through indexing and regular non-property destructuring.

See also [`Z`](@ref), [`P`](@ref).

----

    W(ρ, θ, OPD, orders::Vector{Tuple{Int, Int}})

Fit wavefront errors to specific Zernike polynomials specified in `orders` containing Zernike `(m, n)` tuples.

----

    W(OPD, fit_to; options...)

Fitting method accepting a floating-point matrix of phase data _uniformly_ produced in a polar coordinate system over the pupil.

The matrix is expected to be a polar grid of regularly spaced periodic samples with the first element referring to the value at the origin. The first axis of the matrix (the rows) must correspond to the angular variable `θ` while the second axis (the columns) must correspond to the radial variable `ρ`.

`fit_to` can be either `n_max::Int` or `orders::Vector{Tuple{Int, Int}}`.

----

    W(ρ::Vector, θ::Vector, OPD::Matrix, fit_to; options...)

Fitting method accepting coordinate vectors and a floating-point matrix of corresponding phase data produced in a polar coordinate system over the pupil under the aforementioned dimensional ordering assumption. This method does not assume equally spaced samples.

----

    W(x, y, OPD; fit_to, options...)

Fitting method accepting normalized Cartesian coordinate data.

----

# Positional argument options:

    W(ρ, θ, OPD, n_max, ::Type{Model})

Return the wavefront error function `ΔW(ρ, θ)` corresponding to an `n_max` fit.

# Keyword argument options:

    W(ρ, θ, OPD, n_max; [precision = 3], [finesse::Int])

* `precision`: number of digits to use after the decimal point in computing the expansion coefficients. Results will be rounded according to this precision and any polynomials with zero-valued coefficients will be ignored when pulling in the Zernike functions while constructing the composite wavefront error; this means lower precision values yield faster results.

\u2063\u2063\u2063\u2063 If you want full 64-bit floating-point precision use `precision = "full"`.

* `finesse`: `{1 ≤ finesse ≤ 100}`: multiplicative factor determining the size of the plotted matrix; the total number of elements is capped at 2^20 (~ 1 million).
"""
W

"""
    P(v::Vector{T}, ε::T, [δ::Complex{T}], [ϕ::T], [ω::Tuple{T,T}]) where T <: Float64

Compute a new set of Zernike wavefront error expansion coefficients under a given set of transformation factors and plot the result.

Available transformations are scaling, translation, & rotation for circular and elliptical exit pupils. These are essentially coordinate transformations in the pupil plane over the wavefront map.

# Main arguments

* `v`: vector of full Zernike expansion coefficients ordered in accordance with the ANSI / OSA single index standard. This is the `v` vector returned by the wavefront error fitting function `W(ρ, θ, OPD, n_max)`.
* `ε`: scaling factor `{0 ≤ ε ≤ 1}`;
* `δ`: translational complex coordinates (displacement of the pupil center in the complex plane);
* `ϕ`: rotation of the pupil in radians `(mod 2π)`, defined positive counter-clockwise from the horizontal x-axis;
* `ω`: elliptical pupil transform parameters; 2-tuple where `ω[1]` is the ratio of the minor radius to the major radius of the ellipse and `ω[2]` is the angle defined positive counter-clockwise from the horizontal coordinate axis of the exit pupil to the minor axis of the ellipse.

The order the transformations are applied is:\\
scaling --> translation --> rotation --> elliptical transform.

See also [`Z`](@ref), [`W`](@ref).

----

# Positional argument options:

    P(v, ε, [δ], [ϕ], [ω], ::Type{Model})

Return the wavefront error function `ΔW(ρ, θ)` corresponding to the input transform parameters.

# Keyword argument options:

    P(v, ε, [δ], [ϕ], [ω]; [precision = 3], [finesse::Int])

* `precision`: number of digits to use after the decimal point in computing the expansion coefficients. Results will be rounded according to this precision and any polynomials with zero-valued coefficients will be ignored when pulling in the Zernike functions while constructing the composite wavefront error; this means lower precision values yield faster results.

\u2063\u2063\u2063\u2063 If you want full 64-bit floating-point precision use `precision = "full"`.

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
P

"""
`Zernike.Polynomial`

Callable type: function `Z(ρ, θ)` bound to a given set of Zernike indices `m` and `n`.

Fields:

* `inds`: named tuple containing the Zernike polynomial indices;
* `N`: normalization factor;
* `R`: RadialPolynomial callable type: function `R(ρ)`;
* `M`: Sinusoid callable type: function `M(θ)`.

See also [`Zernike.WavefrontError`](@ref).
"""
Polynomial

"""
`Zernike.WavefrontError`

Callable type: function `ΔW(ρ, θ)` bound to a given set of Zernike polynomial functions `Zᵢ(ρ, θ)` and their corresponding expansion coefficients `aᵢ`.

Specifically, `ΔW(ρ, θ)` = `∑aᵢZᵢ(ρ, θ)`

Fields:

* `i`: vector of named tuples containing the Zernike polynomial indices and the corresponding expansion coefficients;
* `v`: full vector of the unfiltered full-precision standardized expansion coefficients up to `n_max`;
* `n_max`: maximum radial degree fit to;
* `fit_to`: vector of `(m, n)` tuples specifying the polynomials used for the fit;
* `a`: vector of the Zernike expansion coefficients corresponding to each polynomial present;
* `Z`: the respective Zernike polynomial functions.

The `fit_to` field is an empty vector if the default full range up to `n_max` (`0:j_max`) was used with no `orders` specified. Note that these orders could differ from the polynomials determined after the fit; they are simply what was passed to the fitting function and may refer to polynomials not present in the reconstruction if after filtering the corresponding coefficients are zero.

See also [`Zernike.Polynomial`](@ref).
"""
WavefrontError

"""
    noll_to_j(noll::Int)

Convert Noll indices to ANSI standard indices.

See also [`fringe_to_j`](@ref), [`standardize!`](@ref), [`standardize`](@ref).
"""
noll_to_j

"""
    standardize!(noll::Vector)

Re-order a Noll specified Zernike expansion coefficient vector according to the ANSI standard.

This requires a full ordered vector up to n_max.

See also [`standardize`](@ref), [`noll_to_j`](@ref), [`fringe_to_j`](@ref).
"""
standardize!

"""
    fringe_to_j(fringe::Int)

Convert Fringe indices to ANSI standard indices.

Only indices 1:37 are valid.

See also [`noll_to_j`](@ref), [`standardize`](@ref), [`standardize!`](@ref).
"""
fringe_to_j

"""
    standardize(fringe::Vector)

Format a Fringe specified Zernike expansion coefficient vector according to the ANSI standard.

This function expects unnormalized coefficients; the input coefficients will be re-ordered and normalized in line with the orthonormal standard. As Fringe is a 37 polynomial subset of the full set of Zernike polynomials any coefficients in the standard order missing a counterpart in the input vector will be set to zero.

See also [`standardize!`](@ref), [`fringe_to_j`](@ref), [`noll_to_j`](@ref).

----

    standardize(v_sub::Vector, orders::Vector{Tuple{Int, Int}})

Pad a subset Zernike expansion coefficient vector to the full standard length up to `n_max` (`1:j_max+1`).

The tuples in `orders` must be of the form `(m, n)` associated with the respective coefficients at each index in `v_sub`.
"""
standardize

"""
`plotconfig`

`Zernike` plot settings.

# Fields / Options:

* `size`::**Tuple{Float64, Float64}**: window size (DPI scaled resolution);
* `fontsize`::**Float64**: text size;
* `colormap`::**Symbol**: Default: `:oslo`;
* `focus_on_show`::**Bool**: whether the window is focused on generation (default: `true`).

There are two additional properties which trigger a settings refresh:

* `plotconfig.reset = true` will reset all of the settings to their defaults;
* `plotconfig.resize = true` will reset only the `size` and `fontsize` settings. This is useful if your primary monitor changes or you want to return to the automatically determined values.

See also [`zplot`](@ref)
"""
plotconfig

"""
    zplot(args..., kwargs...)

Plot `Polynomial` and `WavefrontError` input function types as well as quantized wavefront errors; for the latter `args...` must be a collection of discretized ρ, θ, ΔWp objects where the radial variables refer to either ranges or vectors and the wavefront error is a matrix.

# Keyword arguments:

* `size`::**Tuple{Float64, Float64}**: window size (DPI scaled resolution);
* `fontsize`::**Float64**: text size;
* `colormap`::**Symbol**: Default: `:oslo`;
* `focus_on_show`::**Bool**: whether the window is focused on generation (default: `true`);
* `window`::**String**: window title;
* `plot`::**Union{String, LaTeXString}: plot title;
* `m`::**Int**: azimuthal order (used to determine matrix size);
* `n`::**Int**: radial order (used to determine matrix size);
* `finesse`::**Int**: `{1 ≤ finesse ≤ 100}`: (used to determine matrix size);
* `high_order`::**Bool**: whether to apply a logarithmic transform (default: `false`).

See also [`plotconfig`](@ref)
"""
zplot
