"""
    Z(m, n)
    Z(j)

Plots a Zernike polynomial of azimuthal order `m` and radial degree `n`. The single index `j` begins at zero and follows the ANSI / OSA standard. Returns three values contained within a Zernike.Output type, with fields:

1. `fig`: the Makie figure;
2. `coeffs`: vector of radial polynomial coefficients;
3. `latex`: LaTeX string of the Zernike polynomial.

These can also be accessed through indexing and regular non-property destructuring.

The coefficients belong to terms with exponent `n - 2(i - 1)` where `i` is the vector's index.

The radial polynomial coefficients are computed using a fast and accurate algorithm suitable for high orders; it is based on a recursive relation presented by Honarvar & Paramesran (2013) doi:10.1364/OL.38.002487.

----

# Positional argument options:

    Z(m, n, Model())

Returns the Zernike polynomial function `Z(ρ, θ)` corresponding to indices `m` and `n`.

# Keyword argument options:

    Z(m, n; scale::Int = 100)

`scale`: `{1 ≤ scale ≤ 100}`: multiplicative factor determining the size of the plotted matrix; the total number of elements is capped at 1 million.
"""
Z

"""
    W(ρ, θ, OPD, n_max)
    W(data::Matrix, n_max; options...)

Estimates wavefront error by expressing optical aberrations as a linear combination of weighted Zernike polynomials using a linear least squares method. This representation as an expanded series is approximate, especially if the phase field is not sufficiently sampled.

`ρ`, `θ`, and `OPD` must be vectors of equal length; at each specific index the values are elements of an ordered triple over the exit pupil.

* `ρ`: normalized radial exit pupil variable `{0 ≤ ρ ≤ 1}`;
* `θ`: angular exit pupil variable in radians `(mod 2π)`, defined positive counter-clockwise from the horizontal x-axis;
* `OPD`: measured optical path difference in waves;
* `n_max`: maximum radial degree to fit to.

The phase profile data can also be input as a 3 column matrix.

Returns four values contained within a WavefrontOutput type, with fields:

1. `a`: vector of named tuples containing the Zernike polynomial indices and the corresponding expansion coefficients rounded according to `precision`;
2. `v`: full vector of Zernike expansion coefficients;
3. `metrics`: 3-tuple with the peak-to-valley error, RMS wavefront error, and Strehl ratio;
4. `fig`: the plotted Makie figure.

These can also be accessed through indexing and regular non-property destructuring.

----

    W(x, y, OPD; n_max, options...)

Method accepting normalized Cartesian coordinate data.

----

# Positional argument options:

    W(ρ, θ, OPD, n_max, Model())

Returns the wavefront error function `ΔW(ρ, θ)` corresponding to an `n_max` fit.

# Keyword argument options:

    W(ρ, θ, OPD, n_max; precision = 3, scale)

* `precision`: number of digits to use after the decimal point in computing the expansion coefficients. Results will be rounded according to this precision and any polynomials with zero-valued coefficients will be ignored when pulling in the Zernike functions while constructing the composite wavefront error; this means lower precision values yield faster results.

\u2063\u2063\u2063\u2063 If you want full 64-bit floating-point precision use `precision = "full"`.

* `scale`: `{1 ≤ scale ≤ 100}`: multiplicative factor determining the size of the plotted matrix; the total number of elements is capped at 1 million.
"""
W

"""
    S(ε::Float64, v::Vector{Float64})

Pupil scaling function which computes a new set of Zernike expansion coefficients under aperture down-scaling and plots the result.

* `ε`: scaling factor `{0 ≤ ε ≤ 1}`;
* `v`: vector of full Zernike expansion coefficients ordered in accordance with the ANSI / OSA single index standard. This is the `v` vector returned by the wavefront error fitting function `W(ρ, θ, OPD, n_max)`.

`ε` = `r₂/r₁` where `r₂` is the new radius, `r₁` the original

In particular the radial variable corresponding to the rescaled exit pupil is normalized such that:\\
`ρ` = `r/r₂`; `{0 ≤ ρ ≤ 1}`\\
`r`: radial pupil position, `r₂`: max. radius\\
`W₂(ρ₂)` = `W₁(ερ₂)`

The rescaled expansion coefficients are computed using a fast and accurate algorithm suitable for high orders; it is based on a formula presented by Janssen & Dirksen (2007) doi:10.2971/jeos.2007.07012.

# Positional argument options:

    S(ε, v, Model())

Returns the wavefront error function `ΔW(ρ, θ)` corresponding to an `ε` aperture scaling.

# Keyword argument options:

    S(ε, v; precision = 3, scale::Int)

* `precision`: number of digits to use after the decimal point in computing the expansion coefficients. Results will be rounded according to this precision and any polynomials with zero-valued coefficients will be ignored when pulling in the Zernike functions while constructing the composite wavefront error; this means lower precision values yield faster results.

\u2063\u2063\u2063\u2063 If you want full 64-bit floating-point precision use `precision = "full"`.

* `scale`: `{1 ≤ scale ≤ 100}`: multiplicative factor determining the size of the plotted matrix; the total number of elements is capped at 1 million.
"""
S

"""
`Zernike.Polynomial`

Callable type: function `Z(ρ, θ)` bound to a given set of Zernike indices `m` and `n`.

Fields:

1. `inds`: named tuple containing the Zernike polynomial indices;
2. `N`: normalization factor;
3. `R`: RadialPolynomial callable type: function `R(ρ)`;
4. `M`: Sinusoid callable type: function `M(θ)`.
"""
Polynomial

"""
`Zernike.WavefrontError`

Callable type: function `ΔW(ρ, θ)` bound to a given set of Zernike polynomial functions `Zᵢ(ρ, θ)` and their corresponding expansion coefficients `aᵢ`.

Specifically, `ΔW(ρ, θ)` = `∑aᵢZᵢ(ρ, θ)`

Fields:

1. `i`: vector of named tuples containing the Zernike polynomial indices and the corresponding expansion coefficients;
2. `n_max`: maximum radial degree fit to;
3. `a`: vector of the Zernike expansion coefficients used in the fit;
4. `Z`: the respective Zernike polynomial functions.
"""
WavefrontError
