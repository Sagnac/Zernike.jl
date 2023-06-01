# Zernike.jl

Generates Zernike polynomials, models wavefront errors, and plots them in Makie.

![Zernike.jl](images/image.png)

This package can be added from the Julia REPL by:
```
using Pkg
Pkg.add(url="https://github.com/Sagnac/Zernike.jl")
```
or entering the package mode by pressing `]` and entering
```
add https://github.com/Sagnac/Zernike.jl
```

The package provides 3 main functions for modelling Zernike polynomials and wavefront errors:

* `Z(m, n)`: Generates a Zernike polynomial, prints its symbolic representation, and plots it in Makie;

* `W(ρ, θ, OPD, n_max)`: Fits wavefront errors up to radial order `n_max` given an input set of data over the pupil, returns the Zernike expansion coefficients & various metrics, and plots the modelled wavefront error using Makie;

* `S(ε, v)`: Aperture scaling function which takes a scaling factor and a vector of Zernike expansion coefficients and returns the new set of expansion coefficients over the smaller pupil; the clipped wavefront error is also plotted using Makie.

----

## `Z(m, n)` | `Z(j)`

Generates a Zernike polynomial.

* `m`: azimuthal order;
* `n`: radial degree;
* `j`: ANSI Z80.28-2004 / ISO 24157:2008 / Optica (OSA) standard single-mode ordering index.

Returns three values contained within a Zernike.Output type, with fields:

1. `fig`: the Makie figure;
2. `coeffs`: vector of radial polynomial coefficients;
3. `latex`: LaTeX formatted string of the Zernike polynomial.

The coefficients belong to terms with exponent `n - 2(i - 1)` where `i` is the vector's index.

The radial polynomial coefficients are computed using a fast and accurate algorithm suitable for high orders; it is based on a recursive relation presented by [Honarvar & Paramesran (2013)](https://doi.org/10.1364/OL.38.002487).

----

## `W(ρ, θ, OPD, n_max)`

Estimates wavefront error by expressing optical aberrations as a linear combination of weighted Zernike polynomials using a linear least squares method. This representation as an expanded series is approximate, especially if the phase field is not sufficiently sampled.

`ρ`, `θ`, and `OPD` must be vectors of equal length; at each specific index the values are elements of an ordered triple over the exit pupil.

* `ρ`: normalized radial exit pupil variable `{0 ≤ ρ ≤ 1}`;
* `θ`: angular exit pupil variable defined counter-clockwise from the horizontal x-axis `{0 ≤ θ ≤ 2π}`;
* `OPD`: measured optical path difference in waves;
* `n_max`: maximum radial degree to fit to.

You can input normalized Cartesian coordinates using the 3-positional argument method:<br>
`W(x, y, OPD; n_max, options...)`.

The phase profile data can also be input as a 3 column matrix.

Returns four values contained within a WavefrontOutput type, with fields:

1. `a`: vector of named tuples containing the Zernike polynomial indices and the corresponding expansion coefficients rounded according to `precision`;
2. `v`: full vector of Zernike expansion coefficients;
3. `metrics`: 3-tuple with the peak-to-valley error, RMS wavefront error, and Strehl ratio;
4. `fig`: the plotted Makie figure.

----

## `S(ε, v)`

Pupil scaling function which computes a new set of Zernike expansion coefficients under aperture down-scaling and plots the result.

* `ε::Float64`: scaling factor `{0 ≤ ε ≤ 1}`;
* `v::Vector{Float64}`: vector of full Zernike expansion coefficients ordered in accordance with the ANSI / OSA single index standard. This is the `v` vector returned by the wavefront error fitting function `W(ρ, θ, OPD, n_max)`.

`ε` = `r₂/r₁` where `r₂` is the new radius, `r₁` the original

In particular the radial variable corresponding to the rescaled exit pupil is normalized such that:<br>
`ρ` = `r/r₂`; `{0 ≤ ρ ≤ 1}`<br>
`r`: radial pupil position, `r₂`: max. radius<br>
`W₂(ρ₂)` = `W₁(ερ₂)`

The rescaled expansion coefficients are computed using a fast and accurate algorithm suitable for high orders; it is based on a formula presented by [Janssen & Dirksen (2007)](https://doi.org/10.2971/jeos.2007.07012).

----

## Options

There are 2 options you can vary using keyword arguments. All 3 main functions support:

* `scale`: `{1 ≤ scale ≤ 100}`: multiplicative factor determining the size of the plotted matrix; the total number of elements is capped at 1 million.

Default: `100` (for `Z`, proportionally scaled according to the number of polynomials for the wavefront errors).

In creating the plot matrix the step size / length of the variable ranges is automatically chosen such that aliasing is avoided up to order ~317 radially and ~499 azimuthally. The `scale` parameter controls how fine the granularity is subsequently at the expense of performance.

Additionally, the wavefront error functions `W(ρ, θ, OPD, n_max)` and `S(ε, v)` support:

* `precision`: number of digits to use after the decimal point in computing the expansion coefficients. Results will be rounded according to this precision and any polynomials with zero-valued coefficients will be ignored when pulling in the Zernike functions while constructing the composite wavefront error; this means lower precision values yield faster results.

If you want full 64-bit floating-point precision use `precision = "full"`.

----

## Model functions

There exists a special method dispatch which avoids plotting and instead returns the `(ρ, θ)` functions as essentially closures. This is done by calling the 3 main functions with the `Model` constructor as the last positional argument. The pupil can then be evaluated using these functions with polar coordinates:

```
Z40 = Z(0, 4, Model())
Z40(0.7, π/4)
```

For the wavefront errors this is equivalent to `ΔW(ρ, θ)` = `∑aᵢZᵢ(ρ, θ)` where `aᵢ` and `Zᵢ` were determined from the fitting process according to `precision`.

## Additional Notes

* The output types for the 3 main functions can also be accessed by indexing them and regular destructuring in addition to property destructuring and getting the fields directly.

* The Zernike polynomials are currently only valid up to degree ~812 at which point the maximum coefficient approaches the maximum for double-precision floating-point numbers (~1e308).

* If you're interested in precompiling the package into a system image in order to speed up load times please see the [precompile directory](precompile) (at the moment PrecompileTools or the like is not used).

* If you're interested in only the full vector of Zernike expansion coefficients obtained through the least squares fit and want to avoid computing extra values and plotting the results you can call:
```
Zernike.Wf(ρ, θ, OPD, n_max)[1]
```
Similarly you can do this for the radial polynomial coefficients and the NA scaled wavefront error expansion coefficients by importing the functions `Φ` and `Π`, respectively.
