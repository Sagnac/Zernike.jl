# v6.3.0

* `grad(W)` has replaced `derivatives(W)` for Wavefronts;
* `∇` is now an alias for `grad`;
* Type instabilities have been fixed for types with `getproperty` overloaded.

----

# v6.2.0

* Seidel aberration conversions have been added through `Wavefront(aberr::Aberration)`;
* The two output structs now hold the `FigureAxisPlot` rather than the individual plot elements.

----

# v6.1.0

* The `⋆` (`\star`) operator for radial wavefronts has been changed to the usual multiplication operator;
* Multiplication of radial wavefronts now accepts a `threshold` keyword argument which zeros out any coefficients less than its value;
* `Wavefront{RadialPolynomial}`s can now be constructed from an additional input set of radial indices.

----

# v6.0.0

## Major features

* Three new algorithms associated with derivatives, based on a 2014 paper by Janssen, have been added:
  * Zernike polynomial derivatives can now be expressed in a Zernike expansion basis using `Zernike.grad(Z::Polynomial)` or `Wavefront(g::Gradient)`;
  * Laplacians can be expressed in a Zernike basis using `Zernike.lap(Z::Polynomial)` or `Wavefront(L::Laplacian)`;
  * Wavefront error derivatives can now be used to compute the coefficients of the original wavefront using `W(∂x, ∂y)`.

* Derivatives of `Wavefront`s can now be taken by calling `derivatives(W)`.

* Zernike radial polynomial product expansions can be computed by using the the `⋆` (`\star`) operator between two constructed `Wavefront{RadialPolynomial}`s; the quadrature based algorithm is based on a 2024 paper by [Cadiot, et al.](https://arxiv.org/abs/2411.18361)

## Breaking changes

* Various types and functions have been renamed:
  * `WavefrontError` --> `Wavefront`;
  * `PartialDerivative` --> `Derivative`;
  * The pupil transform function `P` --> `Y`;
  * The pupil scaling function `J` --> `S`;

## Plotting changes

* The default plotting theme has been changed to a dark one;
* `zplot` now accepts specified polar coordinates for generic functor arguments;
* `zplot` now forwards any additional keyword arguments to `surface!`;
* Phase `Observables` now have a shortcut zero-argument `zplot` method;

## Miscellaneous

* A trimmed down base version without the `GLMakie` dependency and plotting is now available at the [`base`](/../../tree/base) branch;
* Cartesian coordinate wavefront error fitting methods are now available for input phase matrices;
* A `Laplacian` constructor has been added;
* General callable types now have a convenience method for evaluating in Cartesian coordinates by calling them with complex numbers;
* The wavefront fitting error is now accessible;
* The analytical derivatives are now computed more efficiently;
* Derivatives now support `LaTeX` string formatting through `format_strings`;
* `format_strings` can now be called with a `WavefrontError` yielding its short `LaTeXString` representation;
* A new `Zernike.print_strings` function has been added which prints a batch of polynomial symbolic string representations in Unicode;
* An `mnv` function has been added which constructs a `Matrix` of Zernike indices and the associated input vector;
* Phase types can now be plotted in a new window by calling them with a `Screen` argument;
* `standardize` now has a method under the assumption that the weights are ordered;
* Calling `complex` on polynomial and gradient types will return a complex closure;
* `get_m` now has a single argument method for input `j`; `get_n` is now also exported;
* `getindex` & `setindex!` for `Wavefront`s have been extended for `m`, `n`, indices; single index `setindex!` has been fixed.

**Full Changelog**: https://github.com/Sagnac/Zernike.jl/compare/v5.6.0...v6.0.0

----

# v5.6.0

* Resetting `plotconfig` settings is now handled using functions `resize!` & `reset!` instead of setting properties;
* A Cartesian `W(x, y, OPD)` method has been added;
* The discrete `zplot` method now accepts more general 1-dimensional arrays for the phase domain values;
* Calling phase objects with no arguments will now plot them using `zplot`;
* A `Z(String)` shortcut method has been added which extracts the Unicode string for standard polynomials and their derivatives.

**Full Changelog**: https://github.com/Sagnac/Zernike.jl/compare/v5.5.1...v5.6.0

----

# v5.5.1

* The docstring for `radial_coefficients` has been fixed;
* The documentation for the phase matrix method has been clarified;
* The `stable` docs are now built from master.

**Full Changelog**: https://github.com/Sagnac/Zernike.jl/compare/v5.5.0...v5.5.1

----

# v5.5.0

* The labels for the Cartesian coordinates in the x-y plane on plots have been simplified to just `x` and `y`;
* The zero-based indexing methods for the coefficients of polynomial and wavefront error types have been extended; this includes the allowance of vector indices & defined boundary indices for both as well as `setindex!` for wavefront errors;
* The radial coefficients can now be directly computed as `Int`s or in arbitrary precision;
* `plotconfig` is no longer exported and must either be qualified with the module name as in `Zernike.plotconfig` or explicitly imported;
* Various bugs due to typos have been fixed;
* Documenter is now being used to generate documentation.

**Full Changelog**: https://github.com/Sagnac/Zernike.jl/compare/v5.4.0...v5.5.0

----

# v5.4.0

* An interface for working with partial derivatives and gradients has been added.

**Full Changelog**: https://github.com/Sagnac/Zernike.jl/compare/v5.3.1...v5.4.0

----

# v5.3.1

* The aperture transformation function now accepts non-standard length coefficient vectors.

**Full Changelog**: https://github.com/Sagnac/Zernike.jl/compare/v5.3.0...v5.3.1

----

# v5.3.0

* `standardize` has been overhauled; the mutating method `standardize!` no longer exists;
* The index types `Standard`, `Noll`, & `Fringe` have been added for ease of coefficient vector conversion.

**Full Changelog**: https://github.com/Sagnac/Zernike.jl/compare/v5.2.0...v5.3.0

----

# v5.2.2

* Added `get_j(n_max)`.

**Full Changelog**: https://github.com/Sagnac/Zernike.jl/compare/v5.2.1...v5.2.2

----

# v5.2.1

* `WavefrontError(Z::Polynomial)` is now defined and is the preferred method for directly converting Zernike polynomials;
* The property interface for `W.fit_to` now returns the proper value if the wavefront error was the result of a conversion;
* `plotconfig.reset` and `plotconfig.resize` now do nothing instead of error if set to `false`.

**Full Changelog**: https://github.com/Sagnac/Zernike.jl/compare/v5.2.0...v5.2.1

----

# v5.2.0

* Automatic resizing has been adjusted for `zplot` so that now the axis instead of the layout is rescaled;
* Single argument methods which radially evaluate with angle zero have been added to polynomial functions;
* `getindex` and `iterate` are now disallowed on `WavefrontOutput`;
* Convenience methods for indexing `Polynomial` and `WavefrontError` coefficients have been added;
* `j_to_noll` and `j_to_fringe` have been added;
* A few methods for working with reduced `precision` have been added: `reduce_wave`, `sieve`, & `standardize(W::WavefrontError)`;
* A new method `format_strings(m, n)` has been added for direct extraction of the Unicode and Latex strings.

**Full Changelog**: https://github.com/Sagnac/Zernike.jl/compare/v5.1.0...v5.2.0

----

# v5.1.0

* A `public` API has been declared and documented for Julia versions >= 1.11.0;
* Separate methods for the three public coefficient functions are now defined; namely: `radial_coefficients`, `wavefront_coefficents`, & `transform_coefficients`;
* `get_j(m, n)` & `get_mn(j)` are now exported;
* `Superposition` and `Product` public types are no longer exported;
* A couple of `reconstruct` methods have been added to bring functionality in line with `wavefront` and its expected input types;
* Automatic display has been fixed for `zplot` by now returning a `FigureAxisPlot`;
* Plots are now properly resized when their polynomial `Observable` changes;
* `DomainError`'s are now thrown instead of generic ones where appropriate;
* Fixed a corner case where `recap` was an empty vector for null wavefront errors;
* `show` has been fixed to produce concise output for `MixedPhase` containers;
* The `Sinusoid` internal type has been renamed to `Harmonic`.

**Full Changelog**: https://github.com/Sagnac/Zernike.jl/compare/v5.0.0...v5.1.0

----

# v5.0.0

* Performance has been improved for orders less than 21 by implementing the canonical formula for generating the radial polynomial coefficients;
* The `Model`-based dispatch system has been replaced with separate `Z`, `W`, `P` functions;
* Arithmetic between `Polynomial`s and `WavefrontError`s has been defined;
* The `Polynomial` and `WavefrontError` functions are now stored within the output types;
* `WavefrontError`s now store their `precision`;
* `show` now works properly when the object is part of a collection.

**Full Changelog**: https://github.com/Sagnac/Zernike.jl/compare/v4.2.0...v5.0.0

----

# v4.2.0

* Added a new `Zernike.metrics(ΔW::WavefrontError)` method;
* Improved the right click resizing interaction on the plot.

**Full Changelog**: https://github.com/Sagnac/Zernike.jl/compare/v4.1.0...v4.2.0

----

# v4.1.0

* `Zernike.Output` and `Zernike.WavefrontOutput` now properly utilize `show`;
* The pupil label is now placed above the toggle;
* `getindex` and `iterate` are no longer overloaded for `Zernike.Output` as the number of fields have been expanded.

**Full Changelog**: https://github.com/Sagnac/Zernike.jl/compare/v4.0.0...v4.1.0

----

# v4.0.0

* Several functions and variables have been renamed; notably the three main exported functions `Z`, `W`, and `P` are now `zernike`, `wavefront`, and `transform`;
* A few `WavefrontError` constructor methods have been added for use in creating them directly in terms of Zernike polynomials for given indices and coefficients;
* A new `zplot` function is now open to be used independently and accepts `Polynomial` and `WavefrontError` function types as well as quantized wavefront errors;
* The axis and plot is now returned by `zplot` and included in `Output` and `WavefrontOutput`;
* Deprecated functionality in ZernikePlot has been updated in accordance with the new version of GLMakie (`v0.9.0`);
* The ability to configure plot options through `plotconfig` has been added;
* The plot layout has been improved;
* A resizing mouse interaction which trims extra space on the plot after resizing the window is now bound to right click.

**Full Changelog**: https://github.com/Sagnac/Zernike.jl/compare/v3.2.0...v4.0.0

----

# v3.2.0

* A new wavefront error method has been added which accepts input coordinate vectors and a corresponding phase map;
* The more efficient method of `evalpoly` is now used to evaluate the radial polynomials;
* The shading / lighting is now turned off for the plots in order to improve the display.

**Full Changelog**: https://github.com/Sagnac/Zernike.jl/compare/v3.1.0...v3.2.0

----

# v3.1.0

* Added unit tests;
* `map_phase` won't mutate the input coordinate vectors anymore;
* The full standardized expansion coefficient vector is now stored within the `WavefrontError` structure;
* `ZPlot`: `getproperty` is now used instead of `getfield` when setting the pupil.

**Full Changelog**: https://github.com/Sagnac/Zernike.jl/compare/v3.0.0...v3.1.0

----

# v3.0.0

* Fitting to specific polynomials instead of the full range is now allowed using an input vector of `(m, n)` tuples containing the subset;
* Implemented a wavefront error matrix input method;
* The Model methods are now dispatched using parametric type selectors;
* The `scale` keyword option has been renamed to `finesse`;
* The new stack function is now used to efficiently construct the phase error fitting matrix;
* Slightly optimized the radial polynomial coefficient algorithm;
* Fixed conversion and indexing errors for wavefront error fitting under `n_max = 0` (mean piston output);
* Log transforms are now constructed more efficiently using broadcasted assignment;
* Adjusted the allowed input types for the wavefront error modelling;
* The window and plot title keyword arguments are now optional;
* Added a phase mapping function which converts a uniformly sampled set of data to coordinate vectors and a phase matrix;
* Computation of the RMS wavefront error has been improved by using the inner product;
* Validation of the pupil scaling factor is now carried out independent of translation.

**Full Changelog**: https://github.com/Sagnac/Zernike.jl/compare/v2.0.0...v3.0.0

----

# v2.0.0

* Pupil transform algorithm based on [Lundström & Unsbo (2007)](https://doi.org/10.1364/JOSAA.24.000569);
* The previous aperture scaling algorithm has been superseded by the pupil transform algorithm which is not only faster and more accurate for higher orders, but also includes translation and rotation for circular and elliptical pupils;
*  Functions to convert from Noll and Fringe single-index ordering to the standard ANSI order;
* Log transforms of the Zernike polynomials are now plotted for high orders;
* Various improvements and optimizations.

----

# v1.0.0

* Generation of Zernike polynomials;
* Fitting of wavefront errors;
* Aperture scaling algorithm based on [Janssen & Dirksen (2007)](https://doi.org/10.2971/jeos.2007.07012);
* Radial polynomial coefficient algorithm based on [Honarvar & Paramesran (2013)](https://doi.org/10.1364/OL.38.002487).

----
