## HEAD

- IMPORTANT CHANGE: replaced [Dierckx](https://github.com/kbarbary/Dierckx.jl) with our own cubic spline `MySpline`, implemented in the new `src/Spline.jl`, for all the 1D interpolations of the library (`Cosmology`, `InputPS`, `IntegralIPS`, `IPSTools`, `BackgroundData`, `XiMatter`, ...). It supports the `"Natural"`, `"Parabolic"` and `"ThirdDerivative"` initial conditions and provides its own `derivative`. `Dierckx` is no longer used anywhere inside `src/`; it is kept as a dependency only because the test suite uses it as an independent cross-check;

- added the `MySpline` documentation and the derivation of the algorithm to the manual (`docs/src/Spline.md` and `docs/src/SplineTheory.md`), plus `test/test_Spline.jl`;

- NOTE: `MySpline` and `Dierckx.Spline1D` are not bit-identical, so 33 reference files in `test/datatest` were regenerated and a few test tolerances relaxed accordingly. Note also that `MySpline` only supports `bc="error"`, so the `spline_com_H` of `BackgroundData` throws outside its range instead of clamping to the nearest value;

- the `::Float64` annotations of function arguments and keyword arguments are now `::AbstractFloat`, and the `Vector{Float64}` ones are `Vector{T} where {T<:AbstractFloat}`, so that a float type other than `Float64` can flow through the code. This is a pure widening: nothing changes for `Float64` input;

- removed the inert `x::T` type assertions from the `DEFAULT_IPS_OPTS`, `DEFAULT_IPSTOOLS_OPTS`, `DEFAULT_WFI_OPTS`, `DEFAULT_FMAP_OPTS_hcub` and `DEFAULT_FMAP_OPTS_trap` dictionaries. In expression position `x::T` is a runtime assertion that always passed and constrained nothing; the types accepted for these options are enforced by `check_compatible_dicts`, which is unaffected;

- `WindowF` and `WindowFIntegrated` now build their `GridInterpolations.RectangleGrid` once, in the constructor, and store it. `spline_F` and `spline_integrF` used to rebuild it on every call, i.e. once per integrand evaluation;

- the `print_map_*` functions now truncate to 20 characters the keyword-argument values they write into the header of the output files, so that passing a large object no longer makes the header unreadable;

- the unit tests now also run on pull requests, Julia 1.12 was added as a second non-blocking job, Windows was dropped and macOS moved to `aarch64` (an `x86_64` Julia running under Rosetta 2 on the Apple-silicon runners produced spurious `DomainError`s from `cos`);

- `test/runtests.jl` gained the `TEST_BASICS`, `TEST_PP_PNG`, `TEST_LD`, `TEST_GNC`, `TEST_GNCxLD_LDxGNC` and `TEST_TWOSPECIES` switches, to run only a subset of the suite while developing. They must all be `true` on the shared branches;

- IMPORTANT CHANGE: `Δχ = 0` is no longer an error. It is the exactly-collinear, coincident-point configuration, which the quadrature reaches deterministically (`μ = ±1` are nodes of both `:lobatto` and `:trap`, and `suit_sampling` places a dense sub-grid right on `χ = s`), so the 21 `throw(AssertionError(...))` are now a fall-through to `zero(Δχ_square)`, with the `throw` kept as a comment next to each site. The three `√(Δχ_square) > 1e-8 ? √(Δχ_square) : 1e-8` clamps were normalised with the rest: `√` was evaluated before the comparison, so a negative argument raised a `DomainError` before the guard could act;

- every χ-integrated integrand now evaluates its analytic `Δχ → 0` limit instead of the `J * I_l^n` sum when `Δχ < Δχ_min`, and `Δχ_min::AbstractFloat=1e-1` was added to the 26 integrands that did not have it. The derivations of the eight families of limits are in the new "The Δχ → 0 limits" page of the manual (`docs/src/DeltaChiLimits.md`), each one obtained by expanding along `χ2 = χ1 + p Δχ` and checked to be independent of `p`;

- corrected the `Δχ → 0` branch of `integrand_ξ_LD_Lensing`, which read `9/4 * 4/15 * (5σ_2 + 6σ_0 χ2^2)` = `3σ_2 + 18/5 χ2^2 σ_0`: its `J` coefficients are algebraically identical to those of `integrand_ξ_GNC_Lensing`, so the limit must be the same, and the `σ_0` coefficient was a factor 3 too large;

- KNOWN BUG, documented but not yet fixed: the `Δχ → 0` limits assume `y → 1`, which is what `Δχ → 0` forces at *fixed, non-zero* comoving distances. It does not hold in the small-χ corner, where `χ1` and `χ2` go to zero together at any `y`, and every double-χ grid starts at `χ = 1e-6 * s ≃ 4e-4`, so the whole first row and column of the grid take the branch with full weight. There the correct limits are `3 y σ_2 + 6/5 χ1^2 σ_0` (family 1) and `A/3 y σ_2` (family 8) instead of the same expressions with `y = 1`, so the seven double-χ integrands return a value a factor `1/y` too large, and of the wrong sign for `y < 0`. This is what moves `ξ_GNCxLD_Lensing_Lensing` by 2-4% and makes 8 assertions of `test_GNCxLD_SumXiMultipoles_P1.jl` fail. The derivation, the measured effect and the proposed relative guard `Δχ < Δχ_min * max(χ1, χ2)` are in the "A second way to reach Δχ = 0" section of the manual;

- NOTE on the unit tests: `test_GNCxLD_SumXiMultipoles_P1.jl` reports 28 failing assertions out of 672. 20 of them are older than this branch: they appear with the switch to `MySpline` (they are absent at `Abstractfloat (#10)` and present right after the reference data were regenerated), so the `datatest/GNCxLD_SumXiMultipoles` files still do not match what `MySpline` produces, at `rtol = 1.2e-2`. The remaining 8 are the small-χ corner bug above;

- DOCS FIX: `docs/src/DeltaChiLimits.md` and `docs/src/SplineTheory.md` wrote their formulas with `$...$` and `$$...$$`. Documenter uses the Julia Markdown flavour, which wants ` ``x`` ` and ` ```math ` blocks, so it emitted "Unexpected Julia interpolation in the Markdown" for every formula and rendered the pages with the raw LaTeX source visible instead of the equations. Both pages now use the Documenter syntax, and the build is free of those warnings;

- BUG FIX: `Δχ_min` had been added to the *scalar* method of the six `LD` integrands `integrand_ξ_LD_IntegratedGP`, `..._Doppler_IntegratedGP`, `..._Lensing_Doppler`, `..._Lensing_IntegratedGP`, `..._Lensing_LocalGP` and `..._LocalGP_IntegratedGP`, while the body that uses it lives in the `Point` method: every call raised `UndefVarError: Δχ_min not defined`, aborting the whole `LD` half of the test suite. The keyword now sits on the `Point` method, as in the `GNC` and `GNCxLD` families, and the scalar methods forward it through `kwargs...`;

- added the new `theory/` directory, meant to collect the Julia scripts (`.jl`) and notebooks (`.ipynb`) - together with the plots and the data they produce - that reproduce the figures of the manual and help investigating the theoretical behaviour of GaPSE. Its first content is `theory/Iln_terms.jl`/`theory/Iln_terms.ipynb`, which plot in log-log scale all the `I_l^n` stored in `IPSTools` (and the regularized `I~_0^4`), each one against its small-`s` asymptote, plus a figure with all of them together. The directory has its own `Project.toml`, so that `Plots` and `PyPlot` are not added to the GaPSE dependencies;

- added the "The I_l^n integrals" page of the manual (`docs/src/IlnIntegrals.md`), which defines the `I_l^n`, proves that `I_l^n(s) -> sigma_{n-l} s^{l-n} / (2l+1)!!` and `I~_0^4(s) -> -sigma_2 / (6 s^2)` for `s -> 0`, and shows the figures produced by `theory/Iln_terms.jl`;

- BUG FIX in `theory/Iln_terms.jl`: the plotted `I_l^n` did not match their analytic asymptotes, for three compounding reasons, none of them a flaw of the derivation. (i) An `IntegralIPS` is a spline only between its `left` (`= fit_min = 0.05`) and `right` fields: below `left` it returns a power law fitted on `[0.05, 0.5]` and seeded with a negative exponent, so the first three decades of the old plots were an extrapolation and not the integral. (ii) `IPSTools` hard-codes `kmin, kmax = 1e-5, 1e3` for the `xicalc` call that builds the `I_l^n`, while its `k_min`/`k_max` keywords only affect the `sigma_i` it stores; the script was computing the `sigma_i` of the asymptotes over `[1e-6, 10]`, which is harmless for `sigma_2` (0.1%) but wrong by a factor `5e4` for `sigma_-2` and `5e8` for `sigma_-4`. That is exactly why only `I_0^2`, `I_1^3` and `I~_0^4`, whose limits depend on `sigma_2` alone, appeared to be correct. (iii) The expansion needs `s << 1/k_max = 1e-3`, i.e. 50 times below `fit_min`, so the validity window of the limits and that of the spline do not overlap and the limits cannot be seen through `IPSTools` at all;

- `theory/Iln_terms.jl`/`.ipynb` now use the right `sigma_i`, mark in grey the regions where an `IntegralIPS` is an extrapolation, and add `I_direct`/`I04_tilde_direct`, a brute-force quadrature over the same extremes `xicalc` uses, which confirms the analytic limits to four digits for `s <= 1e-4`. A new figure `ratios.png` plots the ratio to the asymptote for all of them. `docs/src/IlnIntegrals.md` gained the corresponding "A warning before looking at the plots" section, with the measured tables;

- added `SpecialFunctions` to `theory/Project.toml`, needed by the direct quadrature;

- DOCS FIX: `docs/src/SphericalBesselFunctions.md` was written with LaTeX that Documenter cannot render (`longtable`, `parbox`, `multirow`, `makecell`, `equation`/`label`, and the personal macros `\deriv`, `\secderiv`, `\dd`, `\versor`). The two-column table is now a sequence of subsections, the macros are expanded and the page renders. Its `j_l` identities were all verified numerically, and one of them was wrong: in `int_0^inf j_l(Kx) j_l(kx) dx = pi/(2(2l+1)) K^l/k^(l+1)` it is the SMALLER argument that goes to the numerator, so the condition is `K < k` and not `K > k`. The validity range `-2l-1 < p < 1` was added to the `int x^p j_l^2` formula, and the regression for the first zero of `j_l` (`x = 4.75 + 1.05 l`) was confirmed (fit gives `4.7466 + 1.0513 l`);

- added `theory/spherical_bessels.jl` and `theory/spherical_bessels.ipynb`, which reproduce the figures of the "Spherical Bessel Functions" page and check numerically the two quantitative claims it makes. Besides the `j_l(x)` figure that page already showed (which had no script behind it), they produce a log-log comparison of each `j_l` with its leading small-`x` term `x^l/(2l+1)!!` - showing that the truncation is legitimate only up to `x ~ 1`, which is the same statement as the `s << 1/k_max` condition of the `I_l^n` - and the first zero of `j_l` for `0 <= l <= 100` against both the linear regression and the exact `l + 1.8557 l^(1/3)`. The zeros are found by bisection and the straight line by least squares, so no new dependency is needed;

- cosmetic fixes to the `theory/Iln_terms.jl` figures: the decade ticks are now generated from the plotted range (`logticks`) instead of once and for all, so they no longer pile up on the left edge of the figures that do not span all the 11 decades; the asymptote is drawn only up to `s = 1`, since being a pure power law it otherwise spans 25 decades and squashes everything else; the vertical range is set by the data alone; and the direct-quadrature curve is drawn on top of the asymptote rather than under it;


## development branch qls

- added `readchoosen` and `readxchoosey` functions in `src/OtherUtils.jl`;

- changed logo of GaPSE;

- added `ipynbs/eBOSS_Window.ipynb` and related files;

- added `Dockerfile` for building the image `gapse-julia-1.9.1:0.8.0a`;

- IMPORTANT CHANGE: NOW YOU CAN GO FURTHER THAN `z=1.5`!

- IMPORTANT CHANGE: TWO TRACERS IMPLEMENTATION!
  Now you have two comsological biases sets `b1`, `s_b1`, `𝑓_evo1` and `b2`, `s_b2`, `𝑓_evo2` in the `CosmoParams` struct (you cannot specify anymore them as `b`, `s_b`, `𝑓_evo`);

- added all the previously missing docstrings in `GNCxLD` and `LD` TPCFs;

- added `ipynb/Computations_b1p5-sb0-fevo0.ipynb`, `Computations_b1p5-sb0-fevo0.jl` and `Generic_Window.jl` for the analysis of the PNG, even with a generic window function

- added `src/PowerSpectraGenWin.jl` and `src/WindowF_QMultipoles.jl`: now it's possible to compute the PS for a generic window!


## VERSION 0.7.0

- huge improvements on docstrings and API.

- renamed `src/PowerSpectrum.jl` to `src/PowerSpectra.jl`

- creation of the package `GaPSE.jl`; it will still take a while

## VERSION 0.6.0

- Made the docstrings for almost all the functions in the code; only the `LD`, `GNCxLD` and `LDxGNC` TPCFs are huge holes in this picture

- Improved the structure of the functions that manage the Integrated Window Function (`print_map_F` and `WindowF` in particular); some options have changed name
  
- Improved the structure of the functions that manage the Integrated Window Function (`print_map_IntegratedF` and `WindowFIntegrated` in particular); some options have changed name

- Improved the code for `PS_multipole` and associated functions; some options have changed name

## VERSION 0.5.0

- Implemented tests.

- Added the `PNG.jl` source code file for the analysis of the Primordial Non-Gaussianites

- Calibrated the integrals on `χ` for some of the GNC terms

- Now you can compute the integral on `μ` (for `GNC`, `LD`,  `GNCxLD` and `LDxGNC` TPCFs) through three algorithms,and you can specify with one to use with the keyword argument `alg`:  `:quad`, `:lobatto`, and `:trap`.

## VERSION 0.4.0

- Added tests for all the previous features.

- Added `PPXiMatter.jl` for calculating the multipoles (with and without window function) of the matter TPCF.

- Added `XiMatter.jl` in order to compute a TPCF from a PS; the case of particular interests concerns the computation of the matter TPCF from the PS deriving from CLASS code.

- Added the file `WindowFIntegrated.jl`: now the computation of the TPCFs with the window function is correct.

- Added the functions `ξ_GNCxLD_multipole`, `map_ξ_GNCxLD_multipole`, `print_map_ξ_GNCxLD_multipole`, their `sum_` extensions and their `LDxGNC` counterparts for the computations of the relativistic GNC effects.

- Implemented all the `GNCxLD` and `LDxGNC` TPCFs.
  
- Added the functions `ξ_GNC_multipole`, `map_ξ_GNC_multipole`, `print_map_ξ_GNC_multipole` ant their `sum_` extensions for the computations of the relativistic GNC effects.

- Implemented all the `GNC` TPCFs.
  
- Changed the name of the following functions (because they refer to the LD perturbations):
  -  `ξ_multipole` in `ξ_LD_multipole`;
  -  `map_ξ_multipole` in `map_ξ_LD_multipole`;
  -  `print_map_ξ_multipole` in `print_map_ξ_LD_multipole`;
  -  `sum_ξ_multipole` in `sum_ξ_LD_multipole`;
  -  `map_sum_ξ_multipole` in `map_sum_ξ_LD_multipole`;
  -  `print_map_sum_ξ_multipole` in `print_map_sum_ξ_LD_multipole`;


## VERSION 0.3.0

- Added the functions for the evaluation of the Doppler auto-CF with the Plane-Parallel approximation and their tests

- Modified default range 10 .^ range(-1, 3, length= N_log) to  10 .^ range(0, 3, length= N_log) 

- Added tests for PS evaluation


## VERSION 0.2.0

- Modified the `CosmoParams` struct: now you should pass dictionaries for `InputPS` and `IPSTools` options

- Added the `ipynbs/TUTORIAL.ipynb` notebook, which is a small tour of how the code should be used

- Changed the keyword argument `s_1` into `s1` (in functions such as `map_ξ_multipole`, `map_sum_ξ_multipole`, ...)

- Removed the functions `integrand_on_mu`, `integral_on_mu`, `map_integral_on_mu` and `print_map_integral_on_mu`; added `integrand_ξ_multipole`. This restructuration let the code more flexible.

- Optimized `F_map` and re-written its keyword argument structure

- Improved the documentation of the code, and the `README.md`



## VERSION 0.1.0

Created the basic structure of the program, with the structs `WindowF`, `InputPS`, `CosmoParams` and `Cosmology`.
It's possible to:
- read an arbitrary input matter power spectrum, produced by CLASS
- read an arbitrary input background data file, produced by CLASS
- write and read a map of the window F function
- calculate all the General Relativistic (GR) effects Two-Point Correlation Functions (TPCFs) multipole functions (both auto-correlation and cross-correlation) arising from the perturbed luminosity distance [[1]](#1) for an arbitrary multipole order 
- calculate the sum of all these TPCFs

Some notebook are already provided, inside the directory `ipynb`, such as:
- `ALL_CF_L.ipynb`
- `ALL_Multipoles.ipynb`
- `PS_Multipoles.ipynb`
- `SUM_CF.ipynb`
  


## References

<a id="1">[1]</a> 
See the equation (C.23) with ``s`` Dalal, Doré et al., _Imprints of primordial non-Gaussianities on large-scale structure_ (2008), American Physical Society, DOI: 10.1103/PhysRevD.77.123514, 
url: https://journals.aps.org/prd/abstract/10.1103/PhysRevD.77.123514
