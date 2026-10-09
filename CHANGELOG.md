## VERSION 0.10.0

Minor bump and not a patch: this release changes results that callers can observe.

### The code

- **IMPORTANT CHANGE**: replaced [Dierckx](https://github.com/kbarbary/Dierckx.jl) with our own cubic spline `MySpline` (new `src/Spline.jl`) for every 1D interpolation; `Dierckx` is now only a *test* dependency, where the suite uses it as an independent cross-check. Note that `MySpline` only supports `bc="error"`, so `spline_com_H` throws outside its range instead of clamping, and that it is not bit-identical to `Spline1D`: 33 reference files in `test/datatest` were regenerated and a few tolerances relaxed;

- **IMPORTANT CHANGE**: the `Δχ → 0` limits
  * `Δχ = 0` is no longer an error. It is the exactly-collinear, coincident-point configuration, which the quadrature reaches deterministically, so the 21 `throw(AssertionError(...))` are now a fall-through to `zero(Δχ_square)`. The three `√(Δχ_square) > 1e-8 ? ... : 1e-8` clamps evaluated `√` before the comparison, so a negative argument raised a `DomainError` before the guard could act;
  * every `χ`-integrated integrand now evaluates its analytic `Δχ → 0` limit instead of the `J * I_l^n` sum when `Δχ` is small, and `Δχ_min::AbstractFloat=1e-1` was added to the 26 integrands that lacked it. The derivations of the eight families of limits are in the new "The Δχ → 0 limits" pages of the manual, each obtained by expanding along `χ2 = χ1 + p Δχ` and checked to be independent of `p`;
  * BUG FIX: those branches use a threshold **relative** to the local comoving distances where those are small, `Δχ ≥ min(Δχ_min, Δχ_min * max(χ1, χ2))`, applied to all thirty of them. The limits assume `y → 1`, which `Δχ → 0` forces only at *fixed, non-zero* distances; in the small-χ corner, where `χ1` and `χ2` vanish together at any `y`, the absolute threshold fired with `y` nowhere near 1 and returned a value a factor `1/y` too large. The `min` keeps the absolute cap, since the expansion needs `Δχ << 1/k_max`: a purely relative threshold would let the branch fire up to `Δχ = 0.1 χ ≃ 100` (measured: errors of 1000-3000%). Measured on `ξ_GNCxLD_Lensing_Lensing`, this removes a uniform 1.5% bias of the windowed multipoles at `s = 1000`;
  * corrected the `Δχ → 0` branch of `integrand_ξ_LD_Lensing`, whose `σ_0` coefficient was a factor 3 too large: its `J` coefficients are algebraically identical to those of `integrand_ξ_GNC_Lensing`, so the limit must be the same;
  * BUG FIX: `Δχ_min` had been added to the *scalar* method of six `LD` integrands instead of the `Vector` one, so it never reached the limit branch;

- NUMERICAL STABILITY: in the three Lensing-Lensing integrands the brackets of `J_00`, `J_02`, `J_22` and `Δχ^2` are now written as expansions around the singular configuration, in `u = χ1²+χ2²`, `v = χ1χ2`, `t = y-1` and `w = (χ1-χ2)²`. All four vanish there while being evaluated as sums of terms of size `~χ⁴`: for `J_22` the bracket is exactly `8Δχ⁴` at `y = 1`, i.e. `8e-4` out of terms of `8e12` at `χ = 1000, Δχ = 0.1`, a ratio below `eps(Float64)`. Against exact rational arithmetic the old form had **0** correct digits there and the new one has 10 to 12. The rewrites are algebraic identities, checked to agree exactly on 3000 random rational `(χ1, χ2, y)`;

- `WindowF` and `WindowFIntegrated` build their `GridInterpolations.RectangleGrid` once, in the constructor, instead of rebuilding it on every integrand evaluation;

- the `::Float64` annotations are now `::AbstractFloat` and the `Vector{Float64}` ones `Vector{T} where {T<:AbstractFloat}`, so another float type can flow through the code; nothing changes for `Float64`;

- removed the inert `x::T` type assertions from the `DEFAULT_*_OPTS` dictionaries: in expression position they always passed and constrained nothing, the types being enforced by `check_compatible_dicts`;

- the `print_map_*` functions truncate to 20 characters the keyword values they write into the output headers, so a large object no longer makes the header unreadable.

- **IMPORTANT CHANGE**: the `k_min`/`k_max` of `IPSTools` now bound the `xicalc` that builds the nine `I_l^n` and the `quadgk` of `Ĩ^4_0` too, not only the five `σ_i` (those two were hard-coded to `1e-5, 1e3`), and `s0`, the first point of the `s` grid `xicalc` returns, became a keyword. They are the same integral - `I_0^0(s) → σ_0` for `s → 0` is what every `Δχ → 0` limit asserts - so they cannot live on two different supports. It is not a cosmetic change: with `k_max = 10` the `I_l^n` lose power below `s ≃ 0.2` (`I_0^0(0.06)` drops to 57%), which is the region the `Δχ → 0` switch feeds on, so the GNC Lensing-Lensing `L = 0` multipole of the test cosmology moves up to 38% at `s = 1000` while the GNC sum moves 0.5%. The reference files store the old mixed convention and have to be regenerated;

- the same work found, but did NOT fix, a defect in `power_law_from_data`: in the `con == true` branch the relative errors are compared *signed* against 0.05, so an arbitrarily bad negative one passes, and the two fall-back `curve_fit` calls are given a 3-element `p0` against 2-parameter models, so `si, b, a = vcat(vals_2, vals_1[3])` silently takes `a = p0[3] = 0.0`. With `fit_min = 0.05` and `k_max = 10` the left fit of `I_0^0` collapses onto a degenerate solution and the extrapolation below `fit_min` jumps by a factor 3e3. Nothing GaPSE currently computes reaches below `fit_min`, so no result moves today, but it is a landmine for a smaller `Δχ_min` or `s_min`. Fixing it moves results of its own (the signed test already fires on the current defaults), so it wants its own branch. Measurements in the new `theory/kmin_kmax.ipynb`;

- THREADS: the four `map_ξ_*_multipole` functions spread their `s` values over the available threads, through the new `map_over_ss` (`src/OtherUtils.jl`). The `s` points are independent and the `Cosmology` is read-only on that path, so the results are bit-for-bit the serial ones; the scheduling is `:dynamic`, the cost of one `s` growing with `s`. Measured 3.56x on 4 threads. With one thread it is an ordinary loop, so Julia must be started with `-t auto` (or `JULIA_NUM_THREADS`) to get anything out of it;

- new `print_log_generic`/`print_log` (`src/OtherUtils.jl`): log on `stdout`, on an open stream or on a file, with an optional `[yyyy-mm-dd HH:MM:SS]` stamp, and with a method that redirects `stdout` around a function that prints on its own. `Dates` is a new (stdlib) dependency.


### Release and infrastructure

- `Project.toml` has now a `[compat]` section with a lower bound per dependency and `julia = "1.12"`, plus the `[workspace]` table declaring `test`;

- the unit tests run on **Julia 1.12 only**: the 1.9 job and the advisory `continue-on-error` 1.12 probe are replaced by one blocking matrix entry per platform, `ubuntu-latest/x64` (where coverage is taken) and `macos-latest/aarch64`. The previous workflow declared `aarch64` in the matrix but hard-coded `arch: x64` in the setup step, so macOS was in fact running x86_64 under Rosetta;

- DOCKERFILE: the `Dockerfile` moves to `quay.io/jupyter/julia-notebook:julia-1.12.7`. The Jupyter Docker Stacks publish on Quay since 2023, so the `jupyter/*` repositories on Docker Hub are the stale ones, not the project. The base image already provides Julia, IJulia and a registered kernel, so only GaPSE and the plotting extras are added, with `PYTHON` pinned to the stack's interpreter so PyCall does not build a private one;

- DEPENDENCIES: the test-only and documentation-only dependencies are out of the package. `Project.toml` loses `ArbNumerics`, `IJulia`, `Documenter`, `Dierckx`, `NPZ`, `Suppressor` and `Test`; the last four, plus `DelimitedFiles` and `QuadGK`, live in the new `test/Project.toml`. `src/GaPSE.jl` drops `using Dierckx`, `using Test` and `using Documenter`, which also removes a latent name clash on `derivative`. `install_gapse.jl` was trimmed to the same list;

- TEST FIX: `test_Spline.jl` drew its evaluation points with an unseeded `rand()` and compared them with a purely relative tolerance, which is meaningless where the derivative crosses zero: `linear range - nu=2` failed for about 6% of the seeds. The draws are now seeded and the six testsets that compare against `Dierckx` use `atol = RTOL * maximum(abs, ...)`.

- `test/runtests.jl` gained the `TEST_BASICS`, `TEST_PP_PNG`, `TEST_LD`, `TEST_GNC`, `TEST_GNCxLD_LDxGNC` and `TEST_TWOSPECIES` switches, to run a subset of the suite while developing. They must all be `true` on the shared branches;

- a third matrix entry runs the whole suite on `ubuntu-latest/x64` with `JULIA_NUM_THREADS: 4`, against the same reference files: the parallel path has to give the single-threaded numbers. 4 is a literal and not `auto` on purpose, `auto` degrading silently to one thread on a smaller runner;

- `test/runtests.jl` prints at the start the thread-related environment variables, `Threads.nthreads()`, `Sys.CPU_THREADS`, the BLAS thread count and `Sys.cpu_summary()`, so a log says on how many threads it ran;

- NOTEBOOK FIX: the notebooks in `ipynbs/` ran `include(PATH_TO_GAPSE * "src/GaPSE.jl")`, which evaluates the sources in `Main`: GaPSE's own `Project.toml` is never read, so its dependencies are looked for in the kernel's active project and the include fails with `Package TwoFAST [...] is required but does not seem to be installed`. They now just `using GaPSE`, out of the new `ipynbs/Project.toml`, which declares GaPSE through a `[sources]` entry - the only place where the path of an unregistered package is written down *and* tracked by git, since `Pkg.develop` records it in the gitignored `Manifest.toml`. The new `ipynbs/README.md` documents how that environment was built and why. The same `[sources]` entry was added to `theory/Project.toml`;

- NOTEBOOK FIX: the `using` lines of the `ipynbs/` notebooks listed six packages that no cell ever calls (`ProgressMeter`, `QuadGK`, `Trapz`, `LegendrePolynomials`, `SpecialFunctions`, `TwoFAST`), left over from the `include` days and all of them dependencies of GaPSE anyway; they are out of both the notebooks and `ipynbs/Project.toml`;

- BACKEND: every notebook selects `gr()`, and `PyPlot` leaves both `Project.toml` without a replacement. Under Julia 1.12 `pythonplot()` warns that it accesses `Plots._py_drawfig` "in a world prior to its definition world" and "will error in future versions of Julia" - `Plots` `include`s its backend file lazily, so its bindings land in a later world age than the code calling them - while `pyplot()` is implemented in `Plots/src/backends/deprecated/`. GR needs no Python at all and draws everything these notebooks use, `st = :surface` and `heatmap` included;


### The theory in the manual

- THEORY DOCS: A lot of documentation about the theory (physical and numerical) of GaPSE has been written on the docs:
  * the Theory pages are all named with a `theory_` prefix inside the `docs/src` dir
  * added the `MySpline` documentation and the derivation of its algorithm (`Spline.md` and `theory_SplineTheory.md`), plus `test/test_Spline.jl`;
  * added derivation of all `Δχ → 0` limits
  * added theory about Power Spectrum

- DOCS FIX: several pages wrote their formulas with `$...$` and `$$...$$`, or with LaTeX that Documenter cannot render (`longtable`, `parbox`, `multirow`); they now use ```` ```math ```` blocks and plain `\begin{align*}`. Also, the hand-written tables of contents and section links, which used GitHub-style anchors that Documenter does not generate and were therefore dead on the published manual, are `@contents` blocks and `@ref` links, which Documenter validates at build time;

- DOCSTRING FIX: the docstrings of `integrand_ξ_GNC_Newtonian_Lensing` and `integrand_ξ_GNCxLD_Newtonian_Lensing` had the wrong sign in the `b_1` part of `J^{δκ}_{02}`.


- added the `theory/` directory, which collects the notebooks that reproduce the figures and the numbers of the Theory pages, together with the plots and data they produce
  * they share `theory/Project.toml`, which declares GaPSE through a `[sources]` entry, so a fresh clone only needs `Pkg.instantiate()`, and the default IJulia kernel activates it by itself (it runs Julia with `--project=@.`);
  * `sigma_i.ipynb` studies the moments every `Δχ → 0` limit reduces to, plotting the five integrands and the fraction of each collected below a given `q`
  * `spherical_bessels.ipynb` reproduces the figures of the "Spherical Bessel Functions" page
  * `Iln_terms.ipynb` studies the Iln integrals
  * `spline_comparison.ipynb` compares `MySpline` with `Dierckx`
  * `deltachi_limits.ipynb` looks at the `Δχ_min` switch between the `J ⋅ I_l^n` sum and the analytic limit, on two GNC auto-correlations
  * `kmin_kmax.ipynb` scans `k_min` and `k_max` over the `σ_i`, the `I_l^n` and their power-law fits, and shows that it is the fit, not the range, that moves a TPCF



## branch oneapi -> should have lead to version 0.9.0

- added `MySpline`

- trying to parallelize the code with `KernelAbstractions`; seems that the GPU offloading is overkill, due to small size of matrixes in the single (Lensing-... and IntegratedGP-...) and double integral terms (Lensing-Lensing, IntegratedGP-IntegratedGP)




## VERSION 0.8.0

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
