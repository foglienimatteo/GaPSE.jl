## VERSION 0.10.0

Minor bump and not a patch: this release changes results that callers can observe.

### The library

- IMPORTANT CHANGE: replaced [Dierckx](https://github.com/kbarbary/Dierckx.jl) with our own cubic spline `MySpline` (new `src/Spline.jl`) for every 1D interpolation - `Cosmology`, `InputPS`, `IntegralIPS`, `IPSTools`, `BackgroundData`, `XiMatter`, ... It supports the `"Natural"`, `"Parabolic"` and `"ThirdDerivative"` initial conditions and provides its own `derivative`. `Dierckx` is now only a *test* dependency, where the suite uses it as an independent cross-check. Note that `MySpline` only supports `bc="error"`, so `spline_com_H` throws outside its range instead of clamping, and that it is not bit-identical to `Spline1D`: 33 reference files in `test/datatest` were regenerated and a few tolerances relaxed;

- IMPORTANT CHANGE: `Δχ = 0` is no longer an error. It is the exactly-collinear, coincident-point configuration, which the quadrature reaches deterministically, so the 21 `throw(AssertionError(...))` are now a fall-through to `zero(Δχ_square)`. The three `√(Δχ_square) > 1e-8 ? ... : 1e-8` clamps evaluated `√` before the comparison, so a negative argument raised a `DomainError` before the guard could act;

- every χ-integrated integrand now evaluates its analytic `Δχ → 0` limit instead of the `J * I_l^n` sum when `Δχ` is small, and `Δχ_min::AbstractFloat=1e-1` was added to the 26 integrands that lacked it. The derivations of the eight families of limits are in the new "The Δχ → 0 limits" pages of the manual, each obtained by expanding along `χ2 = χ1 + p Δχ` and checked to be independent of `p`;

- BUG FIX: those branches use a threshold **relative** to the local comoving distances where those are small, `Δχ ≥ min(Δχ_min, Δχ_min * max(χ1, χ2))`, applied to all thirty of them. The limits assume `y → 1`, which `Δχ → 0` forces only at *fixed, non-zero* distances; in the small-χ corner, where `χ1` and `χ2` vanish together at any `y`, the absolute threshold fired with `y` nowhere near 1 and returned a value a factor `1/y` too large. The `min` keeps the absolute cap, since the expansion needs `Δχ << 1/k_max`: a purely relative threshold would let the branch fire up to `Δχ = 0.1 χ ≃ 100` (measured: errors of 1000-3000%). Measured on `ξ_GNCxLD_Lensing_Lensing`, this removes a uniform 1.5% bias of the windowed multipoles at `s = 1000`;

- NUMERICAL STABILITY: in the three Lensing-Lensing integrands the brackets of `J_00`, `J_02`, `J_22` and `Δχ^2` are now written as expansions around the singular configuration, in `u = χ1²+χ2²`, `v = χ1χ2`, `t = y-1` and `w = (χ1-χ2)²`. All four vanish there while being evaluated as sums of terms of size `~χ⁴`: for `J_22` the bracket is exactly `8Δχ⁴` at `y = 1`, i.e. `8e-4` out of terms of `8e12` at `χ = 1000, Δχ = 0.1`, a ratio below `eps(Float64)`. Against exact rational arithmetic the old form had **0** correct digits there and the new one has 10 to 12. The rewrites are algebraic identities, checked to agree exactly on 3000 random rational `(χ1, χ2, y)`;

- corrected the `Δχ → 0` branch of `integrand_ξ_LD_Lensing`, whose `σ_0` coefficient was a factor 3 too large: its `J` coefficients are algebraically identical to those of `integrand_ξ_GNC_Lensing`, so the limit must be the same;

- BUG FIX: `Δχ_min` had been added to the *scalar* method of six `LD` integrands instead of the `Vector` one, so it never reached the limit branch;

- `WindowF` and `WindowFIntegrated` build their `GridInterpolations.RectangleGrid` once, in the constructor, instead of rebuilding it on every integrand evaluation;

- the `::Float64` annotations are now `::AbstractFloat` and the `Vector{Float64}` ones `Vector{T} where {T<:AbstractFloat}`, so another float type can flow through the code; nothing changes for `Float64`;

- removed the inert `x::T` type assertions from the `DEFAULT_*_OPTS` dictionaries: in expression position they always passed and constrained nothing, the types being enforced by `check_compatible_dicts`;

- the `print_map_*` functions truncate to 20 characters the keyword values they write into the output headers, so a large object no longer makes the header unreadable.

### Release and infrastructure

- the version is `0.10.0`, and `Project.toml` gains a `[compat]` section - it had none - with a lower bound per dependency and `julia = "1.12"`, plus the `[workspace]` table declaring `test`;

- the unit tests run on **Julia 1.12 only**: the 1.9 job and the advisory `continue-on-error` 1.12 probe are replaced by one blocking matrix entry per platform, `ubuntu-latest`/x64 (where coverage is taken) and `macos-latest`/aarch64. The previous workflow declared `aarch64` in the matrix but hard-coded `arch: x64` in the setup step, so macOS was in fact running x86_64 under Rosetta; native aarch64 is also 25% faster. The documentation workflow builds on 1.12 too, and the Julia badge follows;

- the `Dockerfile` moves to `quay.io/jupyter/julia-notebook:julia-1.12.7`. The Jupyter Docker Stacks publish on Quay since 2023, so the `jupyter/*` repositories on Docker Hub are the stale ones, not the project. The base image already provides Julia, IJulia and a registered kernel, so only GaPSE and the plotting extras are added, with `PYTHON` pinned to the stack's interpreter so PyCall does not build a private one;

- MAINTENANCE: the test-only and documentation-only dependencies are out of the package. `Project.toml` loses `ArbNumerics`, `IJulia`, `Documenter`, `Dierckx`, `NPZ`, `Suppressor` and `Test`; the last four, plus `DelimitedFiles` and `QuadGK`, live in the new `test/Project.toml`. `src/GaPSE.jl` drops `using Dierckx`, `using Test` and `using Documenter`, which also removes a latent name clash on `derivative`. `install_gapse.jl` was trimmed to the same list;

- FLAKY TEST FIX: `test_Spline.jl` drew its evaluation points with an unseeded `rand()` and compared them with a purely relative tolerance, which is meaningless where the derivative crosses zero: `linear range - nu=2` failed for about 6% of the seeds. The draws are now seeded and the six testsets that compare against `Dierckx` use `atol = RTOL * maximum(abs, ...)`. Failure rate 0% over 500 seeds, with 25x of headroom;

- `test/runtests.jl` gained the `TEST_BASICS`, `TEST_PP_PNG`, `TEST_LD`, `TEST_GNC`, `TEST_GNCxLD_LDxGNC` and `TEST_TWOSPECIES` switches, to run a subset of the suite while developing. They must all be `true` on the shared branches;

- replaced with four spaces the indentation tabs left in eight `src/GNCxLD_CrossCorrelations/*.jl` files and in `src/OtherUtils.jl`. The only tabs left in `src/` and `test/` are inside comments describing tab-separated input files.

### The theory in the manual

- CORRECTION, the `σ_i` are all finite and `I_0^0` does have an `s → 0` limit. An earlier version of this branch wrote their definitions with `∫_0^∞` and concluded that `σ_0` diverges in the ultraviolet, `σ_4` in the infrared, and that `I_0^0(s) ~ s^(-0.359)` has no limit. The step that produced it - replacing `j_l(qs)` by `(qs)^l` under the integral because `s → 0` - is not valid on an infinite range, since `q` runs to infinity and there is always a region `q > 1/s` where it is false. What the code integrates is the finite `[k_min, k_max]` that `IPSTools` hands to `xicalc`, where `s << 1/k_max` licenses the replacement uniformly. Measured over `[1e-5, 1e3]`, `I_0^0(s)/σ_0` is 0.166, 0.973 and 1.000000 at `s` = 1e-1, 1e-3 and 1e-6: a monotone approach to a finite limit. What remains true is that `σ_0` and `σ_4` *depend on their range* and cannot be quoted without it, while `σ_1`, `σ_2` and `σ_3` converge at both ends;

- `theory_IlnIntegrals.md` gains "Why the cut cannot be dropped", which splits `∫_0^∞` into `[0,k_min] + [k_min,k_max] + [k_max,∞)` and shows that the outer pieces converge but scale as `s^(α-3)`, a different power of `s` from the `s^(l-n)` of the middle one - which is where the spurious divergence came from;

- the section about the two `k` ranges is rewritten as a statement of what each range is for: the `I_l^n` need the wide `xicalc` grid `[1e-5, 1e3]`, while `σ_0` stands in for `I_0^0(Δχ_min)` in the limit branches and is therefore cut at `k_max ≈ 1/Δχ_min = 10`. On `WideA_ZA_pk.dat`, `σ_0(<10) = 18.58` against `I_0^0(0.1) = 23.73`, i.e. 22% below the number it replaces, where the `xicalc` range would give `143.3`, a factor 6 above it;

- the "from below" analysis of `Δχ_min` was measured before the Lensing-Lensing brackets were rewritten and is now re-measured: `J_22 I_2²` used to round to zero at `Δχ = 5e-2` and reach `-9.1e13` at `1e-3`, with the sum `6.6e7` times the limit; it is now `-4.9e5` and `-2.7e6` there, smooth and monotone. What limits `Δχ_min` from below is no longer the conditioning of the brackets but `fit_min = 0.05`, past which the `IntegralIPS` are diverging extrapolations. `Δχ_min = 1e-1` stays close to optimal, now because it is the largest value still of order `1/k_max` while staying above `fit_min`;

- new page "The input Power Spectrum" (`theory_InputPowerSpectrum.md`), which collects the metric and the conventions, the derivation of `P_m` from the primordial curvature spectrum, the asymptotic behaviour of `P(k)` at both ends, the extrapolation `InputPS` uses and which `σ_i` survive the removal of the cuts. Its Poisson equation carries the scale factor `a²(z)`, not the growth factor `D(z)`: the two coincide only in matter domination and differ by `D(0) ≃ 0.79` today, so the substitution would be a factor 1.6 in `P_m`. Every equation of the page is numbered `(P.n)`;

- NOTE on `n_P`: it is the small-`k` slope of the late-time MATTER Power Spectrum, **not** the `n_s - 1` of the dimensionless primordial curvature spectrum. The four powers of `k` of the Poisson equation minus the three of the `Δ²_R ↔ P_R` conversion are exactly what separates them;

- the `s → +∞` behaviour of the `I_l^n` and its Mellin-transform proof live in `theory_Ilnintegrals-mellin.md`: `I_l^n(s) → A/(2π²) M_l(3+n_P-n) s^-(3+n_P)`, an exponent independent of both `l` and `n`. The obvious route - substitute `x = qs` and replace `P(x/s)` by its small-`q` power law - is not legitimate, because the integral runs to `x = ∞` where `x/s` is not small; the Mellin route never forms that object. `I~_0^4` is excluded, its `μ = n_P - 1 = -0.04` sitting on the pole of `M_0(z)` at `z = 0` that the `-1` of its numerator subtracts away. The page is deliberately absent from the navigation menu while the argument is reviewed;

- the `Δχ → 0` limits were rewritten with step-by-step derivations and then split into nine pages, one per family. Every definition carries a `(D.n)` number and is referenced by it; each family page opens with a "Recap: everything this page needs"; every intermediate manipulation is written out instead of asserted, including the five Taylor coefficients of `B_22` around `y = 1` and the vanishing orders of the Families 2 and 3 numerators, which were previously quoted without proof;

- added to the same pages "A trap: the two small parameters are NOT of the same order" - along `χ2 = χ1 + p Δχ` one has `χ2-χ1 = O(Δχ)` but `1-y = O(Δχ²/χ²)`, so a term quadratic in `t` is the same size as one linear in `w` - and "A pattern worth noticing", on the low-order Legendre polynomials that keep appearing as leading coefficients;

- corrected the vanishing orders quoted for two families: the `J_02` and `J_04` numerators of Lensing x Doppler vanish as `Δχ³`, not `Δχ`;

- added the `MySpline` documentation and the derivation of its algorithm (`Spline.md` and `theory_SplineTheory.md`), plus `test/test_Spline.jl`;

- DOCS FIX: several pages wrote their formulas with `$...$` and `$$...$$`, or with LaTeX that Documenter cannot render (`longtable`, `parbox`, `multirow`); they now use ```` ```math ```` blocks and plain `align*`. The Theory pages are renamed with a `theory_` prefix, and the hand-written tables of contents and section links - which used GitHub-style anchors that Documenter does not generate, and were therefore dead on the published manual - are `@contents` blocks and `@ref` links, which Documenter validates at build time;

- DOCSTRING FIX: the docstrings of `integrand_ξ_GNC_Newtonian_Lensing` and `integrand_ξ_GNCxLD_Newtonian_Lensing` had the wrong sign in the `b_1` part of `J^{δκ}_{02}`.

### The `theory/` directory

- added `theory/`, which collects the notebooks that reproduce the figures and the numbers of the Theory pages, together with the plots and data they produce: `Iln_terms.ipynb`, `sigma_i.ipynb`, `input_ps.ipynb` and `spherical_bessels.ipynb`. They are notebooks only - the parallel `.jl` scripts were removed, since keeping the two in sync was manual and they had already drifted;

- all four share one plotting API: `plot_kwargs(kwargs...)` merging overrides onto the shared defaults, `logticks(lo, hi; step)`, `hlabel` for the horizontal `yguidefontrotation = -90` labels, `vspec`/`vlines!` for the markers, and every label, scale and tick range a keyword. Plotting is split from saving, so a figure can be looked at and tweaked before being written out. Sections follow the same layout - machinery collapsed through `jp-MarkdownHeadingCollapsed`, results left open;

- they rebuild their environment themselves: they `Pkg.develop` the GaPSE of this repository when it is not reachable, add what `Project.toml` does not declare yet, then `Pkg.resolve()` and `Pkg.instantiate()`. `resolve` before `instantiate` matters, otherwise a `Project.toml` restored over an older manifest fails. `pyplot()` falls back to `gr()` where matplotlib is missing;

- `sigma_i.ipynb` studies the moments every `Δχ → 0` limit reduces to, plotting the five integrands and the fraction of each collected below a given `q`. `σ_2` and `σ_3` agree to better than 0.1% between any two reasonable cuts, while `σ_0` keeps growing with `k_max` and `σ_4` as `k_min` is lowered: for those two the cut *is* the value. Worth knowing: the input file is tabulated only on `[2.1e-7, 20.2]`, so most of a `σ_0` integrated to `k_max = 1e3` comes from the extrapolation, not from data;

- `spherical_bessels.ipynb` reproduces the figures of the "Spherical Bessel Functions" page and adds, under each of the two log-log ones, a ratio panel with `ylims = (0.9, 1.1)`: the small-`x` expansion is good to 10% up to `x = 0.79` for `ℓ = 0`, drifting only to `1.52` by `ℓ = 4`;

- `Iln_terms.ipynb` adds `I_direct`/`I04_tilde_direct` for the direct quadrature, marks in grey the regions where an `IntegralIPS` is an extrapolation, and carries the large-`s` asymptote on each figure. Plotted against the right `σ_i`, the `I_l^n` match their analytic asymptotes to four digits.

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
