## HEAD

- MAINTENANCE: the test-only and documentation-only dependencies are out of the package. `Project.toml` loses `ArbNumerics` and `IJulia`, which are not used anywhere in the repository, `Documenter`, which belongs to `docs/Project.toml`, and `Dierckx`, `NPZ`, `Suppressor` and `Test`, which only the test suite needs. The four of them, plus `DelimitedFiles` and `QuadGK`, now live in the new `test/Project.toml`: since Julia 1.2 `Pkg.test` activates that file when it exists and adds the package under test to it, which is the portable form of the "workspace" layout (the literal `[workspace]` table of the root `Project.toml` is a Julia 1.12 feature, and the blocking CI job runs 1.9). `src/GaPSE.jl` drops `using Dierckx`, `using Test` and `using Documenter`: none of the three exported a name the library uses - `derivative` is `MySpline`'s own and is always called as `GaPSE.derivative` - so removing `Dierckx` also removes a latent name clash with it. `install_gapse.jl`, which activates the root environment and would have added them all back, was trimmed to the same list, and the dependency sections of `README.md` and `docs/src/index.md` now say which packages are the library's and which are only the tests' or the documentation's;

- replaced with four spaces the indentation tabs left in eight `src/GNCxLD_CrossCorrelations/*.jl` files (`GNCxLD_LocalGPLensing.jl`, `GNCxLD_LensingLocalGP.jl`, `GNCxLD_LensingLensing.jl`, `GNCxLD_LensingDoppler.jl`, `GNCxLD_LensingIntegratedGP.jl`, `GNCxLD_IntegratedGPDoppler.jl`, `GNCxLD_IntegratedGPLocalGP.jl`, `GNCxLD_DopplerLocalGP.jl`) and the stray one on a blank line of `src/OtherUtils.jl`, completing the sweep begun in the previous commit. The only tabs left in `src/` and `test/` are inside the comments that describe the tab-separated columns of an input file, where they are the content;

- IMPORTANT CHANGE: replaced [Dierckx](https://github.com/kbarbary/Dierckx.jl) with our own cubic spline `MySpline`, implemented in the new `src/Spline.jl`, for all the 1D interpolations of the library (`Cosmology`, `InputPS`, `IntegralIPS`, `IPSTools`, `BackgroundData`, `XiMatter`, ...). It supports the `"Natural"`, `"Parabolic"` and `"ThirdDerivative"` initial conditions and provides its own `derivative`. `Dierckx` is no longer used anywhere inside `src/`; it is kept only as a *test* dependency (`test/Project.toml`), where the suite uses it as an independent cross-check;

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

- KNOWN BUG, documented but not yet fixed: the `Δχ → 0` limits assume `y → 1`, which is what `Δχ → 0` forces at *fixed, non-zero* comoving distances. It does not hold in the small-χ corner, where `χ1` and `χ2` go to zero together at any `y`, and every double-χ grid starts at `χ = 1e-6 * s ≃ 4e-4`, so the whole first row and column of the grid take the branch with full weight. There the correct limits are `3 y σ_2 + 6/5 χ1^2 σ_0` (family 1) and `A/3 y σ_2` (family 8) instead of the same expressions with `y = 1`, so the seven double-χ integrands return a value a factor `1/y` too large, and of the wrong sign for `y < 0`. This is what moves `ξ_GNCxLD_Lensing_Lensing` by 2-4%, and `ξ_GNCxLD_Lensing_Lensing` is the only quantity that fails in `test_GNCxLD_SumXiMultipoles_P1.jl` (see the note below). The derivation, the measured effect and the proposed relative guard `Δχ < Δχ_min * max(χ1, χ2)` are in the "A second way to reach Δχ = 0" section of the manual;

- NOTE on the unit tests: `test_GNCxLD_SumXiMultipoles_P1.jl` reports 28 failing assertions out of 672, and every one of them is the `lensing_lensing` entry (index 10 of `GR_EFFECTS_GNCxLD`); no other effect, and no `res_sums_*`, ever fails. Both blocks of the file see the same thing: `map_sum_ξ_GNCxLD_multipole` tests the 20 effects one by one and only `[10]` fails, while `sum_ξ_GNCxLD_multipole` compares the whole 20-element vector at a single `s` and fails at exactly the matching `(L, s)`. The comparison with `064e23f` - same test file, same reference data, `MySpline` already in, but *before* the `Δχ → 0` limit branches - shows what the branches did: there were 20 failures then, all of them `no_window` (`map_sum` for every `L = 0 ... 4`, `sum_` at `s = 500`); now there are 4 `no_window` ones (`L = 3` only, at `s = 10`) and 24 `with_window` ones (every `L`, at `s = 1000`). So the limit branches removed 16 of the 20 pre-existing failures and introduced 24 new ones. Note that `ξ_lensing_lensing` is never more than about 1% of the sum of the 20 effects, which is why the total stays inside `rtol = 1.2e-2` while the single term is out;

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

- added the `s -> +infinity` counterpart to `docs/src/IlnIntegrals.md`, with proof: `I_l^n(s) -> A/(2 pi^2) M_l(3+n_P-n) s^-(3+n_P)`, where `P(q) -> A q^n_P` for `q -> 0` and `M_l(z) = sqrt(pi/2) 2^(z-3/2) Gamma((l+z)/2)/Gamma((l-z+3)/2)` is the Mellin transform of the spherical Bessel function. THE EXPONENT DEPENDS ON NEITHER `l` NOR `n`: the `s^-n` pulled out by the `(qs)^-n` factor exactly compensates the `n` it subtracts from `mu = 3 + n_P - n`. Only the amplitude does, and the limit is always 0. This is also why `IntegralIPS` seeds its right-hand fit with `p0 = [-4.0, 1.0]`: `3 + n_P = 3.96`;

- NOTE on `n_P`: it is the small-`k` slope of the late-time MATTER Power Spectrum, in `(h^-1 Mpc)^3`, and NOT the `n_s - 1` of the dimensionless primordial curvature spectrum `Delta^2_R(k) = A_s (k/k*)^(n_s-1)`. The two differ by `P_m(k) ~ k^4 T^2(k) P_R(k) ~ k^n_s T^2(k)`, with `T -> 1` on large scales, so `n_P = n_s ~ +0.96` and the matter `P(k)` GROWS at small `k`. Measured on `WideA_ZA_pk.dat`: `d ln P / d ln k = +0.9600` over `[1e-6, 1e-5]` and `+0.9599` over `[1e-5, 1e-4]`. `theory/Iln_terms.jl` reads it out of the `InputPS` left fit instead of assuming it, and the constant is now named `N_P` rather than `NS`;

- NOTE on the proof: the obvious route - substitute `x = q s` and replace `P(x/s)` with its small-`q` power law - is NOT legitimate, because the integral runs to `x = +infinity` and `x/s` is an indeterminate form there, with no dominating function. It is not an academic objection: what it leaves behind, `int_0^inf x^(mu-1) j_l(x) dx`, diverges for `mu >= 2`, i.e. for five of the eight `I_l^n`. The page now derives the result from the residue at the rightmost pole of the Mellin-Parseval representation, which never forms that integral and evaluates `M_l(z)` at `z = mu`, where it is analytic (its only poles are at `z = -l-2k`);

- `theory/Iln_terms.jl`/`.ipynb` gained `mellin`, `asymptote_large`, `plot_ratios_large_s` and `save_large_s_data`, plus the large-`s` asymptote on each single figure. The formula reproduces the stored `I_l^n` to better than 0.7% at `s = 1e4` for all eight, the five `mu > 2` cases included, which is the numerical confirmation that the analytic continuation is the right object;

- NOTE: `I~_0^4` is excluded from the large-`s` analysis. Its `mu = n_P - 1 = -0.04` sits on the pole of `M_0(z)` at `z = 0`, which is exactly the `sigma_4/s^4` divergence its subtraction removes; what is left has two powers, `s^-(3+n_P)` and `s^-4`, degenerate up to `1 - n_P = 0.04`, so it never settles on a clean power law (measured local slope still `-3.6` at `s = 1e3`). Note also that its `right` is `9888`, well inside the `Dchi ~ 2 chi_max` a real run reaches, unlike the `96466` of the other `I_l^n`;

- rewrote `docs/src/DeltaChiLimits.md` with step-by-step derivations for all eight families, instead of the terse "the bracket vanishes at the singular point" arguments of the first version. Each family now shows the exact algebraic decomposition of its vanishing brackets, the order-by-order substitution, the per-term limit and the cancellation of the direction parameter `p`. Every one of the eight limits was re-verified symbolically and all of them, and the values coded in `src/`, are confirmed correct;

- added the section "A trap: the two small parameters are NOT of the same order" to the same page. Along `chi2 = chi1 + p Dchi` one has `chi2 - chi1 = O(Dchi)` while `y - 1 = O(Dchi^2)`, so a term QUADRATIC in `(chi1 - chi2)` is the same order as a term LINEAR in `(y-1)`. Setting `chi1 = chi2` inside a bracket before expanding in `y` therefore silently loses a leading contribution. For the `J_00` bracket of the Lensing x Lensing family this is the difference between `8(chi1-chi2)^2 + 8(y-1)(chi1^2+chi2^2) - 9 chi1 chi2 (y^2-1) -> (1+7p^2) Dchi^2` and the wrong `(1-p^2) Dchi^2`, i.e. between `3/4 chi1^2 sigma_0 (1-p^2)(1+7p^2)` and `3/4 chi1^2 sigma_0 (1-p^2)^2`. The wrong form coincides with the right one at `p = 0` and `p = +-1`, so endpoint spot checks do not catch it - only the full `p` cancellation does. The page now carries this as an explicit warning box;

- corrected the vanishing orders quoted for two families: the `J_02` and `J_04` numerators of Lensing x Doppler vanish as `Dchi^3`, not `Dchi`; the `J_04` bracket of Newtonian x Lensing vanishes as `Dchi^4`, not `Dchi^5`, so that `J_04` there is finite rather than `O(Dchi)`. The limits are unchanged;

- added "A pattern worth noticing": the leading coefficients of these expansions are always low-order Legendre polynomials in the direction parameter, `35p^4-30p^2+3 = 8 L_4(p)` and `3p^2-1 = 2 L_2(p)`, which is a quick sanity check when redoing any of them;

- split the `Delta chi -> 0` limits documentation, which had grown to ~1400 lines, into nine pages: `docs/src/DeltaChiLimits.md` now holds only the shared material - definitions, how the limit is taken, the trap about the two orders, the Legendre pattern, the index of the eight families, the integrand-by-integrand summary table, the small-`chi` corner and the discussion of `Dchi_min` - while each family gets its own page, `DeltaChiLimits_1_LensingLensing.md` ... `DeltaChiLimits_8_J22J31.md`, registered as a nested entry of the Theory section in `docs/make.jl`;

- every family page opens with a "Recap: everything this page needs" section repeating the equations it uses - the definition of its own `Delta chi` and its singular point, the `I_l^n` limits (2.1a)/(2.1b) with the explicit table, the parametrisation (2.2) and its consequences (2.3a)/(2.3b), and the warning about the two orders - so that no page has to be read alongside another. Equations keep a single global numbering, the prefix identifying the page they belong to, and whenever a page quotes an equation derived elsewhere that equation is reproduced in full where it is used;

- rewrote the derivations in a more explicit style throughout: `### Term N: J I` headings instead of `### Step N`, every intermediate algebraic manipulation written out rather than summarised, and an explicit order-counting table for the `B_22` Taylor expansion. Where a decomposition is compared with one from another family, the latter is now restated on the spot instead of being referenced;

- restored the cross-page markdown links in the `Delta chi -> 0` documentation (the eight-family index, the Legendre-pattern table, the integrand summary, the small-`chi` corner table and the "derived in ..." lines of every recap), plus the two that link `IlnIntegrals.md` to `SphericalBesselFunctions.md` and back;

- filled in the derivation steps that were only asserted. In Family 1 the five Taylor coefficients of `B_22` around `y = 1` are now computed one by one: the `y` derivative is taken term by term and shown, the result is evaluated at `y = 1`, and the factorisation is done through the symmetric combinations `u = chi1^2 + chi2^2` and `v = chi1 chi2`, for which `chi1^4 + chi2^4 = u^2 - 2v^2` and, crucially, `u - 2v = (chi1 - chi2)^2`. Every `(u - 2v)` is then an `O(Dchi^2)`, so reading the order off a coefficient becomes a matter of counting its powers: `8(u-2v)^2`, `4(u-2v)(7u-2v)`, `2(7u^2-24uv+26v^2)` (which does not contain one, hence `O(1)`), `-4v(4u-11v)` and `11v^2`;

- same treatment for the numerators of Families 2 and 3, whose vanishing orders were previously quoted without proof: `N_02` and `N_04` of Lensing x Doppler are Taylor-expanded in `(y-1)`, their `k=0` coefficients factored as `4(chi1-s2)^3(chi1+2s2)` and half of it, their `k=1` ones as `2(chi1-s2)(...)`, the `k>=2` ones dropped by the order count, and the two surviving contributions summed to `-3 chi1 p (3p^2+1) Dchi1^3` and `3 chi1 p (3-5p^2) Dchi1^3`. Likewise `N_02^(b)`, `N_02^(f)` and `N_04` of Newtonian x Lensing, the last one with the same explicit order-counting table used for `B_22`. All the factorisations and derivatives were verified symbolically;

- every definition in the `Delta chi -> 0` documentation now carries a `(D.n)` number and is referenced by it throughout: `(D.1)`-`(D.3)` the three separations `Dchi`, `Dchi_1`, `Dchi_2`, `(D.4)` the direction of approach `chi2 := chi1 + p Dchi` (these four were previously numbered `(1.1)`-`(1.3)` and `(2.2)`, and are renumbered), then `(D.5)`-`(D.19)` for the symbols each family introduces - `B_00`, `B_02`, `B_22`, the symmetric pair `u`, `v`, the numerators `N_02`, `N_04`, `N_02^(b)`, `N_02^(f)`, the factors `F`, `J_20`, `G`, `J_22`, `J_31`. The introduction gained an "Index of the definitions" table listing all nineteen with a link to the page that uses them, and the numbering note now states that definitions carry the `D` prefix wherever they live, while everything else stays numbered by its page;

- applied to Families 3 to 8 of the `Delta chi -> 0` documentation the style settled on Families 1 and 2: each page now opens with the full GNC TPCF whose limit is being taken (`xi^{delta kappa}`, `xi^{kappa phi}`, `xi^{delta int phi}`, `xi^{int phi int phi}`, `xi^{v_par int phi}`, `xi^{kappa int phi}`), copied from the corresponding docstring, followed by its `J` coefficients carrying the same superscript, an explicit "The limit to be taken is" display and a note on which combinations the analysis covers; the boxed result at the bottom carries the superscripted `J` as well;

- the polynomial groupings are written out step by step everywhere, as in Family 2: the coefficient split (e.g. `-10 = -6-4`, `11 = 3+8`), the regrouping into two bracketed sums, the factoring out of the common linear factor and only then the result; the derivatives are shown as `d/dy[<full expression>]` before being evaluated; each Taylor coefficient is followed by its order; and every expansion restates `N = sum_k (1/k!) d^kN/dy^k|_(y=1) (y-1)^k`, then the surviving terms symbolically, then the explicit factored forms, then the `p` substitution;

- the final-sum tables were replaced by `align*` blocks, each row tagged with the equation number it comes from;

- BUG FIX in the docstring of `integrand_ξ_GNC_Newtonian_Lensing`: the `b_1` part of `J^{delta kappa}_{02}` is written as `[(9y^2+11) f_1 - 7(y^2+3) b_1] s_1^2 chi_2`, with a minus sign, while the code has a plus; the difference is `14 b_1 chi_2 s_1^2 (y^2+3)`. It matters: with the docstring sign the bracket does NOT vanish at the singular point (it gives `-5 s_1^2`), so `J_02` would diverge as `Dchi_2^-2` and the limit of the family would pick up an extra `-b_1 s_1^3 sigma_-2/6`. The code is the correct one; the discrepancy is flagged in a warning box on the Family 3 page, and the docstring still needs fixing;

- DOCSTRING FIX: the docstrings of `integrand_ξ_GNC_Newtonian_Lensing` and `integrand_ξ_GNCxLD_Newtonian_Lensing` wrote the `b_1` part of `J^{delta kappa}_{02}` as `[(9y^2+11) f_1 - 7(y^2+3) b_1] s_1^2 chi_2`, with a minus sign, while the code has a plus (4 occurrences, 2 per file, each docstring appearing once per method signature). The difference is `14 b_1 chi_2 s_1^2 (y^2+3)`: with the minus the bracket does NOT vanish at the singular point (it gives `-5 s_1^2`), `J_02` would diverge as `Dchi_2^-2` and the `Dchi -> 0` limit of that family would pick up a spurious `-b_1 s_1^3 sigma_-2/6`. Only the docstrings were wrong; the code is untouched. Every other `J` reproduced in the `Delta chi -> 0` pages - Families 1, 2, 4, 5, 7 and 8 - was compared symbolically against its code and matches;


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
