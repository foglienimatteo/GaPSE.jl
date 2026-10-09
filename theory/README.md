# theory

Julia scripts and notebooks used to *investigate* the theoretical behaviour of the
quantities GaPSE computes, and to produce the figures shown in the **Theory** section of
the documentation.

This directory is not part of the GaPSE library: nothing under `src/` depends on it, and
it is not exercised by the test suite. It is the place where to put a self-contained piece
of analysis, the plots it produces and the data behind them, so that a result shown in the
manual can be reproduced and re-checked later.

Each analysis is made of

- `<name>.ipynb` : the notebook holding both the computation and the theory behind
  it. It is the only source - there is no parallel `.jl` script to keep in sync;
- `<name>/` : the directory where that analysis saves its plots and data.

## Setup

This directory has its own environment, so that the plotting packages are not added to
the dependencies of GaPSE itself. `Project.toml` is tracked by git and declares the GaPSE of
this repository through a `[sources]` entry, so the first time, from the `theory/` directory:

```julia
julia> using Pkg

julia> Pkg.activate(".")

julia> Pkg.instantiate()          # downloads and installs everything, GaPSE included
```

`Manifest.toml` is gitignored and is rebuilt by `instantiate` on each machine.
`../ipynbs/README.md` explains how that `Project.toml` was built, and why a notebook needs
nothing more than `using GaPSE` to find it.

Then a notebook is opened simply with

```bash
$ jupyter lab Iln_terms.ipynb
```

## Contents

- **`deltachi_limits`** : the `Δχ_min` where a ``\chi``-integrated TPCF switches from the
  ``J \cdot I_\ell^n`` sum to its analytic ``\Delta\chi \rightarrow 0`` limit, and how
  visible that switch is. It plots both branches for the two GNC auto-correlations,
  Lensing-Lensing and IntegratedGP-IntegratedGP, explains the step of the first one
  through the ``k_\mathrm{max}`` dependence of ``\sigma_0``, and closes by measuring how
  much the integrated TPCF really moves when `Δχ_min` is changed by four orders of
  magnitude (it does not).

- **`Iln_terms`** : the ``I_\ell^n`` integrals that build every TPCF, plotted in log-log
  scale together with their small-``s`` asymptotes. It produces the figures used by the
  "The ``I_\ell^n`` integrals" page of the documentation, and, with `SAVE_TO_DOCS = true`,
  it writes a copy of them directly into `docs/src/assets/Iln_terms/`.

- **`input_ps`** : the input matter Power Spectrum itself, plotted over twelve decades
  together with the two power laws `InputPS` extrapolates with outside the tabulated
  range, plus its local slope and a table of which ``\sigma_i`` converge. It produces
  the figures of the "The input Power Spectrum" page of the documentation, writing them
  into `docs/src/assets/input_ps/`.

- **`sigma_i`** : the moments ``\sigma_i = \int \mathrm{d}q \, q^{2-i} P(q) / 2\pi^2`` of
  the input Power Spectrum, which is what every ``\Delta\chi \rightarrow 0`` limit of the
  TPCFs reduces to. It plots the five integrands and the fraction of each moment already
  collected below a given ``q``, and prints the moments over a grid of integration
  extremes, so that one can see which of them converge and which ones are defined by the
  cut itself. It also compares the two ranges GaPSE actually uses - the `k_min`/`k_max`
  given to `IPSTools`, which fix the stored ``\sigma_i``, and the `1e-5, 1e3` that
  `IPSTools` hard-codes for the `xicalc` call that builds the ``I_\ell^n``. Figures go to
  `docs/src/assets/sigma_i/`.

- **`kmin_kmax`** : what the ``[k_\mathrm{min}, k_\mathrm{max}]`` pair of `IPSTools`
  does to the five ``\sigma_i`` and to the nine ``I_\ell^n``, which it now bounds
  together (it used to bound only the ``\sigma_i``, the `xicalc` behind the
  ``I_\ell^n`` being hard-coded to `1e-5, 1e3`). It scans each extreme in turn, shows
  that the ``I_\ell^n`` only move below ``s \simeq 2 \; h^{-1}``Mpc, and measures the
  one fragile piece: the left power-law fit collapses when the fit window straddles the
  saturation scale, i.e. unless `fit_min` ``\gg 1/k_\mathrm{max}``. It closes on three
  cosmologies that differ only by ``k_\mathrm{max}``, and the GNC Lensing-Lensing
  multipole they give.

- **`spherical_bessels`** : the spherical Bessel functions ``j_\ell(x)``, the region where
  their small-``x`` expansion ``x^\ell/(2\ell+1)!!`` is legitimate, and the first zero of
  ``j_\ell`` as a function of ``\ell``. It produces the figures used by the
  "Spherical Bessel Functions" page of the documentation, writing them into
  `docs/src/assets/misc/`.

- **`spline_comparison`** : `MySpline`, the cubic spline of `src/Spline.jl`, against
  [Dierckx](https://github.com/kbarbary/Dierckx.jl), in accuracy and in speed, on
  ``j_0(x)``, on the input ``P(q)`` and on ``I_0^0(s)``. Each function is plotted with
  the two interpolants on top of each other and their ratio underneath, the error
  against the exact ``j_0`` separates the three `ic` options, and the benchmarks time
  both the construction and the evaluation from ``N = 10^2`` to ``N = 10^5``.
