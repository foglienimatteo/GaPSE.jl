# theory

Julia scripts and notebooks used to *investigate* the theoretical behaviour of the
quantities GaPSE computes, and to produce the figures shown in the **Theory** section of
the documentation.

This directory is not part of the GaPSE library: nothing under `src/` depends on it, and
it is not exercised by the test suite. It is the place where to put a self-contained piece
of analysis, the plots it produces and the data behind them, so that a result shown in the
manual can be reproduced and re-checked later.

Each analysis is made of

- `<name>.jl` : the script, runnable from the terminal;
- `<name>.ipynb` : the same computation as a notebook, with the theory written out;
- `<name>/` : the directory where that analysis saves its plots and data.

## Setup

This directory has its own environment, so that the plotting packages are not added to
the dependencies of GaPSE itself. The first time, from the `theory/` directory:

```julia
julia> using Pkg

julia> Pkg.activate(".")

julia> Pkg.develop(path = "..")   # use the GaPSE of this repository

julia> Pkg.instantiate()
```

Then a script is run simply with

```bash
$ julia --project=. Iln_terms.jl
```

## Contents

- **`Iln_terms`** : the ``I_\ell^n`` integrals that build every TPCF, plotted in log-log
  scale together with their small-``s`` asymptotes. It produces the figures used by the
  "The ``I_\ell^n`` integrals" page of the documentation, and, with `SAVE_TO_DOCS = true`,
  it writes a copy of them directly into `docs/src/assets/Iln_terms/`.

- **`sigma_i`** : the moments ``\sigma_i = \int \mathrm{d}q \, q^{2-i} P(q) / 2\pi^2`` of
  the input Power Spectrum, which is what every ``\Delta\chi \rightarrow 0`` limit of the
  TPCFs reduces to. It plots the five integrands and the fraction of each moment already
  collected below a given ``q``, and prints the moments over a grid of integration
  extremes, so that one can see which of them converge and which ones are defined by the
  cut itself. It also compares the two ranges GaPSE actually uses - the `k_min`/`k_max`
  given to `IPSTools`, which fix the stored ``\sigma_i``, and the `1e-5, 1e3` that
  `IPSTools` hard-codes for the `xicalc` call that builds the ``I_\ell^n``. Figures go to
  `docs/src/assets/sigma_i/`.

- **`spherical_bessels`** : the spherical Bessel functions ``j_\ell(x)``, the region where
  their small-``x`` expansion ``x^\ell/(2\ell+1)!!`` is legitimate, and the first zero of
  ``j_\ell`` as a function of ``\ell``. It produces the figures used by the
  "Spherical Bessel Functions" page of the documentation, writing them into
  `docs/src/assets/misc/`.
