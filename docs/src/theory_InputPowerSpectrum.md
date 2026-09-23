# The input Power Spectrum

Everything GaPSE computes starts from a tabulated matter Power Spectrum ``P(q)``, read
by `InputPS`. Its two **asymptotic power laws** — the one at small ``q`` and the one at
large ``q`` — are not a detail: they decide whether each moment

```math
    \sigma_i = \int_{k_\mathrm{min}}^{k_\mathrm{max}} \frac{\mathrm{d}q}{2\pi^2}
    \, q^{\,2-i} \, P(q)
    \quad \quad (1)
```

is a number, or is set by where the integral happens to be cut. Since the
``\Delta\chi \rightarrow 0`` limits of every TPCF reduce to combinations of ``\sigma_i``
(see [The ``\Delta\chi \rightarrow 0`` limits](theory_DeltaChiLimits.md)), and the
``I_\ell^n`` reduce to them as ``s \rightarrow 0``
(see [The ``I_\ell^n`` integrals](theory_IlnIntegrals.md)), the answer matters.

```@contents
Pages = ["theory_InputPowerSpectrum.md"]
Depth = 3
```


## The metter Power Spectrum at present day


In the ``\Lambda``-CDM cosmology, the primordial curvature perturbations are a pure power law (see the [Planck 2018 results, A&A 641, A10 (2020)](https://doi.org/10.1051/0004-6361/201833887)), usually quoted through the dimensionless power spectrum ``\Delta^2_{\mathcal{R}}``:

```math
    \Delta^2_{\mathcal{R}}(k) := \frac{k^3}{2\pi^2} P_\mathcal{R}(k) \quad \quad (\mathrm{P}.2)
```

What enters the ``I_\ell^n`` is the matter Power Spectrum ``P_m(k,z)`` at present day (``z=0``), a dimensional quantity in ``(h^{-1}\mathrm{Mpc})^3``. The two are related by the Poisson equation and the transfer function ``T(k)``:

```math
    P_m(k, z)  \propto  k^4 \, T^2(k) \, D^2(z) \, P_\mathcal{R}(k)
    \quad \quad (\mathrm{P}.5)
```

The transfer function is asymptotic to ``k^{-2}\ln k`` at small scales, and constant at large ones It is normalized such that at large scales it goes to ``1``:

```math
    T(k) \xrightarrow[k \rightarrow 0^{+}]{} 1 \quad \quad (\mathrm{P}.6) \\[10pt]
    T(k) \underset{k \rightarrow +\infty}{\sim} k^{-2}\ln k
```



## The small-``k`` slope

At large scales, it's known that:

```math
\begin{align*}
    \Delta^2_{\mathcal{R}}(k) &\underset{k \rightarrow 0^{+}}{\sim} \; A_s \, \left( \frac{k}{k^*}\right)^{n_s-1} &&(\mathrm{P}.1)\\[10pt]

    n_s&\approx 0.965&&\mathrm{primordial \; spectral \; index}\\[8pt]
    \ln(10^{10} A_s) &\approx 3.043 &&\mathrm{ log \; power \; of \; primordial \; curvature \; perturbations} \Rightarrow A_s \approx 2.1 \times 10^{-9}\\[8pt]
    k^* &\approx 0.05 \; \mathrm{Mpc}^{-1}  &&\mathrm{arbitrary \; pivot \; scale}
\end{align*}
```

```math
    \Rightarrow P_\mathcal{R}(k) = \frac{2\pi^2}{k^3}\,\Delta^2_{\mathcal{R}}(k) \underset{k \rightarrow 0^{+}}{\sim} k^{\,n_s-4} \quad \quad (\mathrm{P}.2) \\[10pt]
```

First of all, we start from the Poisson equation in real comoving space

```math
\begin{align*}
    \nabla^2 \phi(\mathbf{s}, z) = - \frac{4 \pi G}{c^2}\langle\rho\rangle D^2(z) \, \Delta(\mathbf{s}) \, ,
\end{align*}
```

The transfer function is normalized such that at large scales it goes to ``1``:
```math
    T(k) \xrightarrow[k \rightarrow 0^{+}]{} 1 \quad \quad (\mathrm{P}.6) \\[10pt]
```

so at present day:

```math
\begin{align*}
    (\mathrm{P}.5)\; \mathrm{with} \; z=0 : \quad \quad P_m(k, z=0)  & \propto  k^4 \, T^2(k) \, D^2(0) \, P_\mathcal{R}(k)\\[10pt]
    \mathrm{Inserting} \; (\mathrm{P}.6) \rightarrow  \quad \quad
    &\underset{ k\rightarrow 0^{+}}{\sim} \; k^{4}  \; P_\mathcal{R}(k) \\[10pt]
    \mathrm{Inserting} \; (\mathrm{P}.4) \rightarrow  \quad \quad
    &\underset{ k\rightarrow 0^{+}}{\sim} \; k^{4}  \, k^{\,n_s-4}  \\[10pt]
    &= \; k^{n_s}\\[10pt]
\end{align*}
```

```math
    \quad \Rightarrow \quad
    P(k) \; \underset{k \rightarrow 0^{+}}{\sim}  \, k^{n_s} \; ,
    \quad \quad n_s \simeq 0.96
    \quad \quad (\mathrm{P}.7)
```

The matter Power Spectrum therefore *grows* as ``k^{+0.96}`` at small ``k``. 
The dimensionless curvature goes as ``k^{n_s-1} \simeq k^{-0.035}`` and the four powers of ``k`` supplied by Poisson, minus the three of the ``\Delta^2 \leftrightarrow P`` conversion, are exactly what separates them.

NOTE: this is directly visible in `data/WideA_ZA_pk.dat`, whose local slope ``\mathrm{d}\ln P / \mathrm{d}\ln k`` is ``+0.9600`` over the first decade of the tabulated range, with an amplitude of order ``10^{6} \, (h^{-1}\mathrm{Mpc})^3``.



## The large-``k`` tail

At the other end the transfer function is no longer ``1``. For CDM it falls as
``T(k) \propto k^{-2}\ln k``, so

```math
    P_m(k) \; \underset{k \rightarrow +\infty}{\sim} \; k^{4} \, T^2(k) \, k^{\,n_s-4}
    \; \sim \; k^{\,n_s - 4} \, \ln^2 k \; \simeq \; k^{-3} \, \ln^2 k
    \quad \quad (\mathrm{P}.8)
```

That ``k^{-3}`` is the number that matters for (P.1), and it is **not** what GaPSE
actually uses. `InputPS` continues the tabulated spectrum with a power law
``P(q) = a + b \, q^{\,s}`` fitted on `[fit_right_min, fit_right_max]`, and on
`data/WideA_ZA_pk.dat` that fit gives

```math
    P(q) \; \underset{q \, > \, 20.2}{=} \; 91.60 \; q^{-2.641}
    \quad \quad (\mathrm{P}.9)
```

i.e. a tail *shallower* than the asymptotic ``k^{-3}``, simply because at
``k \simeq 20 \, h\,\mathrm{Mpc}^{-1}`` the true spectrum has not reached its asymptote yet.

## Which ``\sigma_i`` converge

Inserting ``P \sim q^{\,s}`` into (P.1), the integrand goes as ``q^{\,2-i+s}``. The
integral converges at ``q \rightarrow 0`` when ``2-i+s > -1``, and at
``q \rightarrow +\infty`` when ``2-i+s < -1``. With ``s = +0.960`` on the left and
``s = -2.641`` on the right:

| | IR exponent | IR | UV exponent | UV |
|:--|--:|:--|--:|:--|
| ``\sigma_0`` | ``+2.960`` | converges | ``-0.641`` | **diverges** |
| ``\sigma_1`` | ``+1.960`` | converges | ``-1.641`` | converges |
| ``\sigma_2`` | ``+0.960`` | converges | ``-2.641`` | converges |
| ``\sigma_3`` | ``-0.040`` | converges | ``-3.641`` | converges |
| ``\sigma_4`` | ``-1.040`` | **diverges** | ``-4.641`` | converges |

Two of the five do not exist as numbers:

- **``\sigma_0`` diverges in the ultraviolet.** With the true tail (P.8) the exponent
  would be exactly ``-1``, i.e. ``\sigma_0`` would diverge *logarithmically*; with the
  fitted tail (P.9) it diverges as a power,
  ``\sigma_0(<K) \sim K^{\,0.359}``. This is not a defect of the code: ``\sigma_0`` is
  the density variance ``\langle \delta^2 \rangle`` smoothed on *zero* scale, which is
  genuinely infinite in CDM. What makes it finite in GaPSE is ``k_\mathrm{max}``.
- **``\sigma_4`` diverges in the infrared**, as ``\sigma_4(>k) \sim k^{-0.040}``, i.e.
  slowly but without bound.

Measured on `data/WideA_ZA_pk.dat`, with ``k_\mathrm{min} = 10^{-5}``:

| ``k_\mathrm{max}`` | ``1`` | ``10`` | ``10^2`` | ``10^3`` | ``10^4`` |
|:--|--:|--:|--:|--:|--:|
| ``\sigma_0`` | ``3.41`` | ``18.6`` | ``56.5`` | ``143`` | ``341`` |
| ``\sigma_2`` | ``98.8`` | ``101.06`` | ``101.13`` | ``101.13`` | ``101.13`` |

``\sigma_2`` has converged by ``k \simeq 1``; ``\sigma_0`` has not converged anywhere.

### What this means for the ``\Delta\chi \rightarrow 0`` limits

``\sigma_0`` enters five limit branches — the three of the Lensing-Lensing family
(``3\sigma_2 + \frac{6}{5}\chi_1^2\sigma_0``) and the two of Newtonian ``\times``
Lensing (``-s_1(f_1+5b_1)\sigma_0/5``). The derivation of those limits assumed
``I_0^0(s) \rightarrow \sigma_0 = \mathrm{const}``, which is exactly the step that
fails.

What actually happens is that ``j_0(qs)`` cuts the ``q`` integral at ``q \sim 1/s``, so

```math
    I_0^0(s) \; \underset{s \rightarrow 0^{+}}{\simeq} \; \sigma_0(< 1/s)
    \; \sim \; s^{-0.359} \; ,
```

as the direct integration confirms — the ratio ``I_0^0(s)/\sigma_0(<1/s)`` is
``1.28``, ``1.22``, ``0.97``, ``1.00``, ``1.00`` at ``s = 10^{-1} \ldots 10^{-5}``, the
saturation at the end being only the effect of the ``k_\mathrm{max} = 10^3`` cut.

So ``I_0^0`` has **no** finite ``s \rightarrow 0`` limit, and the ``\sigma_0`` appearing
in those five branches is the ``\sigma_0`` regulated at some ``k_\mathrm{max}``. The
consistent choice is ``k_\mathrm{max} \simeq 1/\Delta\chi_\mathrm{min}``: the limit
branch replaces ``I_0^0(\Delta\chi)`` for ``\Delta\chi < \Delta\chi_\mathrm{min}``, and
``I_0^0(\Delta\chi_\mathrm{min}) \simeq \sigma_0(< 1/\Delta\chi_\mathrm{min})``. With the
default ``\Delta\chi_\mathrm{min} = 0.1`` that is ``k_\mathrm{max} = 10``, which is the
`IPSTools` default — so the code is, by luck rather than design, close to consistent:
``I_0^0(0.1) = 23.7`` against ``\sigma_0(<10) = 18.6``.

!!! warning "Do not unify the σ_i onto the xicalc range"
    `IPSTools` hard-codes ``[10^{-5}, 10^{3}]`` for the `xicalc` call that builds the
    ``I_\ell^n``, while the stored ``\sigma_i`` use its `k_min`/`k_max` keywords.
    Making the ``\sigma_i`` use the `xicalc` range too would give
    ``\sigma_0 = 143`` against ``I_0^0(0.1) = 23.7``, i.e. a factor ``6`` jump at the
    branch boundary instead of the present ``1.28``. For the *converged* moments
    (``\sigma_2``, ``\sigma_3``) the choice is irrelevant; for ``\sigma_0`` the range
    that matters is the one set by ``\Delta\chi_\mathrm{min}``, not the one used by
    `xicalc`.

## The figures

```@raw html
<img src="../assets/input_ps/input_ps.png" alt="The input matter Power Spectrum"/>
```

``P(q)`` over twelve decades. Outside the tabulated range (grey lines) the solid curve
*is* the dashed power law, which is the point of the figure: what GaPSE integrates
beyond ``q \simeq 20`` is the extrapolation (P.9), not data.

```@raw html
<img src="../assets/input_ps/input_ps_slope.png" alt="Local slope of the input Power Spectrum"/>
```

The same information as a local slope ``\mathrm{d}\ln P/\mathrm{d}\ln q``: flat at
``+0.960`` on the left, flat at ``-2.641`` on the right, with the turnover around the
equality scale in between. The dotted line marks the ``-3`` of (P.8), which the fitted
tail does not reach.

## Reproducing the figures

Both figures, and the tables above, are produced by `theory/input_ps.jl` (or the
equivalent `theory/input_ps.ipynb`):

```bash
$ cd theory
$ julia --project=. input_ps.jl
```

The companion analysis `theory/sigma_i.jl` plots the five ``\sigma_i`` integrands and
the fraction of each moment collected below a given ``q``, and lets any pair of
integration extremes be tried.
