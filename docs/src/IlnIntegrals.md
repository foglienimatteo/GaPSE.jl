# The ``I_\ell^n`` integrals

Every Two-Point Correlation Function (TPCF) that GaPSE computes is, in the end, a sum of
terms of the form ``J(\chi, s, y) \, I_\ell^n(\Delta\chi)``. This page collects the
definition of these ``I_\ell^n``, proves their behaviour for small separations, and shows
what they look like.

The plots are produced by the script `theory/Iln_terms.jl` (or, equivalently, by the
notebook `theory/Iln_terms.ipynb`); see the end of this page.

## Definitions

The integrals are

```math
    I_\ell^n(s) = \int_0^{+\infty} \frac{\mathrm{d}q}{2\pi^2} \, q^2 \, P(q) \,
        \frac{j_\ell(qs)}{(qs)^n} \; ,
```

with ``P(q)`` the matter Power Spectrum at ``z=0`` stored inside the `Cosmology`, and
``j_\ell`` the spherical Bessel function of order ``\ell``. The eight combinations
``(\ell, n)`` that appear in the code are

```math
    (0,0) \, , \quad (2,0) \, , \quad (4,0) \, , \quad (0,2) \, , \quad
    (2,2) \, , \quad (3,1) \, , \quad (1,3) \, , \quad (1,1) \; ,
```

and they are stored in `IPSTools` as `I00`, `I20`, `I40`, `I02`, `I22`, `I31`, `I13` and
`I11` respectively, i.e. the name is `I` followed by ``\ell`` and then ``n``. Alongside
them there is the auxiliary integral

```math
    \tilde{I}_0^4(s) = \frac{1}{s^4} \int_0^{+\infty}
        \frac{\mathrm{d}q}{2\pi^2} \, \frac{P(q)}{q^2} \, \left[ j_0(qs) - 1 \right] \; ,
```

stored as `I04_tilde` and computed by [`GaPSE.func_I04_tilde`](@ref).

We also need the moments of the Power Spectrum

```math
    \sigma_i = \int_{k_\mathrm{min}}^{k_\mathrm{max}}
        \frac{\mathrm{d}q}{2 \pi^2} \, q^{2-i} \, P(q) \; ,
```

of which `IPSTools` stores ``\sigma_0, \, \sigma_1, \, \sigma_2, \, \sigma_3`` and
``\sigma_4``. Note that ``\sigma_i`` with ``i<0`` also appear below; they are perfectly
finite as long as ``k_\mathrm{max}`` is finite, but they are not stored, so the script
recomputes them when needed.

## The small-``s`` behaviour

!!! note "Result"
    ```math
        I_\ell^n(s) \; \xrightarrow[s \rightarrow 0]{} \;
            \frac{\sigma_{n-\ell}}{(2\ell+1)!!} \, s^{\,\ell - n} \; ,
        \qquad
        \tilde{I}_0^4(s) \; \xrightarrow[s \rightarrow 0]{} \; - \frac{\sigma_2}{6 \, s^2} \; .
    ```

**Proof.** The spherical Bessel function of order ``\ell`` has the everywhere-convergent
Taylor series

```math
    j_\ell(x) = \sum_{k=0}^{+\infty} \frac{(-1)^k \, x^{\ell + 2k}}
        {2^k \, k! \, (2\ell + 2k + 1)!!} \; ,
```

as can be checked for ``\ell = 0``, where it reproduces
``j_0(x) = \sin(x)/x = 1 - x^2/6 + x^4/120 - \dots`` Dividing by ``x^n`` and setting
``x = qs``,

```math
    \frac{j_\ell(qs)}{(qs)^n} = \sum_{k=0}^{+\infty} \frac{(-1)^k \, (qs)^{\ell + 2k - n}}
        {2^k \, k! \, (2\ell + 2k + 1)!!} \; .
```

Inserting this into the definition of ``I_\ell^n`` and exchanging the sum with the
integral, which is legitimate because the series converges uniformly on the compact
integration range ``[k_\mathrm{min}, k_\mathrm{max}]``,

```math
    I_\ell^n(s) = \sum_{k=0}^{+\infty} \frac{(-1)^k \; s^{\,\ell - n + 2k}}
        {2^k \, k! \, (2\ell + 2k + 1)!!}
        \int \frac{\mathrm{d}q}{2\pi^2} \, q^{\,2 + \ell + 2k - n} \, P(q) \; .
```

The remaining integral is by definition ``\sigma_i`` with ``2 - i = 2 + \ell + 2k - n``,
that is ``i = n - \ell - 2k``, so that

```math
    I_\ell^n(s) = \sum_{k=0}^{+\infty} \frac{(-1)^k \; \sigma_{n - \ell - 2k}}
        {2^k \, k! \, (2\ell + 2k + 1)!!} \; s^{\,\ell - n + 2k} \; .
```

Every term carries two more powers of ``s`` than the previous one, so for
``s \rightarrow 0`` the ``k=0`` term dominates, which is the claimed result.

For ``\tilde{I}_0^4`` the same argument applies to ``j_0(x) - 1``, whose series starts at
``k=1``:

```math
    j_0(x) - 1 = \sum_{k=1}^{+\infty} \frac{(-1)^k \, x^{2k}}{(2k+1)!} \; ,
```

where we used ``2^k \, k! \, (2k+1)!! = (2k+1)!``. Then

```math
    \tilde{I}_0^4(s) = \frac{1}{s^4} \sum_{k=1}^{+\infty} \frac{(-1)^k \, s^{2k}}{(2k+1)!}
        \int \frac{\mathrm{d}q}{2\pi^2} \, q^{\,2k-2} \, P(q)
    = \sum_{k=1}^{+\infty} \frac{(-1)^k \, \sigma_{4-2k}}{(2k+1)!} \, s^{\,2k-4} \; ,
```

whose first terms are ``-\sigma_2/(6 s^2) + \sigma_0/120 - \sigma_{-2} s^2/5040 + \dots``
``\blacksquare``

### The three regimes

The sign of ``\ell - n`` decides everything:

- ``\ell > n`` : the integral **vanishes** as ``s^{\ell-n}``;
- ``\ell = n`` : the integral tends to the **finite**, non-zero value
  ``\sigma_0/(2\ell+1)!!``;
- ``\ell < n`` : the integral **diverges** as ``s^{-(n-\ell)}``.

| | ``\ell`` | ``n`` | ``s \rightarrow 0`` | regime |
|:--|:-:|:-:|:--|:--|
| ``I_0^0`` | 0 | 0 | ``\sigma_0`` | finite |
| ``I_2^0`` | 2 | 0 | ``\sigma_{-2} \, s^2 / 15`` | vanishing |
| ``I_4^0`` | 4 | 0 | ``\sigma_{-4} \, s^4 / 945`` | vanishing |
| ``I_0^2`` | 0 | 2 | ``\sigma_2 / s^2`` | diverging |
| ``I_2^2`` | 2 | 2 | ``\sigma_0 / 15`` | finite |
| ``I_3^1`` | 3 | 1 | ``\sigma_{-2} \, s^2 / 105`` | vanishing |
| ``I_1^3`` | 1 | 3 | ``\sigma_2 / (3 s^2)`` | diverging |
| ``I_1^1`` | 1 | 1 | ``\sigma_0 / 3`` | finite |
| ``\tilde{I}_0^4`` | - | - | ``-\sigma_2 / (6 s^2)`` | diverging |

The diverging ones are never a problem in practice: inside a TPCF they always come
multiplied by a ``J`` carrying the matching positive power of ``\Delta\chi``, so that the
product stays finite. How the cancellation works, term by term and for every TPCF, is the
subject of [The ``\Delta\chi \rightarrow 0`` limits](@ref).

## The plots

All the curves show ``|I_\ell^n(s)|`` rather than ``I_\ell^n(s)``: these integrals
oscillate and change sign at large ``s``, where a logarithmic vertical axis would not be
defined. At small ``s``, the regime this page is about, they keep a constant sign, so the
absolute value is immaterial there and the asymptote (dashed, black) can be read off
directly.

```@raw html
<img src="assets/Iln_terms/all_Iln.png" alt="all the I_l^n" width="100%"/>
```

The three regimes are visible at a glance: the curves that flatten out
(``\ell = n``), the ones that fall off as a power law (``\ell > n``) and the ones that
blow up (``\ell < n``).

### One by one

| | |
|:-:|:-:|
| ![I00](assets/Iln_terms/I00.png) | ![I20](assets/Iln_terms/I20.png) |
| ![I40](assets/Iln_terms/I40.png) | ![I02](assets/Iln_terms/I02.png) |
| ![I22](assets/Iln_terms/I22.png) | ![I31](assets/Iln_terms/I31.png) |
| ![I13](assets/Iln_terms/I13.png) | ![I11](assets/Iln_terms/I11.png) |
| ![I04_tilde](assets/Iln_terms/I04_tilde.png) | |

## Reproducing the figures

The figures are not built by the documentation: they are committed under
`docs/src/assets/Iln_terms/`, and regenerated on demand by the script in the `theory/`
directory. From `theory/`, after the one-time setup described in its `README.md`:

```bash
$ julia --project=. Iln_terms.jl
```

which writes the plots and a table of the numerical values in `theory/Iln_terms/`, and,
since `SAVE_TO_DOCS = true`, also refreshes the copies used by this page. The same
computation is available step by step in the notebook `theory/Iln_terms.ipynb`.

See also: [`GaPSE.IPSTools`](@ref), [`GaPSE.IntegralIPS`](@ref),
[`GaPSE.func_I04_tilde`](@ref).
