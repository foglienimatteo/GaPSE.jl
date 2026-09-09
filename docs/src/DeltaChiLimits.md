# The ``\Delta\chi \rightarrow 0`` limits

All the Two-Point Correlation Functions (TPCFs) that require an integration over one or two
comoving distances ``\chi`` share the same structure: their integrand is a sum of terms

```math
    \sum_k J^{(k)}(\chi, s, y) \; I_{\ell_k}^{n_k}(\Delta\chi) \; ,
```

where ``\Delta\chi`` is the distance between the two points that are being correlated and the
``J^{(k)}`` are rational functions of the comoving distances and of the angle cosine ``y``.
Several of these ``I_\ell^n`` diverge for ``\Delta\chi \rightarrow 0``, and several of the
``J^{(k)}`` diverge as well: individually the terms are singular, while their sum is finite.

This page derives that finite value for every affected TPCF. The results are the ones used in
the `Δχ < Δχ_min` branch of the corresponding `integrand_ξ_...` functions.

## Why the limit is needed

``\Delta\chi = 0`` is not a pathological configuration that can be avoided: it is reached
deterministically by the quadrature.

- ``\Delta\chi^2 = \chi_1^2 + \chi_2^2 - 2\,\chi_1\chi_2\,y`` vanishes if and only if
  ``y = 1`` **and** ``\chi_1 = \chi_2``, since ``\Delta\chi^2 = (\chi_1-\chi_2)^2 + 2\chi_1\chi_2(1-y)``
  and both terms are non-negative for ``|y| \leq 1``.
- ``y = 1`` exactly whenever ``\mu = 1``, because
  ``y(s_1, s, \mu) = (\mu s + s_1)/s_2`` and ``s_2(s_1,s,1) = s_1 + s``. The nodes ``\mu = \pm 1``
  **belong** to the Gauss-Lobatto grid (`alg = :lobatto`) and to the trapezoidal grid
  (`alg = :trap`); only `alg = :quad` avoids them.
- ``\chi_1 = \chi_2`` is hit whenever a sampling point of the ``\chi`` grid coincides with the
  other distance. With `suit_sampling = true` this is not even accidental: 
  `sample_subdivision_middle` deliberately places a dense, symmetric sub-grid around the
  singular point, because that is where the integrand varies fastest.

## The ``\sigma_i`` and the small-argument expansion of ``I_\ell^n``

We recall the definitions used throughout GaPSE:

```math
    I_\ell^n(s) = \int_0^{+\infty} \frac{\mathrm{d}q}{2\pi^2} \, q^2 \, P(q) \,
        \frac{j_\ell(qs)}{(qs)^n} \; , \qquad
    \sigma_i = \int_0^{+\infty} \frac{\mathrm{d}q}{2\pi^2} \, q^{2-i} \, P(q) \; ,
```

with ``P(q)`` the matter Power Spectrum at ``z=0`` and ``j_\ell`` the spherical Bessel function of
order ``\ell``. Its Taylor series is

```math
    j_\ell(x) = \sum_{k=0}^{+\infty} \frac{(-1)^k \, x^{\ell + 2k}}
        {2^k \, k! \, (2\ell + 2k + 1)!!} \; ,
```

so that, inserting it in the definition and integrating term by term,

```math
    I_\ell^n(s) = \sum_{k=0}^{+\infty} \frac{(-1)^k \; \sigma_{n - \ell - 2k}}
        {2^k \, k! \, (2\ell + 2k + 1)!!} \; s^{\,\ell - n + 2k} \; .
```

The leading behaviour is therefore

```math
    I_\ell^n(s) \; \xrightarrow[s \rightarrow 0]{} \;
        \frac{\sigma_{n-\ell}}{(2\ell+1)!!} \, s^{\,\ell - n} \; ,
```

which gives three regimes:

- ``\ell > n`` : the integral vanishes as ``s^{\ell-n}``;
- ``\ell = n`` : the integral tends to the finite, non-zero value ``\sigma_0/(2\ell+1)!!``;
- ``\ell < n`` : the integral diverges as ``s^{-(n-\ell)}``.

Explicitly, for the eight ``I_\ell^n`` used in the code:

```math
\begin{split}
    &I_0^0 = \sigma_0 - \frac{\sigma_{-2}}{6}s^2 + O(s^4) \; , \qquad
     I_1^1 = \frac{\sigma_0}{3} + O(s^2) \; , \qquad
     I_2^2 = \frac{\sigma_0}{15} + O(s^2) \; , \\
    &I_2^0 = \frac{\sigma_{-2}}{15}s^2 + O(s^4) \; , \qquad
     I_4^0 = \frac{\sigma_{-4}}{945}s^4 + O(s^6) \; , \qquad
     I_3^1 = \frac{\sigma_{-2}}{105}s^2 + O(s^4) \; , \\
    &I_0^2 = \frac{\sigma_2}{s^2} - \frac{\sigma_0}{6} + O(s^2) \; , \qquad
     I_1^3 = \frac{\sigma_2}{3\,s^2} - \frac{\sigma_0}{30} + O(s^2) \; .
\end{split}
```

Note that the ``\sigma_i`` with ``i<0`` never survive in the final results: they always appear
multiplied by a positive power of ``\Delta\chi``.

## The small-argument expansion of ``\tilde{I}_0^4``

The `I04_tilde` integral, computed by `func_I04_tilde`, is

```math
    \tilde{I}_0^4(s) = \frac{1}{s^4}\int_0^{+\infty} \frac{\mathrm{d}q}{2\pi^2} \,
        \frac{P(q)}{q^2} \, \left[ j_0(qs) - 1 \right] \; .
```

Using ``j_0(x) - 1 = \sum_{k \geq 1} (-1)^k x^{2k}/(2k+1)!`` we get

```math
    \tilde{I}_0^4(s) = \sum_{k=1}^{+\infty} \frac{(-1)^k \, \sigma_{4-2k}}{(2k+1)!} \,
        s^{\,2k-4}
    = -\frac{\sigma_2}{6\,s^2} + \frac{\sigma_0}{120} - \frac{\sigma_{-2}}{5040}s^2 + O(s^4) \; .
```

It diverges as ``s^{-2}``, and it always appears in the code multiplied by ``\Delta\chi^4``, so
that ``\Delta\chi^4 \, \tilde{I}_0^4(\Delta\chi) \rightarrow 0``.

## How the limit is taken

Since ``\Delta\chi \rightarrow 0`` forces both ``\chi_2 \rightarrow \chi_1`` and ``y \rightarrow 1``,
the limit is a joint one and must be checked to be independent of the direction of approach.
We therefore parametrise the approach with a single parameter ``p``,

```math
    \chi_2 = \chi_1 + p \, \Delta\chi \; , \qquad
    y = \frac{\chi_1^2 + \chi_2^2 - \Delta\chi^2}{2 \, \chi_1 \chi_2} \; ,
    \qquad |p| \leq 1 \; ,
```

where the second relation is just the definition of ``\Delta\chi`` solved for ``y``, and the bound
on ``p`` follows from ``|\chi_1 - \chi_2| \leq \Delta\chi``. We then expand each
``J^{(k)} I_{\ell_k}^{n_k}`` in powers of ``\Delta\chi`` and keep the ``\Delta\chi^0`` coefficient.
**A result that still depends on ``p`` would mean that the limit does not exist**; in all the
cases below the ``p``-dependence cancels in the sum, which is a strong consistency check.

## Family 1: Lensing ``\times`` Lensing

Concerned functions: `integrand_ξ_GNC_Lensing`, `integrand_ξ_LD_Lensing`,
`integrand_ξ_GNCxLD_Lensing_Lensing`. Here ``\Delta\chi^2 = \chi_1^2 + \chi_2^2 - 2\chi_1\chi_2 y``
and the four coefficients are

```math
\begin{split}
    J_{00} &= -\frac{3}{4}\frac{\chi_1^2\chi_2^2}{\Delta\chi^4}(y^2-1)
        \left[8y(\chi_1^2+\chi_2^2) - \chi_1\chi_2(9y^2+7)\right] \; , \\
    J_{02} &= -\frac{3}{2}\frac{\chi_1^2\chi_2^2}{\Delta\chi^4}(y^2-1)
        \left[4y(\chi_1^2+\chi_2^2) - \chi_1\chi_2(3y^2+5)\right] \; , \\
    J_{31} &= 9 \, y \, \Delta\chi^2 \; , \\
    J_{22} &= \frac{9}{4}\frac{\chi_1\chi_2}{\Delta\chi^4}
        \left[2(\chi_1^4+\chi_2^4)(7y^2-3) - 16y\chi_1\chi_2(y^2+1)(\chi_1^2+\chi_2^2)
        + \chi_1^2\chi_2^2(11y^4+14y^2+23)\right] \; ,
\end{split}
```

and the sum is ``J_{00}I_0^0 + J_{02}I_2^0 + J_{31}I_1^3 + J_{22}I_2^2``.

All three square brackets vanish at the singular point. Setting ``y=1`` and
``\chi_1 = \chi_2 = \chi``:

```math
\begin{split}
    &8(2\chi^2) - \chi^2 (9+7) = 16\chi^2 - 16\chi^2 = 0 \; , \\
    &4(2\chi^2) - \chi^2 (3+5) = 8\chi^2 - 8\chi^2 = 0 \; , \\
    &2(2\chi^4)(4) - 16\chi^2(2)(2\chi^2) + \chi^4(11+14+23)
        = 16\chi^4 - 64\chi^4 + 48\chi^4 = 0 \; ,
\end{split}
```

so the leading order is not enough and the expansion must be pushed one order further.
Doing so, the individual contributions are

| term | limit |
|:--|:--|
| ``J_{00} I_0^0`` | ``\dfrac{3}{4}\chi_1^2\sigma_0 \left(-7p^4 + 6p^2 + 1\right)`` |
| ``J_{02} I_2^0`` | ``0`` |
| ``J_{31} I_1^3`` | ``3\,\sigma_2`` |
| ``J_{22} I_2^2`` | ``\dfrac{3}{20}\chi_1^2\sigma_0 \left(35p^4 - 30p^2 + 3\right)`` |

The two ``p``-dependent pieces cancel exactly:

```math
    \frac{-7p^4+6p^2+1}{4} + \frac{35p^4-30p^2+3}{20}
    = \frac{-35p^4+30p^2+5 + 35p^4-30p^2+3}{20} = \frac{8}{20} = \frac{2}{5} \; ,
```

leaving

```math
    \boxed{\;
    \lim_{\Delta\chi \rightarrow 0}
    \left(J_{00}I_0^0 + J_{02}I_2^0 + J_{31}I_1^3 + J_{22}I_2^2\right)
    = 3\,\sigma_2 + \frac{6}{5}\,\chi_1^2\,\sigma_0 \; . }
```

## Family 2: Lensing ``\times`` Doppler

Concerned functions: `integrand_ξ_GNC_Lensing_Doppler`,
`integrand_ξ_GNCxLD_Lensing_Doppler`, `integrand_ξ_GNCxLD_Doppler_Lensing`,
`integrand_ξ_LD_Lensing_Doppler`. The sum is
``J_{00}I_0^0 + J_{02}I_2^0 + J_{04}I_4^0 + J_{20}I_0^2`` with ``J_{20} = y\,\Delta\chi^2``.

| term | limit | why |
|:--|:--|:--|
| ``J_{00} I_0^0`` | ``0`` | ``J_{00}`` is regular and vanishes at ``y=1,\;\chi_1 = s_2`` |
| ``J_{02} I_2^0`` | ``0`` | ``J_{02} \sim \Delta\chi^{-2}`` against ``I_2^0 \sim \Delta\chi^{2}``, numerator vanishes |
| ``J_{04} I_4^0`` | ``0`` | ``J_{04} \sim \Delta\chi^{-2}`` against ``I_4^0 \sim \Delta\chi^{4}`` |
| ``J_{20} I_0^2`` | ``\sigma_2`` | ``y\,\Delta\chi^2 \left(\sigma_2/\Delta\chi^2 - \sigma_0/6 + \dots\right)`` |

```math
    \boxed{\;
    \lim_{\Delta\chi \rightarrow 0}
    \left(J_{00}I_0^0 + J_{02}I_2^0 + J_{04}I_4^0 + J_{20}I_0^2\right) = \sigma_2 \; . }
```

## Family 3: Newtonian ``\times`` Lensing

Concerned functions: `integrand_ξ_GNC_Newtonian_Lensing`,
`integrand_ξ_GNCxLD_Newtonian_Lensing`. The sum is
``J_{00}I_0^0 + J_{02}I_2^0 + J_{04}I_4^0``, with

```math
    J_{00} = \frac{1}{5}\left[ f \chi_2 (3y^2-1) - 3 y s_1 f - 5 y s_1 b \right] \; .
```

Only the first term survives. At ``y=1``, ``\chi_2 = s_1`` we have
``J_{00} \rightarrow \frac{1}{5}\left(2 f s_1 - 3 f s_1 - 5 b s_1\right) = -\frac{s_1(f+5b)}{5}``,
while ``I_0^0 \rightarrow \sigma_0``. For the other two, the numerators vanish at the singular
point:

```math
\begin{split}
    &\text{(} J_{02} \text{, } b \text{ part)} \quad
        -2\chi_2^2 y + \chi_2 s_1 (y^2+3) - 2 y s_1^2 \; \rightarrow \;
        (-2 + 4 - 2)\,s_1^2 = 0 \; , \\
    &\text{(} J_{02} \text{, } f \text{ part)} \quad
        4\chi_2^3(3y^2-1) - 2\chi_2^2 y s_1(3y^2+8) + \chi_2 s_1^2 (9y^2+11) - 6 y s_1^3
        \; \rightarrow \; (8 - 22 + 20 - 6)\,s_1^3 = 0 \; , \\
    &\text{(} J_{04} \text{)} \quad
        \chi_2^5(6y^2-2) + 6\chi_2^4 y s_1(y^2-3) - \chi_2^3 s_1^2 (y^4 + 12y^2 - 21)
        + 2\chi_2^2 y s_1^3 (y^2+3) - 12 \chi_2 s_1^4 + 4 y s_1^5 \\
    &\phantom{\text{(} J_{04} \text{)} \quad}
        \rightarrow \; (4 - 12 + 8 + 8 - 12 + 4)\,s_1^5 = 0 \; ,
\end{split}
```

so that

```math
    \boxed{\;
    \lim_{\Delta\chi \rightarrow 0}
    \left(J_{00}I_0^0 + J_{02}I_2^0 + J_{04}I_4^0\right)
    = -\frac{s_1 \, (f + 5b) \, \sigma_0}{5} \; . }
```

## Family 4: Lensing ``\times`` Local GP

Concerned functions: `integrand_ξ_GNC_Lensing_LocalGP`,
`integrand_ξ_GNCxLD_LocalGP_Lensing`. The sum is

```math
    F \left(\frac{I_0^0}{60} + \frac{I_2^0}{42} + \frac{I_4^0}{140}\right)
    + \frac{y\,\Delta\chi^2}{2} \, I_0^2 \; ,
    \qquad F = 2y\chi_1^2 - \chi_1 s_2 (y^2+3) + 2 y s_2^2 \; .
```

``F`` vanishes at the singular point, ``F \rightarrow (2 - 4 + 2)\chi^2 = 0``, so the first three
terms give zero, and only the last one survives:

```math
    \boxed{\;
    \lim_{\Delta\chi \rightarrow 0} \left[ \; \cdots \; \right] = \frac{\sigma_2}{2} \; . }
```

## Family 5: Newtonian ``\times`` Integrated GP

Concerned functions: `integrand_ξ_GNC_Newtonian_IntegratedGP`,
`integrand_ξ_GNCxLD_Newtonian_IntegratedGP`. The sum is

```math
    F \left(\frac{I_0^0}{15} + \frac{2\,I_2^0}{21} + \frac{I_4^0}{35}\right)
    - \Delta\chi^2 \, (3b + f) \, I_0^2 \; ,
    \qquad F = f \left[(3y^2-1)\chi_2^2 - 4 y s_1 \chi_2 + 2 s_1^2\right] \; .
```

Again ``F \rightarrow f\,(2 - 4 + 2)\,s_1^2 = 0`` and only the ``I_0^2`` term survives:

```math
    \boxed{\;
    \lim_{\Delta\chi \rightarrow 0} \left[ \; \cdots \; \right]
    = -(3b + f) \, \sigma_2 \; . }
```

## Family 6: the ``\Delta\chi^4 \, \tilde{I}_0^4`` terms

Concerned functions: `integrand_ξ_GNC_IntegratedGP`, `integrand_ξ_LD_IntegratedGP`,
`integrand_ξ_GNCxLD_IntegratedGP_IntegratedGP`, `integrand_ξ_GNC_LocalGP_IntegratedGP`,
`integrand_ξ_GNCxLD_IntegratedGP_LocalGP`, `integrand_ξ_GNCxLD_LocalGP_IntegratedGP`,
`integrand_ξ_LD_LocalGP_IntegratedGP`.

All of them contain ``\Delta\chi^4`` multiplying ``\tilde{I}_0^4(\Delta\chi)``, and from the
expansion given above

```math
    \Delta\chi^4 \, \tilde{I}_0^4(\Delta\chi)
    = -\frac{\sigma_2}{6}\Delta\chi^2 + \frac{\sigma_0}{120}\Delta\chi^4 + O(\Delta\chi^6)
    \; \xrightarrow[\Delta\chi \rightarrow 0]{} \; 0 \; .
```

```math
    \boxed{\; \lim_{\Delta\chi \rightarrow 0} \left[ \; \cdots \; \right] = 0 \; . }
```

## Family 7: the ``\Delta\chi^2 \times (\text{vanishing factor})`` terms

Concerned functions: `integrand_ξ_GNC_Doppler_IntegratedGP`,
`integrand_ξ_GNCxLD_IntegratedGP_Doppler`, `integrand_ξ_GNCxLD_Doppler_IntegratedGP`,
`integrand_ξ_LD_Doppler_IntegratedGP`.

These carry an overall ``\Delta\chi^2`` together with a geometric factor of the form
``(\chi \, y - s)`` or ``(s - \chi \, y)``. The ``\Delta\chi^2`` cancels the ``\Delta\chi^{-2}`` of
``I_0^2`` or ``I_1^3``, but the geometric factor vanishes at the singular point, since
``y \rightarrow 1`` and ``\chi \rightarrow s`` give ``\chi y - s \rightarrow 0``. Hence

```math
    \boxed{\; \lim_{\Delta\chi \rightarrow 0} \left[ \; \cdots \; \right] = 0 \; . }
```

## Family 8: the ``J_{22} I_2^2 + J_{31} I_1^3`` terms

Concerned functions: `integrand_ξ_GNC_Lensing_IntegratedGP`,
`integrand_ξ_GNCxLD_IntegratedGP_Lensing`, `integrand_ξ_GNCxLD_Lensing_IntegratedGP`,
`integrand_ξ_GNCxLD_Lensing_LocalGP`, `integrand_ξ_LD_Lensing_IntegratedGP`,
`integrand_ξ_LD_Lensing_LocalGP`.

They all have the structure

```math
    J_{22} \, I_2^2 + J_{31} \, I_1^3 \; , \qquad
    J_{22} = \frac{A}{2}\,\chi_a \chi_b \, (y^2-1) \; , \qquad
    J_{31} = A \, y \, \Delta\chi^2 \; ,
```

with ``A = 1`` for `GNC_Lensing_IntegratedGP`, ``A = 2`` for `GNCxLD_IntegratedGP_Lensing` and
`GNCxLD_Lensing_IntegratedGP`, and ``A = -2`` for the remaining three (where the coefficients are
written as ``-2y\Delta\chi^2`` and ``\chi_a\chi_b(1-y^2)``).

Since ``I_2^2 \rightarrow \sigma_0/15`` is finite while
``(y^2-1) = -(1-y)(1+y) \rightarrow -\Delta\chi^2/\chi_a\chi_b`` vanishes, the ``J_{22}`` term
gives zero; the ``J_{31}`` one gives ``A \, y \Delta\chi^2 \cdot \sigma_2/(3\Delta\chi^2)``:

```math
    \boxed{\; \lim_{\Delta\chi \rightarrow 0}
    \left(J_{22} I_2^2 + J_{31} I_1^3\right) = \frac{A}{3}\,\sigma_2 \; . }
```

## A second way to reach ``\Delta\chi = 0``: the small-``\chi`` corner

!!! warning "This is a known bug, not yet fixed in the code"
    The `Δχ < Δχ_min` branches implemented in the integrands use the limits of the previous
    sections, which are derived under the assumption that ``y \rightarrow 1``. That assumption
    fails in the corner described here, and the seven integrands listed at the end of this
    section return a wrong value there. The derivation below gives the correct expression and
    the one-line change that fixes it.

The statement "``\Delta\chi^2 = 0`` if and only if ``y = 1`` and ``\chi_1 = \chi_2``" is true for
*fixed, non-zero* comoving distances. It is not true uniformly: writing

```math
    \Delta\chi^2 = (\chi_1-\chi_2)^2 + 2\,\chi_1\chi_2\,(1-y) \; ,
```

both terms also vanish when ``\chi_1`` and ``\chi_2`` go to zero **together, at any fixed ``y``**.
This corner is not academic. Every double-``\chi`` TPCF builds its grid as

```julia
χ1s = P1.comdist .* range(1e-6, 1, length = N_χs_2)
χ2s = P2.comdist .* range(1e-6, 1, length = N_χs_2)
```

so the very first node sits at ``\chi \simeq 10^{-6} s \simeq 4 \cdot 10^{-4}\,h^{-1}\mathrm{Mpc}``,
and the whole first row and first column of the grid have
``\Delta\chi \ll \Delta\chi_\mathrm{min} = 10^{-1}`` **for every** ``y``, with full trapezoidal
weight.

### The corner limit

Parametrise the corner as ``\chi_1 = a\,\epsilon``, ``\chi_2 = b\,\epsilon`` with ``y`` fixed, so
that ``\Delta\chi = c\,\epsilon`` with ``c = \sqrt{a^2+b^2-2ab\,y}``, and let
``\epsilon \rightarrow 0``. For family 1, every ``J^{(k)}`` except ``J_{31}`` carries a positive
power of ``\epsilon`` once the ``\epsilon^{-4}`` of ``\Delta\chi^4`` is accounted for, so only
``J_{31} I_1^3 = 9y\,\Delta\chi^2 \cdot \sigma_2/(3\Delta\chi^2)`` survives:

```math
    \boxed{\;
    \lim_{\epsilon \rightarrow 0}
    \left(J_{00}I_0^0 + J_{02}I_2^0 + J_{31}I_1^3 + J_{22}I_2^2\right)
    = 3\,y\,\sigma_2 \; , }
```

independent of ``a`` and ``b``, as it must be. The same argument on family 8, where
``J_{22} = \frac{A}{2}\chi_a\chi_b(y^2-1) = O(\epsilon^2)`` against a finite ``I_2^2``, gives

```math
    \boxed{\;
    \lim_{\epsilon \rightarrow 0}\left(J_{22} I_2^2 + J_{31} I_1^3\right)
    = \frac{A}{3}\,y\,\sigma_2 \; . }
```

Families 2, 3, 4, 5 and 7 correlate one integration variable ``\chi`` against a *fixed* comoving
distance ``s_1`` or ``s_2``, so ``\Delta\chi \rightarrow 0`` still forces ``\chi \rightarrow s > 0``
and ``y \rightarrow 1``: they have no corner. Family 6 is double-``\chi`` but its limit is zero in
both regimes, since ``\Delta\chi^4\tilde{I}_0^4 = -\sigma_2\Delta\chi^2/6 + O(\Delta\chi^4)``
vanishes however ``\Delta\chi`` is made small.

### The fix

The two limits are the two iterated limits of the same function, and a single expression
covers both, because ``3y\sigma_2 \rightarrow 3\sigma_2`` as ``y \rightarrow 1`` while
``\frac{6}{5}\chi_1^2\sigma_0 \rightarrow 0`` as ``\chi_1 \rightarrow 0``:

| family | currently in the code | correct in both regimes |
|:-:|:--|:--|
| 1 | `3 * σ_2 + 6/5 * χ1^2 * σ_0` | `3 * y * σ_2 + 6/5 * χ1^2 * σ_0` |
| 8 | `A/3 * σ_2` | `A/3 * y * σ_2` |

Multiplying the ``\sigma_2`` coefficient by ``y`` is exact to leading order in both limits and
costs nothing. The seven integrands that need it are the double-``\chi`` ones:

| integrand | family |
|:--|:-:|
| `integrand_ξ_GNC_Lensing` | 1 |
| `integrand_ξ_LD_Lensing` | 1 |
| `integrand_ξ_GNCxLD_Lensing_Lensing` | 1 |
| `integrand_ξ_GNC_Lensing_IntegratedGP` | 8 |
| `integrand_ξ_GNCxLD_IntegratedGP_Lensing` | 8 |
| `integrand_ξ_GNCxLD_Lensing_IntegratedGP` | 8 |
| `integrand_ξ_LD_Lensing_IntegratedGP` | 8 |

(`integrand_ξ_GNCxLD_Lensing_LocalGP` and `integrand_ξ_LD_Lensing_LocalGP` are family 8 but
single-``\chi``, so they are not affected; giving them the `y` anyway keeps the six of the
family uniform and changes nothing.)

### Measured effect

`integrand_ξ_GNCxLD_Lensing_Lensing` at ``s_1 = 435.37``, ``s_2 = 1000``, ``y = 0.7``, with
``\chi_1 = 0.9\,\epsilon`` and ``\chi_2 = 1.1\,\epsilon``:

| ``\epsilon`` | ``J\,I`` sum | current limit branch | ratio |
|--:|--:|--:|--:|
| ``10^{-1}`` | ``4.7167\cdot10^{-13}`` | ``6.7404\cdot10^{-13}`` | ``0.6998`` |
| ``10^{-2}`` | ``4.7200\cdot10^{-13}`` | ``6.7388\cdot10^{-13}`` | ``0.7004`` |
| ``10^{-3}`` | ``4.7214\cdot10^{-13}`` | ``6.7389\cdot10^{-13}`` | ``0.7006`` |
| ``10^{-4}`` | ``4.7227\cdot10^{-13}`` | ``6.7389\cdot10^{-13}`` | ``0.7008`` |

The ratio is ``y``, exactly as predicted, and the ``J\,I`` sum is perfectly well conditioned here:
``\Delta\chi/\chi = O(1)`` in the corner, so the bracket cancellation that ruins the sum near
the singular configuration simply does not occur. The current branch is therefore replacing a
good value by one that is a factor ``1/y`` too large — and, for ``y < 0``, of the wrong sign.

Integrated up, this moves `ξ_GNCxLD_Lensing_Lensing` by ``2\%`` at ``\mu = 0.5`` and ``4\%`` at
``s = 10``, ``\mu = 1``, which is what makes 8 assertions of
`test_GNCxLD_SumXiMultipoles_P1.jl` fail against reference data generated before the limit
branches existed.

### A better guard

The clean way to keep both regimes apart is to make the threshold **relative to the local
comoving distances** rather than an absolute length,

```julia
Δχ < Δχ_min * max(χ1, χ2)     # Δχ_min ≃ 1e-4 ≃ eps()^(1/4)
```

which is what the roundoff analysis of the previous section asks for anyway
(``\Delta\chi_\mathrm{break} \sim \chi\,\varepsilon^{1/4}``). In the corner
``\Delta\chi/\chi = O(1)``, so the branch is simply not taken and the well-conditioned ``J\,I``
sum is used; near the singular configuration ``\Delta\chi/\chi \rightarrow 0`` and it is. This
is the role the commented-out `func_Δχ_min` was meant to play, except that it scales with the
separation ``s``, which stays of order hundreds while ``\chi \rightarrow 0`` and so does not
separate the two cases.

## Summary

| integrand | family | limit of the ``J\,I`` sum |
|:--|:-:|:--|
| `integrand_ξ_GNC_Lensing` | 1 | ``3\sigma_2 + \frac{6}{5}\chi_1^2\sigma_0`` |
| `integrand_ξ_LD_Lensing` | 1 | ``3\sigma_2 + \frac{6}{5}\chi_1^2\sigma_0`` |
| `integrand_ξ_GNCxLD_Lensing_Lensing` | 1 | ``3\sigma_2 + \frac{6}{5}\chi_1^2\sigma_0`` |
| `integrand_ξ_GNC_Lensing_Doppler` | 2 | ``\sigma_2`` |
| `integrand_ξ_GNCxLD_Lensing_Doppler` | 2 | ``\sigma_2`` |
| `integrand_ξ_GNCxLD_Doppler_Lensing` | 2 | ``\sigma_2`` |
| `integrand_ξ_LD_Lensing_Doppler` | 2 | ``\sigma_2`` |
| `integrand_ξ_GNC_Newtonian_Lensing` | 3 | ``-\frac{1}{5}s_1(f+5b)\sigma_0`` |
| `integrand_ξ_GNCxLD_Newtonian_Lensing` | 3 | ``-\frac{1}{5}s_1(f+5b)\sigma_0`` |
| `integrand_ξ_GNC_Lensing_LocalGP` | 4 | ``\frac{1}{2}\sigma_2`` |
| `integrand_ξ_GNCxLD_LocalGP_Lensing` | 4 | ``\frac{1}{2}\sigma_2`` |
| `integrand_ξ_GNC_Newtonian_IntegratedGP` | 5 | ``-(3b+f)\sigma_2`` |
| `integrand_ξ_GNCxLD_Newtonian_IntegratedGP` | 5 | ``-(3b+f)\sigma_2`` |
| `integrand_ξ_GNC_IntegratedGP` | 6 | ``0`` |
| `integrand_ξ_LD_IntegratedGP` | 6 | ``0`` |
| `integrand_ξ_GNCxLD_IntegratedGP_IntegratedGP` | 6 | ``0`` |
| `integrand_ξ_GNC_LocalGP_IntegratedGP` | 6 | ``0`` |
| `integrand_ξ_GNCxLD_IntegratedGP_LocalGP` | 6 | ``0`` |
| `integrand_ξ_GNCxLD_LocalGP_IntegratedGP` | 6 | ``0`` |
| `integrand_ξ_LD_LocalGP_IntegratedGP` | 6 | ``0`` |
| `integrand_ξ_GNC_Doppler_IntegratedGP` | 7 | ``0`` |
| `integrand_ξ_GNCxLD_IntegratedGP_Doppler` | 7 | ``0`` |
| `integrand_ξ_GNCxLD_Doppler_IntegratedGP` | 7 | ``0`` |
| `integrand_ξ_LD_Doppler_IntegratedGP` | 7 | ``0`` |
| `integrand_ξ_GNC_Lensing_IntegratedGP` | 8 | ``\frac{1}{3}\sigma_2`` |
| `integrand_ξ_GNCxLD_IntegratedGP_Lensing` | 8 | ``\frac{2}{3}\sigma_2`` |
| `integrand_ξ_GNCxLD_Lensing_IntegratedGP` | 8 | ``\frac{2}{3}\sigma_2`` |
| `integrand_ξ_GNCxLD_Lensing_LocalGP` | 8 | ``-\frac{2}{3}\sigma_2`` |
| `integrand_ξ_LD_Lensing_IntegratedGP` | 8 | ``-\frac{2}{3}\sigma_2`` |
| `integrand_ξ_LD_Lensing_LocalGP` | 8 | ``-\frac{2}{3}\sigma_2`` |

In every case the limit is multiplied by the same prefactor (`common`, `factor`, `denomin`,
...) that multiplies the ``J\,I`` sum in the `Δχ ≥ Δχ_min` branch, so only the bracket needs to
be replaced.

## On the value of `Δχ_min`

The threshold is squeezed between two errors that grow in opposite directions, so it cannot
be made arbitrarily small.

**From above**, the expansions of the previous sections are controlled by the dimensionless
combination ``q \, \Delta\chi``, and the ``q`` integration runs up to ``k_\mathrm{max}``, which is
``10 \, h \, \mathrm{Mpc}^{-1}`` by default. The leading term is accurate only for
``\Delta\chi \ll 1/k_\mathrm{max} = 0.1 \, h^{-1}\mathrm{Mpc}``, so the relative error of the
limit is ``O\left[(k_\mathrm{max}\Delta\chi)^2\right]``: at the default
`Δχ_min = 1e-1` it is of order unity.

**From below**, the ``J\,I`` sum becomes unusable, but not for the reason one might expect: the
four products ``J^{(k)} I_{\ell_k}^{n_k}`` do *not* nearly cancel against each other. The loss
of significance happens one level down, **inside each ``J^{(k)}``**. Take ``J_{22}`` of family 1:
its square bracket is a sum of terms of size ``O(\chi^6)`` whose value at the singular point is
zero (that is exactly what was shown in the derivation), so for small ``\Delta\chi`` the bracket
is ``O(\chi^4\Delta\chi^2)`` — a relative cancellation of ``(\Delta\chi/\chi)^2`` — and the result
is then divided by ``\Delta\chi^4``. The absolute rounding error of the bracket,
``\varepsilon\,\chi^6``, therefore reaches ``J_{22}`` multiplied by ``\chi/\Delta\chi^4``, while
``J_{22}`` itself is ``O(\chi^3)``:

```math
    \frac{\delta J_{22}}{J_{22}} \; \sim \; \varepsilon
        \left(\frac{\chi}{\Delta\chi}\right)^{4} \; ,
    \qquad\text{so}\qquad
    \Delta\chi_\mathrm{break} \; \sim \; \chi \, \varepsilon^{1/4}
        \; \simeq \; 1.2 \cdot 10^{-4} \, \chi \; .
```

The same argument applies to ``J_{00}`` and ``J_{02}``, whose brackets vanish too. Evaluating the
four terms of `integrand_ξ_GNC_Lensing` separately at ``\chi_1 = 250``,
``\chi_2 = \chi_1 + 0.4\,\Delta\chi`` (all values in units of the enhancer):

| ``\Delta\chi`` | ``J_{00}I_0^0`` | ``J_{02}I_2^0`` | ``J_{31}I_1^3`` | ``J_{22}I_2^2`` | sum | limit |
|--:|--:|--:|--:|--:|--:|--:|
| ``1`` | ``4.28\cdot10^{5}`` | ``-2.37\cdot10^{4}`` | ``3.00\cdot10^{2}`` | ``-6.44\cdot10^{4}`` | ``3.40\cdot10^{5}`` | ``1.17\cdot10^{6}`` |
| ``3\cdot10^{-1}`` | ``9.28\cdot10^{5}`` | ``-3.64\cdot10^{4}`` | ``3.02\cdot10^{2}`` | ``-1.26\cdot10^{5}`` | ``7.65\cdot10^{5}`` | ``1.17\cdot10^{6}`` |
| ``1\cdot10^{-1}`` | ``1.61\cdot10^{6}`` | ``-5.33\cdot10^{4}`` | ``3.03\cdot10^{2}`` | ``-7.12\cdot10^{4}`` | ``1.49\cdot10^{6}`` | ``1.17\cdot10^{6}`` |
| ``5\cdot10^{-2}`` | ``2.20\cdot10^{6}`` | ``-6.78\cdot10^{4}`` | ``3.03\cdot10^{2}`` | ``0`` | ``2.13\cdot10^{6}`` | ``1.17\cdot10^{6}`` |
| ``3\cdot10^{-2}`` | ``2.73\cdot10^{6}`` | ``-8.09\cdot10^{4}`` | ``3.03\cdot10^{2}`` | ``1.45\cdot10^{7}`` | ``1.72\cdot10^{7}`` | ``1.17\cdot10^{6}`` |
| ``1\cdot10^{-2}`` | ``4.26\cdot10^{6}`` | ``-1.18\cdot10^{5}`` | ``3.03\cdot10^{2}`` | ``-1.81\cdot10^{9}`` | ``-1.81\cdot10^{9}`` | ``1.17\cdot10^{6}`` |
| ``1\cdot10^{-3}`` | ``1.02\cdot10^{7}`` | ``-2.61\cdot10^{5}`` | ``3.03\cdot10^{2}`` | ``-4.26\cdot10^{13}`` | ``-4.26\cdot10^{13}`` | ``1.17\cdot10^{6}`` |

``J_{22}I_2^2`` rounds to exactly zero at ``\Delta\chi = 5\cdot10^{-2}`` and is pure noise below
it — right at the predicted ``\chi\,\varepsilon^{1/4} \simeq 3\cdot 10^{-2}`` — after which it
runs away by four orders of magnitude per decade. A single quadrature node landing at
``\Delta\chi = 10^{-3}`` contributes an integrand ``10^{7}`` times too large: enough to destroy
the whole ``\chi`` integral.

The two errors cross between ``\Delta\chi = 10^{-1}`` and ``3\cdot10^{-1}``: the sum still tracks
the true bracket at ``3\cdot10^{-1}`` (where the limit is ``35\%`` off), already overshoots it by
``27\%`` at ``10^{-1}``, and is meaningless below ``5\cdot10^{-2}``. So `Δχ_min = 1e-1` sits close
to the optimum, and **lowering it is not safer, it is dangerous**.

Two consequences worth keeping in mind:

- The breakdown scale is ``\chi\,\varepsilon^{1/4}``, so the right threshold is **proportional
  to the comoving distances**, not an absolute length. A fixed `Δχ_min = 1e-1` is tuned for
  ``\chi \sim`` a few hundred ``h^{-1}\mathrm{Mpc}`` and is too small at large ``\chi``, too large
  at small ``\chi``. The commented-out `func_Δχ_min(s1, s2, y; frac)` in
  `GNC_LensingIntegratedGP.jl` is exactly the relative threshold this argument calls for; it
  wants ``\mathrm{frac} \simeq \varepsilon^{1/4} \simeq 10^{-4}`` against ``s``, and reviving it
  would remove the tuning.
- The residual error of the limit is confined to an interval of length `Δχ_min` out of a
  ``\chi`` range of hundreds of ``h^{-1}\mathrm{Mpc}``, so its effect on the integrated TPCF stays
  at the ``10^{-3}`` level for a generic ``y``; it grows to a few per cent at ``y = 1``, where the
  quadrature deliberately samples the singular point.

!!! warning "`integrand_ξ_LD_Lensing`"
    This is the one function that still uses `Δχ_min = 1e-4`, inherited from before these
    limits were derived. That is three orders of magnitude inside the region where the
    ``J^{(k)}`` are noise, so it should be aligned with the `1e-1` used everywhere else.
