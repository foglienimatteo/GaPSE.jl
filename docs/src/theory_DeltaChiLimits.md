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

This page sets up the problem, fixes the notation and collects **all** the results, which are
the ones used in the `Δχ < Δχ_min` branch of the corresponding `integrand_ξ_...` functions.
The derivations themselves are one per page, grouped into eight *families* of integrands that
share the same singular structure — see [The eight families](#The-eight-families) below.

Each family page is self-contained: it repeats at the top every equation of this one that it
needs, so that it can be read without jumping back and forth.

## Definitions and why the limit is needed

We report here the important definitions:

![Positions of the observer $\mathbf{O}$ and of the galaxies $\vs_1$ and $\vs_2$, together with their separation $\vs = \vs_2 - \vs_2$](assets/sketches/sketch_s1-s2-s-1.png)

```math
\Delta\chi := \sqrt{\chi_1^2 + \chi_2^2 - 2 \, \chi_1 \, \chi_2 \,y} \quad \quad (D.1) \\[10pt]

\Delta\chi_1 := \sqrt{\chi_1^2 + s_2^2 - 2 \, \chi_1 \,s_2 \,y} \quad \quad (D.2) \\[10pt]
\Delta\chi_2 := \sqrt{s_1^2 + \chi_2^2 - 2 \,s_1 \,\ \chi_2 \,y} \quad \quad (D.3)  \\[10pt]

(1.4\mathrm{a}): \quad y = \cos{\theta} := \hat{\mathbf{s}}_1 \cdot \hat{\mathbf{s}}_2 \quad \Rightarrow \quad \; -1 \leq y \leq 1 \\[16pt]
```


```math
\begin{align*}
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    y &= \cos \theta = \hat{\mathbf{s}}_1 \cdot  \hat{\mathbf{s}}_2 
        \quad &(1.4\mathrm{a}) \quad
    && \mu &= \cos \varphi = \hat{\mathbf{s}}_1 \cdot  \hat{\mathbf{s}} 
        \quad &(1.4\mathrm{b}) \\[15pt]
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    s(s_1, s_2, y) &= \sqrt{s_1^2 + s_2^2 - 2 \, s_1 \, s_2 \, y} 
        \quad &(1.5\mathrm{a}) \quad
    && s_2(s_1, s, \mu) &= \sqrt{s_1^2 + s^2 + 2 \, s_1 \, s \, \mu} 
        \quad &(1.5\mathrm{b})\\[15pt]
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    \mu(s_1, s_2, y) &= \frac{y \, s_2 - s_1}{s(s_1, s_2, y)} 
        \quad &(1.6\mathrm{a})\quad
    && y(s_1, s, \mu) &= \frac{\mu \, s + s_1}{s_2(s_1, s, \mu)}
        \quad &(1.6\mathrm{b})
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\end{align*}
```




We will consider only ``\Delta\chi`` in the following steps, but the analysis is the same:

```math
\begin{align*}
\Delta\chi^2 &= \chi_1^2 + \chi_2^2 - 2\,\chi_1\chi_2 \,y  \\[10pt]
&= \chi_1^2 + \chi_2^2 - 2\,\chi_1\chi_2 + 2\,\chi_1\chi_2 - 2\,\chi_1\chi_2 \,y  \\[10pt]
&= (\chi_1-\chi_2)^2 + 2\chi_1\chi_2(1-y)\\[10pt]
&\quad\quad\quad (\chi_1-\chi_2)^2 \geq 0 \;, \quad \forall\chi_1,\chi_2\in \mathbb{R}\\[10pt]
&\quad\quad\quad 2\chi_1\chi_2(1-y)\geq 0 \quad \forall\chi_1,\chi_2\in \mathbb{R}\;, \quad \forall y \in [-1,1]\\[10pt]

\quad\quad \Delta\chi = 0\, &\iff y=1 \; \land \; \chi_1=\chi_2 \quad\quad (1.7\mathrm{a})\\[10pt]
\Rightarrow \quad \Delta\chi_1 = 0 &\iff y=1 \; \land \; \chi_1=s_2 \quad\quad (1.7\mathrm{b})\\[10pt]
\quad\quad \Delta\chi_2 = 0 &\iff y=1 \; \land \; s_1=\chi_2 \quad\quad (1.7\mathrm{c})\\[10pt]
\end{align*}
```

Putting the condition ``y=1`` into Eq. (1.6b) and Eq. (1.5b):

```math
\begin{align*}
(1.6\mathrm{b}) \; \mathrm{and} \; (1.6\mathrm{a})  \; : \quad y = 1 &\iff 1 = y(s_1, s, \mu) = \frac{\mu \, s + s_1}{s_2(s_1, s, \mu)} =  \frac{\mu \, s + s_1}{\sqrt{s_1^2 + s^2 + 2 \, s_1 \, s \, \mu} }\\[10pt]
&\iff s_1^2 + s^2 + 2 \, s_1 \, s \, \mu = \left(\mu \, s + s_1 \right)^2 \\[10pt]
&\iff s_1^2 + s^2 + 2 \, s_1 \, s \, \mu = \mu^2 \, s^2 + s_1^2 + 2 \, s_1 \, s \, \mu \\[10pt]
&\iff s^2 (1 - \mu^2) = 0 \\[10pt]
&\iff \mu = \pm 1 \\[15pt]
\Rightarrow \quad y=1&\iff \mu = \pm 1 \quad \quad (1.8)\\[10pt]
\end{align*}
```

We then understand that ``\Delta\chi = 0`` is not a pathological configuration that can be avoided: it is reached deterministically by the quadrature.

The nodes ``\mu = \pm 1`` belong to the Gauss-Lobatto grid (`alg = :lobatto`) and to the trapezoidal grid (`alg = :trap`); only `alg = :quad` avoids them (which is a problem, because are the parts that contributes the most).

``\chi_1 = \chi_2`` is hit whenever a sampling point of the ``\chi`` grid coincides with the other distance. With `suit_sampling = true` this is not even accidental: `sample_subdivision_middle` deliberately places a dense, symmetric sub-grid around the singular point, because that is where the integrand varies fastest.

## How the limit is taken

We got the limit in the Iln Integral page:

```math
\boxed{
    I_\ell^n(s) \; \underset{s \rightarrow 0}{\sim}  \;
        \frac{\sigma_{n-\ell}}{(2\ell+1)!!} \, s^{\,\ell - n} 
}
\quad \quad (2.1\mathrm{a})
\qquad\qquad
\boxed{
    \tilde{I}_0^4(s) \; \underset{s \rightarrow 0}{\sim} \; - \frac{\sigma_2}{6 \, s^2}
}
\quad \quad (2.1\mathrm{b})
```

So, written explicitly:
```math
\begin{align*}
I_0^0           &\sim \sigma_0                         &&\rightarrow \mathrm{const}^{(*)} \quad\quad\quad
    &I_2^2           &\sim \frac{\sigma_0}{15}              &&\rightarrow \mathrm{const} \\[10pt]
I_2^0           &\sim \frac{\sigma_{-2}}{15} \, s^2    &&\rightarrow 0              
    & I_3^1          &\sim \frac{\sigma_{-2}}{105} \, s^2   &&\rightarrow 0              \\[10pt]
I_4^0           &\sim \frac{\sigma_{-4}}{945} \, s^4   &&\rightarrow 0              
    &I_1^3           &\sim \frac{\sigma_2}{3} \, s^{-2}     &&\rightarrow +\infty        \\[10pt]
I_0^2           &\sim \sigma_2 \, s^{-2}               &&\rightarrow +\infty        
    &I_1^1           &\sim \frac{\sigma_0}{3}               &&\rightarrow \mathrm{const} \\[10pt]
\end{align*}
```

```math
\tilde{I}_0^4   \sim -\frac{\sigma_2}{6}\, s^{-2}     \rightarrow +\infty\\[10pt]
```

``{}^{(*)}`` ``\sigma_0`` depends on the integration range over which it is computed,
unlike ``\sigma_2`` and ``\sigma_3``. See
[What the ``\sigma_i`` depend on](#what-the-sigma_i-depend-on) below: the results do not
change, but the range has to be quoted with them.

Since ``\Delta\chi \rightarrow 0`` forces both ``\chi_2 \rightarrow \chi_1`` and ``y \rightarrow 1``,
the limit is a joint one and must be checked to be independent of the direction of approach.
We therefore parametrise the approach with a single parameter ``p``,

```math
    \chi_2 := \chi_1 + p \, \Delta\chi  \quad \quad (D.4)\\[10pt]
    |\chi_1 - \chi_2| \leq \Delta\chi \quad \Rightarrow \quad |p| \leq 1
```
NOTE: the bound on ``p`` follows from ``|\chi_1 - \chi_2| \leq \Delta\chi``

```math
\begin{align*}
\mathrm{Inverting \; }(D.1)\; : \quad 
    y &= \frac{\chi_1^2 + \chi_2^2 - \Delta\chi^2}{2 \, \chi_1 \chi_2} \\[10pt]
      &= \frac{\chi_1^2 + \chi_1^2 + p^2 \Delta\chi^2 + 2 \chi_1 p \Delta \chi - \Delta\chi^2}{2 \, \chi_1 (\chi_1 + p \Delta \chi)} \\[10pt]
      &= \frac{2\chi_1^2 + 2 \chi_1 p \Delta \chi - (1-p^2)\Delta\chi^2}{2 \, \chi_1 (\chi_1 + p \Delta \chi)} \\[10pt]

\Rightarrow \quad 
    y-1  &= \frac{2\chi_1^2 + 2 \chi_1 p \Delta \chi - (1-p^2)\Delta\chi^2 }{2 \, \chi_1 (\chi_1 + p \Delta \chi)} -1 \\[10pt]
        &= \frac{2\chi_1^2 + 2 \chi_1 p \Delta \chi - (1-p^2)\Delta\chi^2 - 2 \chi_1^2 - 2 \chi_1 p \Delta \chi}{2 \, \chi_1 (\chi_1 + p \Delta \chi)} \\[10pt]
        &= - \frac{(1-p^2)\Delta\chi^2}{2 \, \chi_1 (\chi_1 + p \Delta \chi)} \\[10pt]
        &\underset{\Delta\chi \rightarrow 0^{+}}{\sim} - \frac{(1-p^2)}{2\chi_1^2}\Delta\chi^2 \\[10pt]

\Rightarrow \quad 
    y^2-1  &= (y-1)(y+1)\\[10pt]
        &= - \frac{(1-p^2)\Delta\chi^2}{2 \, \chi_1 (\chi_1 + p \Delta \chi)} 
        \left[ \frac{2\chi_1^2 + 2 \chi_1 p \Delta \chi - (1-p^2)\Delta\chi^2 }{2 \, \chi_1 (\chi_1 + p \Delta \chi)} +1 \right] \\[10pt]
        &= - \frac{(1-p^2)\Delta\chi^2 [4\chi_1^2 + 4 \chi_1 p \Delta \chi + (1+p^2)\Delta\chi^2]}{4 \, \chi_1^2 (\chi_1 + p \Delta \chi)^2} 
        \\[10pt]
        &\underset{\Delta\chi \rightarrow 0^{+}}{\sim}- \frac{(1-p^2)}{\chi_1^2}\Delta\chi^2 \\[15pt]

\end{align*}
```

```math
\Rightarrow \quad 
    y-1  \underset{\Delta\chi\rightarrow 0^{+}}{\sim}  - \frac{(1-p^2)}{2\chi_1^2}\Delta\chi^2
    \quad \quad (2.3\mathrm{a})
    \quad , \quad \quad \quad 
    y^2-1 \underset{\Delta\chi \rightarrow 0^{+}}{\sim}- \frac{(1-p^2)}{\chi_1^2}\Delta\chi^2
    \quad \quad (2.3\mathrm{b}) \\[10pt]
```

We then expand each ``J\, I_{\ell}^{n}`` in powers of ``\Delta\chi`` and keep the ``\Delta\chi^0`` coefficient.
**If the result still depends on ``p``, it means that the limit does not exist**; in all the
eight families the ``p``-dependence cancels in the sum, which is a strong consistency check.

### A trap: the two small parameters are NOT of the same order

This is the single easiest way to get these limits wrong, so it deserves to be spelled out.

Along the path ``\chi_2 = \chi_1 + p\,\Delta\chi`` the two quantities that vanish do so at
**different rates**:

```math
    \chi_2 - \chi_1 = p \, \Delta\chi = \mathcal{O}(\Delta\chi) \; ,
    \qquad \qquad
    y - 1 \; \sim \; -\frac{(1-p^2)}{2\chi_1^2}\Delta\chi^2 = \mathcal{O}(\Delta\chi^2) \; .
    \quad \quad (2.4)
```

Consequently, inside a bracket that vanishes at the singular point:

```math
    \underbrace{(\chi_1-\chi_2)^2}_{\mathcal{O}(\Delta\chi^2)}
    \quad \mathrm{and} \quad
    \underbrace{(y-1)}_{\mathcal{O}(\Delta\chi^2)}
    \quad \mathrm{contribute \; at \; the \; SAME \; order.}
    \quad \quad (2.5)
```

**So one must never set ``\chi_1 = \chi_2`` first and expand in ``y`` afterwards.** Doing so
silently throws away every term quadratic in ``(\chi_1 - \chi_2)``, which is exactly as large
as the ``(y-1)`` terms that are being kept.

A concrete example, the one that matters most below. Take the square bracket of
``J^{\kappa\kappa}_{00}``, i.e. Eq.(D.5) of
[Family 1](theory_DeltaChiLimits_1_LensingLensing.md):

```math
    B_{00} := 8 y (\chi_1^2 + \chi_2^2) - \chi_1\chi_2 (9y^2+7) \; .
```

It is true that ``B_{00} \rightarrow 0`` when ``y \rightarrow 1`` **and**
``\chi_2 \rightarrow \chi_1`` together. But if we only set ``y = 1``, keeping
``\chi_1 \neq \chi_2``, we get

```math
    B_{00}\big|_{y=1} = 8(\chi_1^2+\chi_2^2) - 16\chi_1\chi_2 = 8(\chi_1-\chi_2)^2
    = 8 p^2 \Delta\chi^2 \; \neq \; 0 \; ,
    \quad \quad (2.6)
```

which is ``\mathcal{O}(\Delta\chi^2)``: precisely the order we are trying to extract.
Conversely, setting ``\chi_1 = \chi_2`` first gives
``B_{00}\big|_{\chi_1=\chi_2} = \chi_1^2\,(16y - 9y^2 - 7) = -\chi_1^2 (y-1)(9y-7)``,
which at ``y \rightarrow 1`` is ``2\chi_1^2(1-y)``, i.e. it keeps **only** the
``(y-1)`` half of the answer and loses the ``8p^2\Delta\chi^2`` of Eq.(2.6).

The safe recipe, used systematically in every family page, is:

> Rewrite the bracket **exactly** — no approximation — as a combination of the vanishing
> quantities ``(\chi_1-\chi_2)``, ``(y-1)`` and ``(y^2-1)``, with coefficients that are
> regular at the singular point. Only then substitute the leading orders (D.4), (2.3a),
> (2.3b).




### What the ``\sigma_i`` depend on

The table above, and every boxed result of the eight family pages, is written as if the
``\sigma_i`` were numbers. They are — but numbers *attached to an integration range*, and
that is worth one paragraph because it is easy to get wrong.

The ``\sigma_i`` and the ``I_\ell^n`` are defined over the same finite
``[k_\mathrm{min}, k_\mathrm{max}]`` (see
[The ``I_\ell^n`` integrals](theory_IlnIntegrals.md)). Over that range every ``\sigma_i``
is finite, and Eq.(2.1a) follows rigorously: ``q s \leq k_\mathrm{max} s`` uniformly, so
the small-argument form of ``j_\ell`` may be inserted under the integral. Three of the
five moments, ``\sigma_1``, ``\sigma_2`` and ``\sigma_3``, would remain finite even if the
cuts were removed, so their value does not really depend on the range. The other two do:

- ``\sigma_0`` grows with ``k_\mathrm{max}`` (as ``K^{0.359}`` with the fitted tail of the
  input spectrum, logarithmically with the asymptotic CDM one) — it is the density
  variance smoothed on zero scale;
- ``\sigma_4`` grows as ``k_\mathrm{min}`` is lowered.

So ``\sigma_0`` cannot be quoted without saying over which range it was computed, and the
five branches that use it — the three of Family 1 and the two of Family 3 — inherit that
dependence.

**Where the approximation actually sits.** Eq.(2.1a) is exact only for
``\Delta\chi \ll 1/k_\mathrm{max}``, while the limit branch switches over at
``\Delta\chi < \Delta\chi_\mathrm{min} = 0.1``. With the `xicalc` range
``k_\mathrm{max} = 10^{3}`` the two are three decades apart:
``I_0^0(0.1) = 23.7`` against ``\sigma_0 = 143.3``. This is a real approximation in the
code, not a bookkeeping subtlety, and it is the reason the stored ``\sigma_i`` and the
``I_\ell^n`` should be computed consistently — the discussion, with both options written
out, is in
[What this means for the limits](theory_InputPowerSpectrum.md#what-this-means-for-the-delta-chi-rightarrow-0-limits).

## A pattern worth noticing

Across all eight families the same two polynomials keep appearing, and they are the
Legendre ones evaluated at the direction parameter ``p``:

| where                             | leading coefficient                     |
| :-------------------------------- | :-------------------------------------- |
| ``B_{22}`` , Eq.(3.9) of [Family 1](theory_DeltaChiLimits_1_LensingLensing.md) | ``35p^4-30p^2+3 = 8\,\mathcal{L}_4(p)`` |
| ``N_{04}`` , Eq.(5.4) of [Family 3](theory_DeltaChiLimits_3_NewtonianLensing.md) | ``35p^4-30p^2+3 = 8\,\mathcal{L}_4(p)`` |
| ``N_{02}`` , Eq.(5.3) of [Family 3](theory_DeltaChiLimits_3_NewtonianLensing.md) | ``3p^2-1 = 2\,\mathcal{L}_2(p)``        |
| ``F`` , Eq.(6.2) of [Family 4](theory_DeltaChiLimits_4_LensingLocalGP.md)      | ``3p^2-1 = 2\,\mathcal{L}_2(p)``        |
| ``F`` , Eq.(7.2) of [Family 5](theory_DeltaChiLimits_5_NewtonianIntegratedGP.md)      | ``3p^2-1 = 2\,\mathcal{L}_2(p)``        |

This is not a coincidence. The ``J`` coefficients come from expanding the TPCF kernels in
Legendre polynomials of ``y = \cos\theta``, and along the path ``\chi_2 = \chi_1 + p\Delta\chi``
the angle and the radial separation are locked together by Eq.(D.1) in such a way that
``p`` inherits the role of ``\cos\theta``. It is a useful check when redoing any of these
expansions: a leading coefficient that is *not* a low-order Legendre polynomial in ``p`` is
a good sign that a term has been lost.

## The eight families

The thirty affected integrands fall into eight groups that share the same singular
structure, and therefore the same derivation. One page each:

|     | family                                               | ``J\,I`` structure                                                                  | limit                                       |
| :-: | :--------------------------------------------------- | :---------------------------------------------------------------------------------- | :------------------------------------------ |
|  1  | [Lensing ``\times`` Lensing](theory_DeltaChiLimits_1_LensingLensing.md) | ``J_{00}I_0^0 + J_{02}I_2^0 + J_{31}I_1^3 + J_{22}I_2^2``                           | ``3\sigma_2 + \frac{6}{5}\chi_1^2\sigma_0`` |
|  2  | [Lensing ``\times`` Doppler](theory_DeltaChiLimits_2_LensingDoppler.md) | ``J_{00}I_0^0 + J_{02}I_2^0 + J_{04}I_4^0 + J_{20}I_0^2``                           | ``\sigma_2``                                |
|  3  | [Newtonian ``\times`` Lensing](theory_DeltaChiLimits_3_NewtonianLensing.md) | ``J_{00}I_0^0 + J_{02}I_2^0 + J_{04}I_4^0``                                         | ``-\frac{1}{5}s_1(f+5b)\sigma_0``           |
|  4  | [Lensing ``\times`` Local GP](theory_DeltaChiLimits_4_LensingLocalGP.md) | ``F\left(\frac{I_0^0}{60}+\frac{I_2^0}{42}+\frac{I_4^0}{140}\right) + J_{20}I_0^2`` | ``\frac{1}{2}\sigma_2``                     |
|  5  | [Newtonian ``\times`` Integrated GP](theory_DeltaChiLimits_5_NewtonianIntegratedGP.md) | ``F\left(\frac{I_0^0}{15}+\frac{2I_2^0}{21}+\frac{I_4^0}{35}\right) + J_{20}I_0^2`` | ``-(3b+f)\sigma_2``                         |
|  6  | [the ``\Delta\chi^4\,\tilde{I}_0^4`` terms](theory_DeltaChiLimits_6_Ichi4Tilde.md) | ``\Delta\chi^4 \, \tilde{I}_0^4``                                                   | ``0``                                       |
|  7  | [the ``\Delta\chi^2 \times`` (vanishing factor) terms](theory_DeltaChiLimits_7_VanishingFactor.md) | ``\Delta\chi^2 \, G \left(\cdots + I_0^2\right)``                                   | ``0``                                       |
|  8  | [the ``J_{22}I_2^2 + J_{31}I_1^3`` terms](theory_DeltaChiLimits_8_J22J31.md) | ``J_{22}I_2^2 + J_{31}I_1^3``                                                       | ``\frac{A}{3}\sigma_2``                     |

## Index of the definitions

Every definition carries a ``(\mathrm{D}.n)`` number, so that it can be referred to from
anywhere without being restated. The first four are shared by all the families and are
repeated in the recap of each page; the others are local to the family that uses them.

| | definition | where |
|:-:|:--|:--|
| ``(\mathrm{D}.1)`` | ``\Delta\chi := \sqrt{\chi_1^2+\chi_2^2-2\chi_1\chi_2 y}`` | this page, Families [1](theory_DeltaChiLimits_1_LensingLensing.md), [6](theory_DeltaChiLimits_6_Ichi4Tilde.md), [8](theory_DeltaChiLimits_8_J22J31.md) |
| ``(\mathrm{D}.2)`` | ``\Delta\chi_1 := \sqrt{\chi_1^2+s_2^2-2\chi_1 s_2 y}`` | this page, Families [2](theory_DeltaChiLimits_2_LensingDoppler.md), [4](theory_DeltaChiLimits_4_LensingLocalGP.md) |
| ``(\mathrm{D}.3)`` | ``\Delta\chi_2 := \sqrt{s_1^2+\chi_2^2-2 s_1\chi_2 y}`` | this page, Families [3](theory_DeltaChiLimits_3_NewtonianLensing.md), [5](theory_DeltaChiLimits_5_NewtonianIntegratedGP.md), [7](theory_DeltaChiLimits_7_VanishingFactor.md) |
| ``(\mathrm{D}.4)`` | ``\chi_2 := \chi_1 + p \, \Delta\chi`` (the direction of approach) | this page and every family |
| ``(\mathrm{D}.5)`` | ``B_{00} := 8y(\chi_1^2+\chi_2^2) - \chi_1\chi_2(9y^2+7)`` | [Family 1](theory_DeltaChiLimits_1_LensingLensing.md) |
| ``(\mathrm{D}.6)`` | ``B_{02} := 4y(\chi_1^2+\chi_2^2) - \chi_1\chi_2(3y^2+5)`` | [Family 1](theory_DeltaChiLimits_1_LensingLensing.md) |
| ``(\mathrm{D}.7)`` | ``B_{22} := 2(\chi_1^4+\chi_2^4)(7y^2-3) - \dots`` | [Family 1](theory_DeltaChiLimits_1_LensingLensing.md) |
| ``(\mathrm{D}.8)`` | ``u := \chi_1^2 + \chi_2^2`` , ``v := \chi_1\chi_2`` | [Family 1](theory_DeltaChiLimits_1_LensingLensing.md) |
| ``(\mathrm{D}.9)`` | ``N_{02}`` , the numerator of ``J_{02}`` | [Family 2](theory_DeltaChiLimits_2_LensingDoppler.md) |
| ``(\mathrm{D}.10)`` | ``N_{04}`` , the numerator of ``J_{04}`` | [Family 2](theory_DeltaChiLimits_2_LensingDoppler.md) |
| ``(\mathrm{D}.11)`` | ``N_{02}^{(b)}`` , the bias part of the ``J_{02}`` numerator | [Family 3](theory_DeltaChiLimits_3_NewtonianLensing.md) |
| ``(\mathrm{D}.12)`` | ``N_{02}^{(f)}`` , its growth-rate part | [Family 3](theory_DeltaChiLimits_3_NewtonianLensing.md) |
| ``(\mathrm{D}.13)`` | ``N_{04}`` , the numerator of ``J_{04}`` | [Family 3](theory_DeltaChiLimits_3_NewtonianLensing.md) |
| ``(\mathrm{D}.14)`` | ``F := 2y\chi_1^2 - \chi_1 s_2(y^2+3) + 2y s_2^2`` | [Family 4](theory_DeltaChiLimits_4_LensingLocalGP.md) |
| ``(\mathrm{D}.15)`` | ``J_{20} := y\,\Delta\chi_1^2/2`` | [Family 4](theory_DeltaChiLimits_4_LensingLocalGP.md) |
| ``(\mathrm{D}.16)`` | ``F := f\left[(3y^2-1)\chi_2^2 - 4y s_1\chi_2 + 2 s_1^2\right]`` | [Family 5](theory_DeltaChiLimits_5_NewtonianIntegratedGP.md) |
| ``(\mathrm{D}.17)`` | ``J_{20} := -\Delta\chi_2^2 (3b+f)`` | [Family 5](theory_DeltaChiLimits_5_NewtonianIntegratedGP.md) |
| ``(\mathrm{D}.18)`` | ``G := \chi_2 y - s_1`` | [Family 7](theory_DeltaChiLimits_7_VanishingFactor.md) |
| ``(\mathrm{D}.19)`` | ``J_{22} := \frac{A}{2}\chi_a\chi_b(y^2-1)`` , ``J_{31} := A y \Delta\chi^2`` | [Family 8](theory_DeltaChiLimits_8_J22J31.md) |

!!! note "Equation numbering\"
    Equations are numbered **globally** across these nine pages. Definitions carry the
    ``(\mathrm{D}.n)`` prefix of the index above, wherever they live; everything else is
    numbered by the page it belongs to — ``(1.x)`` and ``(2.x)`` on this page, ``(3.x)`` on
    Family 1, ``(4.x)`` on Family 2, and so on up to ``(10.x)`` on Family 8. A family page
    sometimes quotes an equation derived on another one — in that case the equation is
    reproduced in full where it is used, and its prefix tells where it comes from.


## Summary

| integrand                                      | family | limit of the ``J\,I`` sum                   |
| :--------------------------------------------- | :----: | :------------------------------------------ |
| `integrand_ξ_GNC_Lensing`                      |   [1](theory_DeltaChiLimits_1_LensingLensing.md)    | ``3\sigma_2 + \frac{6}{5}\chi_1^2\sigma_0`` |
| `integrand_ξ_LD_Lensing`                       |   [1](theory_DeltaChiLimits_1_LensingLensing.md)    | ``3\sigma_2 + \frac{6}{5}\chi_1^2\sigma_0`` |
| `integrand_ξ_GNCxLD_Lensing_Lensing`           |   [1](theory_DeltaChiLimits_1_LensingLensing.md)    | ``3\sigma_2 + \frac{6}{5}\chi_1^2\sigma_0`` |
| `integrand_ξ_GNC_Lensing_Doppler`              |   [2](theory_DeltaChiLimits_2_LensingDoppler.md)    | ``\sigma_2``                                |
| `integrand_ξ_GNCxLD_Lensing_Doppler`           |   [2](theory_DeltaChiLimits_2_LensingDoppler.md)    | ``\sigma_2``                                |
| `integrand_ξ_GNCxLD_Doppler_Lensing`           |   [2](theory_DeltaChiLimits_2_LensingDoppler.md)    | ``\sigma_2``                                |
| `integrand_ξ_LD_Lensing_Doppler`               |   [2](theory_DeltaChiLimits_2_LensingDoppler.md)    | ``\sigma_2``                                |
| `integrand_ξ_GNC_Newtonian_Lensing`            |   [3](theory_DeltaChiLimits_3_NewtonianLensing.md)    | ``-\frac{1}{5}s_1(f+5b)\sigma_0``           |
| `integrand_ξ_GNCxLD_Newtonian_Lensing`         |   [3](theory_DeltaChiLimits_3_NewtonianLensing.md)    | ``-\frac{1}{5}s_1(f+5b)\sigma_0``           |
| `integrand_ξ_GNC_Lensing_LocalGP`              |   [4](theory_DeltaChiLimits_4_LensingLocalGP.md)    | ``\frac{1}{2}\sigma_2``                     |
| `integrand_ξ_GNCxLD_LocalGP_Lensing`           |   [4](theory_DeltaChiLimits_4_LensingLocalGP.md)    | ``\frac{1}{2}\sigma_2``                     |
| `integrand_ξ_GNC_Newtonian_IntegratedGP`       |   [5](theory_DeltaChiLimits_5_NewtonianIntegratedGP.md)    | ``-(3b+f)\sigma_2``                         |
| `integrand_ξ_GNCxLD_Newtonian_IntegratedGP`    |   [5](theory_DeltaChiLimits_5_NewtonianIntegratedGP.md)    | ``-(3b+f)\sigma_2``                         |
| `integrand_ξ_GNC_IntegratedGP`                 |   [6](theory_DeltaChiLimits_6_Ichi4Tilde.md)    | ``0``                                       |
| `integrand_ξ_LD_IntegratedGP`                  |   [6](theory_DeltaChiLimits_6_Ichi4Tilde.md)    | ``0``                                       |
| `integrand_ξ_GNCxLD_IntegratedGP_IntegratedGP` |   [6](theory_DeltaChiLimits_6_Ichi4Tilde.md)    | ``0``                                       |
| `integrand_ξ_GNC_LocalGP_IntegratedGP`         |   [6](theory_DeltaChiLimits_6_Ichi4Tilde.md)    | ``0``                                       |
| `integrand_ξ_GNCxLD_IntegratedGP_LocalGP`      |   [6](theory_DeltaChiLimits_6_Ichi4Tilde.md)    | ``0``                                       |
| `integrand_ξ_GNCxLD_LocalGP_IntegratedGP`      |   [6](theory_DeltaChiLimits_6_Ichi4Tilde.md)    | ``0``                                       |
| `integrand_ξ_LD_LocalGP_IntegratedGP`          |   [6](theory_DeltaChiLimits_6_Ichi4Tilde.md)    | ``0``                                       |
| `integrand_ξ_GNC_Doppler_IntegratedGP`         |   [7](theory_DeltaChiLimits_7_VanishingFactor.md)    | ``0``                                       |
| `integrand_ξ_GNCxLD_IntegratedGP_Doppler`      |   [7](theory_DeltaChiLimits_7_VanishingFactor.md)    | ``0``                                       |
| `integrand_ξ_GNCxLD_Doppler_IntegratedGP`      |   [7](theory_DeltaChiLimits_7_VanishingFactor.md)    | ``0``                                       |
| `integrand_ξ_LD_Doppler_IntegratedGP`          |   [7](theory_DeltaChiLimits_7_VanishingFactor.md)    | ``0``                                       |
| `integrand_ξ_GNC_Lensing_IntegratedGP`         |   [8](theory_DeltaChiLimits_8_J22J31.md)    | ``\frac{1}{3}\sigma_2``                     |
| `integrand_ξ_GNCxLD_IntegratedGP_Lensing`      |   [8](theory_DeltaChiLimits_8_J22J31.md)    | ``\frac{2}{3}\sigma_2``                     |
| `integrand_ξ_GNCxLD_Lensing_IntegratedGP`      |   [8](theory_DeltaChiLimits_8_J22J31.md)    | ``\frac{2}{3}\sigma_2``                     |
| `integrand_ξ_GNCxLD_Lensing_LocalGP`           |   [8](theory_DeltaChiLimits_8_J22J31.md)    | ``-\frac{2}{3}\sigma_2``                    |
| `integrand_ξ_LD_Lensing_IntegratedGP`          |   [8](theory_DeltaChiLimits_8_J22J31.md)    | ``-\frac{2}{3}\sigma_2``                    |
| `integrand_ξ_LD_Lensing_LocalGP`               |   [8](theory_DeltaChiLimits_8_J22J31.md)    | ``-\frac{2}{3}\sigma_2``                    |

In every case the limit is multiplied by the same prefactor (`common`, `factor`, `denomin`,
...) that multiplies the ``J\,I`` sum in the `Δχ ≥ Δχ_min` branch, so only the bracket needs to
be replaced.

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
| `integrand_ξ_GNC_Lensing` | [1](theory_DeltaChiLimits_1_LensingLensing.md) |
| `integrand_ξ_LD_Lensing` | [1](theory_DeltaChiLimits_1_LensingLensing.md) |
| `integrand_ξ_GNCxLD_Lensing_Lensing` | [1](theory_DeltaChiLimits_1_LensingLensing.md) |
| `integrand_ξ_GNC_Lensing_IntegratedGP` | [8](theory_DeltaChiLimits_8_J22J31.md) |
| `integrand_ξ_GNCxLD_IntegratedGP_Lensing` | [8](theory_DeltaChiLimits_8_J22J31.md) |
| `integrand_ξ_GNCxLD_Lensing_IntegratedGP` | [8](theory_DeltaChiLimits_8_J22J31.md) |
| `integrand_ξ_LD_Lensing_IntegratedGP` | [8](theory_DeltaChiLimits_8_J22J31.md) |

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
``s = 10``, ``\mu = 1``. In the multipoles it shows up as a **uniform** ``1.5\%`` bias of the
windowed ``\xi_{\kappa\kappa}`` at ``s = 1000``, the same for every ``L`` and every quadrature,
which is what made 24 assertions of `test_GNCxLD_SumXiMultipoles_P1.jl` fail.

### The guard that is implemented

What separates the two regimes is a threshold **relative to the local comoving distances**
rather than an absolute length: in the corner ``\Delta\chi/\chi = O(1)``, so the branch must
not be taken and the well-conditioned ``J\,I`` sum must be used, while near the singular
configuration ``\Delta\chi/\chi \rightarrow 0`` and it must. This is the role the commented-out
`func_Δχ_min` was meant to play, except that it scales with the separation ``s``, which stays
of order hundreds while ``\chi \rightarrow 0`` and so does not separate the two cases.

A *purely* relative threshold is not enough on its own, because the limit has its own upper
validity bound: the next section shows it is accurate only for
``\Delta\chi \ll 1/k_\mathrm{max} = 0.1 \, h^{-1}\mathrm{Mpc}``. Taking
``\Delta\chi < \Delta\chi_\mathrm{min} \, \max(\chi_1,\chi_2)`` while keeping
`Δχ_min = 1e-1` would allow ``\Delta\chi`` up to ``0.1\,\chi``, i.e. up to ``100`` for
``\chi \sim 10^3`` — the branch would fire far outside the regime where the expansion holds.
(Measured, that combination is catastrophic: errors of ``10^3`` to ``3\cdot10^3`` percent.) It
becomes usable only if ``\Delta\chi_\mathrm{min}`` is simultaneously lowered to
``\simeq \varepsilon^{1/4} \simeq 10^{-4}``, which changes the behaviour everywhere, not just
in the corner.

The threshold therefore has to be bounded **both** ways, and what the code implements is the
minimal form that keeps the absolute cap and only tightens it where ``\chi`` is small:

```julia
Δχ ≥ min(Δχ_min, Δχ_min * max(χ1, χ2))     # Δχ_min = 1e-1 unchanged
```

For ``\max(\chi_1,\chi_2) \geq 1`` the `min` selects ``\Delta\chi_\mathrm{min}`` and nothing
changes; below that it scales down with ``\chi`` and the corner is excluded. The same form is
applied to all thirty limit branches, each with its own pair of distances — ``(\chi_1,\chi_2)``
for the double-``\chi`` families, ``(s_1,\chi_2)`` or ``(\chi_1,s_2)`` for the others, where it
is a no-op in practice since the ``s`` are of order hundreds.

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
