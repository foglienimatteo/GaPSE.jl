# Family 4: Lensing ``\times`` Local GP

## Recap: everything this page needs

This page is self-contained. All the results quoted here are derived in
[The ``\Delta\chi \rightarrow 0`` limits](DeltaChiLimits.md), and the ``I_\ell^n``
asymptotics in [The ``I_\ell^n`` integrals](IlnIntegrals.md).

**The separation and its singular point.** For this family the two competing distances are
``\chi_1`` and ``s_2``, so the relevant separation is

```math
    \Delta\chi_1 := \sqrt{\chi_1^2 + s_2^2 - 2 \, \chi_1 \, s_2 \, y} \quad \quad (D.2)
```

```math
    \Delta\chi_1 = 0 \quad \iff \quad y = 1 \; \land \; \chi_1 = s_2 \quad \quad (1.7b)
```

i.e. the limit is a **joint** one: ``y \rightarrow 1`` and ``s_2 \rightarrow \chi_1`` together.

**The small-``s`` behaviour of the ``I_\ell^n``.**

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

written explicitly:

```math
\begin{align*}
I_0^0           &\sim \sigma_0                         &&\rightarrow \mathrm{const} \quad\quad\quad
    &I_2^2           &\sim \frac{\sigma_0}{15}              &&\rightarrow \mathrm{const} \\[10pt]
I_2^0           &\sim \frac{\sigma_{-2}}{15} \, s^2    &&\rightarrow 0              
    & I_3^1          &\sim \frac{\sigma_{-2}}{105} \, s^2   &&\rightarrow 0              \\[10pt]
I_4^0           &\sim \frac{\sigma_{-4}}{945} \, s^4   &&\rightarrow 0              
    &I_1^3           &\sim \frac{\sigma_2}{3} \, s^{-2}     &&\rightarrow +\infty        \\[10pt]
I_0^2           &\sim \sigma_2 \, s^{-2}               &&\rightarrow +\infty        
    &I_1^1           &\sim \frac{\sigma_0}{3}               &&\rightarrow \mathrm{const} \\[10pt]
\tilde{I}_0^4   &\sim -\frac{\sigma_2}{6}\, s^{-2}     &&\rightarrow +\infty
\end{align*}
```

**The parametrisation of the approach.** Since the limit is joint, it must be checked to be
independent of the direction of approach; we parametrise the latter with a single ``p``,

```math
    s_2 := \chi_1 + p \, \Delta\chi_1  \quad \quad (D.4) \\[10pt]
    |\chi_1 - s_2| \leq \Delta\chi_1 \quad \Rightarrow \quad |p| \leq 1
```

**and its two consequences**, obtained by inverting (D.2):

```math
\Rightarrow \quad 
    y-1  \underset{\Delta\chi_1\rightarrow 0^{+}}{\sim}  - \frac{(1-p^2)}{2\chi_1^2}\Delta\chi_1^2
    \quad \quad (2.3\mathrm{a})
    \quad , \quad \quad \quad 
    y^2-1 \underset{\Delta\chi_1 \rightarrow 0^{+}}{\sim}- \frac{(1-p^2)}{\chi_1^2}\Delta\chi_1^2
    \quad \quad (2.3\mathrm{b}) \\[10pt]
```

**The trap to avoid.** The two vanishing quantities have **different orders**,

```math
    s_2 - \chi_1 = p \, \Delta\chi_1 = \mathcal{O}(\Delta\chi_1)
    \qquad \qquad
    y - 1 = \mathcal{O}(\Delta\chi_1^2)
    \quad \quad (2.4)
```

so a term **quadratic** in ``(\chi_1 - s_2)`` contributes at the **same order** as a term
**linear** in ``(y-1)``. One must therefore never set ``\chi_1 = s_2`` first and expand in ``y``
afterwards: the safe recipe is to rewrite every vanishing bracket **exactly** as a
combination of ``(\chi_1-s_2)``, ``(y-1)`` and ``(y^2-1)``, with coefficients that are regular
at the singular point, and only then substitute (D.4), (2.3a), (2.3b).

## The integrand

Concerned functions: `integrand_ξ_GNC_Lensing_LocalGP`,
`integrand_ξ_GNCxLD_LocalGP_Lensing`.
The GNC Lensing-LocalGP function is:

```math
\begin{split}
    \xi^{\kappa \phi} ( s_1 , s_2, y ) &= 
    D_2 \int_0^{s_1}\mathrm{d} \chi_1 \; 
    J^{\kappa \phi}_{\alpha}\left[ 
        J^{\kappa \phi}_{20} I_0^2 ( \Delta \chi_1 ) +
        \right.\nonumber \\
        &\left.\qquad
        J^{\kappa \phi}_{\beta}
        \left(
            \frac{1}{60} I_0^0 ( \Delta \chi_1 ) +
            \frac{1}{42} I_2^0 ( \Delta \chi_1 ) +
            \frac{1}{140} I_4^0 ( \Delta \chi_1 ) 
        \right)
    \right] \, ,
\end{split}
```

with

```math
\begin{align*}
    J^{\kappa \phi}_{\alpha}  &= 
    \frac{\mathcal{H}_0^2 \Omega_{\mathrm{M}0} s_2 D(\chi_1)}{a(\chi_1) a_2 s_1}
    (\chi_1 - s_1)  (5s_{\mathrm{b}, 1} - 2) \times \\
    &\qquad\qquad
    \left[
       2 f_2 a_2 \mathcal{H}_2^2 (\mathit{f}_{\mathrm{evo}, 2} - 3) + 
       3 \mathcal{H}_0^2 \Omega_{\mathrm{M}0} (f_2 + \mathcal{R}_2 + 5s_{\mathrm{b}, 2} - 2)
    \right]
    \, , \\[6pt]
    J^{\kappa \phi}_{\beta} &=
    \underbrace{2 y \chi_1^2 - \chi_1 s_2 (y^2 + 3) + 2 y s_2^2}_{=: \; F \; , \quad (\mathrm{D}.14)}
    \, ,\\[6pt]
    J^{\kappa \phi}_{20} &= \underbrace{\frac{1}{2} y \Delta\chi_1^2}_{(\mathrm{D}.15)} \, .
\end{align*}
```

The limit to be taken is:

```math
    \lim_{\Delta\chi_1\rightarrow 0^{+}}\left[
        F \left(\frac{I_0^0}{60} + \frac{I_2^0}{42} + \frac{I_4^0}{140}\right)
        + J_{20} \, I_0^2 \right] \; .
```

This analysis is valid for Lensing-LocalGP in both the combinations (GNC and GNCxLD).

### Term 1: the ``F`` group

``F``, Eq.(D.14), vanishes at the singular point,

```math
    F\big|_{y=1,\,s_2=\chi_1} = \left(2 - (1+3) + 2\right)\chi_1^2 = (2 - 4 + 2)\,\chi_1^2 = 0 \; ,
```

so it must be decomposed exactly. The structure is the **same** as the ``B_{00}`` of
Family 1, Eq.(D.5) there — which for reference reads

```math
    B_{00} := 8 y (\chi_1^2 + \chi_2^2) - \chi_1\chi_2 (9y^2+7)
    = 8(\chi_1-\chi_2)^2 + 8(y-1)(\chi_1^2+\chi_2^2) - 9\chi_1\chi_2(y^2-1) \; ,
    \quad \quad (3.2)
```

only with different coefficients and with ``s_2`` in place of ``\chi_2``:

```math
\begin{align*}
    F &= 2y\chi_1^2 - \chi_1 s_2 (y^2+3) + 2 y s_2^2 \\[10pt]
    &\quad\quad 2y = 2 + 2(y-1) \\[10pt]
        &= \left[2 + 2(y-1)\right]\chi_1^2 - \chi_1 s_2 (y^2+3) + \left[2 + 2(y-1)\right] s_2^2 \\[10pt]
        &= 2\chi_1^2 + 2 s_2^2 + 2(y-1)\left(\chi_1^2 + s_2^2\right) - \chi_1 s_2 (y^2+3) \\[10pt]
    &\quad\quad y^2+3 = 4 + (y^2-1) \\[10pt]
        &= 2\chi_1^2 + 2 s_2^2 + 2(y-1)\left(\chi_1^2 + s_2^2\right)
            - \chi_1 s_2 \left[4 + (y^2-1)\right] \\[10pt]
        &= 2\chi_1^2 - 4\chi_1 s_2 + 2 s_2^2
            + 2(y-1)\left(\chi_1^2 + s_2^2\right) - (y^2-1)\chi_1 s_2 \\[10pt]
    &\quad\quad 2\chi_1^2 - 4\chi_1 s_2 + 2 s_2^2 = 2(\chi_1 - s_2)^2 \\[10pt]
        &= 2(\chi_1 - s_2)^2 + 2(y-1)\left(\chi_1^2 + s_2^2\right)
            - (y^2-1)\chi_1 s_2 \; .
    \quad \quad (6.1)
\end{align*}
```

Here — exactly as in ``B_{00}``, and unlike Family 2 — the ``(\chi_1 - s_2)`` factor is
**quadratic**, so by (2.4) it is of the same order as the ``(y-1)`` terms and must be kept.
Substituting (D.4), (2.3a), (2.3b) with ``s_2 \rightarrow \chi_1`` in the regular
coefficients:

```math
\begin{align*}
    F &\underset{\Delta\chi_1\rightarrow 0^{+}}{\sim}
        \underbrace{2 (-p\Delta\chi_1)^2}_{(\chi_1-s_2)^2}
        + \underbrace{2\left[- \frac{(1-p^2)}{2\chi_1^2}\Delta\chi_1^2\right]\left(2\chi_1^2\right)}_{(y-1)}
        - \underbrace{\left[- \frac{(1-p^2)}{\chi_1^2}\Delta\chi_1^2\right]\chi_1^2}_{(y^2-1)} \\[10pt]
    &= 2p^2\Delta\chi_1^2 - 2(1-p^2)\Delta\chi_1^2 + (1-p^2)\Delta\chi_1^2 \\[10pt]
    &= \left[2p^2 - 2(1-p^2) + (1-p^2)\right]\Delta\chi_1^2 \\[10pt]
    &= \left[2p^2 - (1-p^2)\right]\Delta\chi_1^2 \\[10pt]
    &= \left[2p^2 - 1 + p^2\right]\Delta\chi_1^2 \\[10pt]
    &= \left(3p^2 - 1\right)\Delta\chi_1^2 \; = \; 2\,\mathcal{L}_2(p)\,\Delta\chi_1^2 \; ,
    \quad \quad (6.2)
\end{align*}
```

``\mathcal{L}_2`` being the second Legendre polynomial. So ``F`` is **not** zero: it is
``\mathcal{O}(\Delta\chi_1^2)``. It multiplies a **finite** combination, since
``I_0^0 \rightarrow \sigma_0`` while ``I_2^0, I_4^0 \rightarrow 0``:

```math
\begin{align*}
    F \left(\frac{I_0^0}{60} + \frac{I_2^0}{42} + \frac{I_4^0}{140}\right)
    &\underset{\Delta\chi_1\rightarrow 0^{+}}{\sim}
        \left(3p^2-1\right)\Delta\chi_1^2 \left(
            \frac{\sigma_0}{60} + \frac{1}{42}\frac{\sigma_{-2}}{15}\Delta\chi_1^2
            + \frac{1}{140}\frac{\sigma_{-4}}{945}\Delta\chi_1^4 \right) \\[10pt]
    &\underset{\Delta\chi_1\rightarrow 0^{+}}{\sim}
        \underbrace{\left(3p^2-1\right)\Delta\chi_1^2}_{\mathcal{O}(\Delta\chi_1^2)}
        \cdot \underbrace{\frac{\sigma_0}{60}}_{\mathcal{O}(1)} \\[10pt]
    &= \frac{\left(3p^2-1\right)\sigma_0}{60}\,\Delta\chi_1^2 \\[10pt]
    &= \mathcal{O}(\Delta\chi_1^2)
    \; \xrightarrow[\Delta\chi_1\rightarrow 0^{+}]{} \; 0 \; .
    \quad \quad (6.3)
\end{align*}
```

### Term 2: ``J_{20} I_0^2``

With ``J_{20}`` of Eq.(D.15):

```math
\begin{align*}
    J_{20}^{\kappa \phi} I_0^2 &= \frac{y\,\Delta\chi_1^2}{2} \cdot I_0^2(\Delta\chi_1) \\[10pt]
    (2.1\mathrm{a}) \; \rightarrow \quad
        &\underset{\Delta\chi_1\rightarrow 0^{+}}{\sim}
            \frac{\overbrace{y}^{\rightarrow 1}\cancel{\Delta\chi_1^2}}{2} \,
            \frac{\sigma_2}{\cancel{\Delta\chi_1^2}} \\[10pt]
        &= \frac{\sigma_2}{2} \; .
    \quad \quad (6.4)
\end{align*}
```

### The sum

```math
\begin{align*}
(6.3) : \quad J_{\beta}^{\kappa \phi}\left(\tfrac{I_0^0}{60} + \tfrac{I_2^0}{42} + \tfrac{I_4^0}{140}\right)
    &\underset{\Delta\chi_1\rightarrow 0^{+}}{\sim} \mathcal{O}(\Delta\chi_1^2)
    &&\rightarrow 0 \\[10pt]
(6.4) : \quad J_{20}^{\kappa \phi} I_0^2
    &\underset{\Delta\chi_1\rightarrow 0^{+}}{\sim} \mathcal{O}(1)
    &&\rightarrow \frac{\sigma_2}{2} \\[10pt]
\end{align*}
```

```math
    \boxed{\;
    \lim_{\Delta\chi_1 \rightarrow 0} \left[
        J_{\beta}^{\kappa \phi} \left(\frac{I_0^0}{60} + \frac{I_2^0}{42} + \frac{I_4^0}{140}\right)
        + J_{20}^{\kappa \phi} I_0^2
    \right] = \frac{\sigma_2}{2} \; . }
    \quad \quad (6.5)
```

The ``p`` of Eq.(6.2) disappears because it multiplies a vanishing ``\Delta\chi_1^2``, so
again nothing is left to cancel.
