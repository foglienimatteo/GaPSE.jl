# Family 5: Newtonian ``\times`` Integrated GP

## Recap: everything this page needs

This page is self-contained. All the results quoted here are derived in
[The ``\Delta\chi \rightarrow 0`` limits](DeltaChiLimits.md), and the ``I_\ell^n``
asymptotics in [The ``I_\ell^n`` integrals](IlnIntegrals.md).

**The separation and its singular point.** For this family the two competing distances are
``s_1`` and ``\chi_2``, so the relevant separation is

```math
    \Delta\chi_2 := \sqrt{s_1^2 + \chi_2^2 - 2 \, s_1 \, \chi_2 \, y} \quad \quad (D.3)
```

```math
    \Delta\chi_2 = 0 \quad \iff \quad y = 1 \; \land \; s_1 = \chi_2 \quad \quad (1.7c)
```

i.e. the limit is a **joint** one: ``y \rightarrow 1`` and ``\chi_2 \rightarrow s_1`` together.

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
    \chi_2 := s_1 + p \, \Delta\chi_2  \quad \quad (D.4) \\[10pt]
    |s_1 - \chi_2| \leq \Delta\chi_2 \quad \Rightarrow \quad |p| \leq 1
```

**and its two consequences**, obtained by inverting (D.3):

```math
\Rightarrow \quad 
    y-1  \underset{\Delta\chi_2\rightarrow 0^{+}}{\sim}  - \frac{(1-p^2)}{2s_1^2}\Delta\chi_2^2
    \quad \quad (2.3\mathrm{a})
    \quad , \quad \quad \quad 
    y^2-1 \underset{\Delta\chi_2 \rightarrow 0^{+}}{\sim}- \frac{(1-p^2)}{s_1^2}\Delta\chi_2^2
    \quad \quad (2.3\mathrm{b}) \\[10pt]
```

**The trap to avoid.** The two vanishing quantities have **different orders**,

```math
    \chi_2 - s_1 = p \, \Delta\chi_2 = \mathcal{O}(\Delta\chi_2)
    \qquad \qquad
    y - 1 = \mathcal{O}(\Delta\chi_2^2)
    \quad \quad (2.4)
```

so a term **quadratic** in ``(s_1 - \chi_2)`` contributes at the **same order** as a term
**linear** in ``(y-1)``. One must therefore never set ``s_1 = \chi_2`` first and expand in ``y``
afterwards: the safe recipe is to rewrite every vanishing bracket **exactly** as a
combination of ``(s_1-\chi_2)``, ``(y-1)`` and ``(y^2-1)``, with coefficients that are regular
at the singular point, and only then substitute (D.4), (2.3a), (2.3b).

## The integrand

Concerned functions: `integrand_ξ_GNC_Newtonian_IntegratedGP`,
`integrand_ξ_GNCxLD_Newtonian_IntegratedGP`.
The GNC Newtonian-IntegratedGP function is:

```math
\begin{split}
    \xi^{\delta \int\!\phi} ( s_1 , s_2, y ) &=
    D_1 \int_0^{s_2}\mathrm{d} \chi_2 \;
    J^{\delta \int\!\phi}_{\alpha}\left[
        J^{\delta \int\!\phi}_{20} I_0^2 ( \Delta \chi_2 ) +
        \right.\nonumber \\
        &\left.\qquad
        J^{\delta \int\!\phi}_{\beta}
        \left(
            \frac{1}{15} I_0^0 ( \Delta \chi_2 ) +
            \frac{2}{21} I_2^0 ( \Delta \chi_2 ) +
            \frac{1}{35} I_4^0 ( \Delta \chi_2 )
        \right)
    \right] \, ,
\end{split}
```

with (``f_1 = f(s_1)`` the growth rate and ``b_1 = b(s_1)`` the bias)

```math
\begin{align*}
    J^{\delta \int\!\phi}_{\alpha} &=
    \frac{\mathcal{H}_0^2 \Omega_{\mathrm{M}0} D(\chi_2)}{3 a(\chi_2) s_2} 
    \left[ 
        s_2 \mathcal{R}_2 \mathcal{H}(\chi_2) ( f(\chi_2) - 1) - 5 s_{\mathrm{b}, 2} + 2
    \right] 
    \, , \\[6pt]
    J^{\delta \int\!\phi}_{\beta} &=
    \underbrace{f_1 \left[ 
        (3 y^2 - 1) \chi_2^2 - 4 y s_1 \chi_2 + 2 s_1^2
    \right]}_{=: \; F \; , \quad (\mathrm{D}.16)}
    \, , \\[6pt]
    J^{\delta \int\!\phi}_{20} &=
    \underbrace{- \Delta\chi_2^2 ( 3 b_1 + f_1)}_{(\mathrm{D}.17)}
    \, .
\end{align*}
```

The limit to be taken is:

```math
    \lim_{\Delta\chi_2\rightarrow 0^{+}}\left[
        F \left(\frac{I_0^0}{15} + \frac{2\,I_2^0}{21} + \frac{I_4^0}{35}\right)
        + J_{20} \, I_0^2 \right] \; .
```

This analysis is valid for Newtonian-IntegratedGP in both the combinations (GNC and GNCxLD).

### Term 1: the ``F`` group

The structure is identical to Family 4, only with a different pair of distances and
different numerical coefficients. ``F``, Eq.(D.16), vanishes at the singular point,

```math
    \frac{F}{f_1}\bigg|_{y=1,\,\chi_2=s_1} = \left((3-1) - 4 + 2\right)s_1^2 = (2 - 4 + 2)\,s_1^2 = 0 \; ,
```

so we decompose it exactly. For reference, the two decompositions already met are

```math
    B_{00} = 8(\chi_1-\chi_2)^2 + 8(y-1)(\chi_1^2+\chi_2^2) - 9\chi_1\chi_2(y^2-1) \; ,
    \quad \quad (3.2)
```
```math
    F^{\,(\mathrm{Family\;4})} = 2(\chi_1 - s_2)^2 + 2(y-1)\left(\chi_1^2 + s_2^2\right)
        - (y^2-1)\chi_1 s_2 \; ,
    \quad \quad (6.1)
```

and this one comes out the same way:

```math
\begin{align*}
    \frac{F}{f_1} &= (3y^2-1)\chi_2^2 - 4 y s_1 \chi_2 + 2 s_1^2 \\[10pt]
    &\quad\quad 3y^2-1 = 2 + 3(y^2-1) \\[10pt]
        &= \left[2 + 3(y^2-1)\right]\chi_2^2 - 4 y s_1 \chi_2 + 2 s_1^2 \\[10pt]
        &= 2\chi_2^2 + 2 s_1^2 + 3(y^2-1)\chi_2^2 - 4 y s_1 \chi_2 \\[10pt]
    &\quad\quad 4y = 4 + 4(y-1) \\[10pt]
        &= 2\chi_2^2 + 2 s_1^2 + 3(y^2-1)\chi_2^2
            - \left[4 + 4(y-1)\right] s_1 \chi_2 \\[10pt]
        &= 2\chi_2^2 - 4 s_1\chi_2 + 2 s_1^2
            + 3(y^2-1)\chi_2^2 - 4(y-1) s_1 \chi_2 \\[10pt]
    &\quad\quad 2\chi_2^2 - 4 s_1\chi_2 + 2 s_1^2 = 2(\chi_2 - s_1)^2 \\[10pt]
        &= 2(\chi_2 - s_1)^2 + 3(y^2-1)\chi_2^2 - 4(y-1) s_1 \chi_2 \; .
    \quad \quad (7.1)
\end{align*}
```

Again the ``(\chi_2 - s_1)`` factor is **quadratic**, hence of the same order as the
``(y-1)`` terms. Substituting (D.4), (2.3a), (2.3b) with ``\chi_2 \rightarrow s_1``:

```math
\begin{align*}
    \frac{F}{f_1} &\underset{\Delta\chi_2\rightarrow 0^{+}}{\sim}
        \underbrace{2 (-p\Delta\chi_2)^2}_{(\chi_2-s_1)^2}
        + \underbrace{3\left[- \frac{(1-p^2)}{s_1^2}\Delta\chi_2^2\right]s_1^2}_{(y^2-1)}
        - \underbrace{4\left[- \frac{(1-p^2)}{2 s_1^2}\Delta\chi_2^2\right]s_1^2}_{(y-1)} \\[10pt]
    &= 2p^2\Delta\chi_2^2 - 3(1-p^2)\Delta\chi_2^2 + 2(1-p^2)\Delta\chi_2^2 \\[10pt]
    &= \left[2p^2 - 3(1-p^2) + 2(1-p^2)\right]\Delta\chi_2^2 \\[10pt]
    &= \left[2p^2 - (1-p^2)\right]\Delta\chi_2^2 \\[10pt]
    &= \left[2p^2 - 1 + p^2\right]\Delta\chi_2^2 \\[10pt]
    &= \left(3p^2 - 1\right)\Delta\chi_2^2 \; = \; 2\,\mathcal{L}_2(p)\,\Delta\chi_2^2 \; ,
    \quad \quad (7.2)
\end{align*}
```

the very same ``2\mathcal{L}_2(p)`` of Eq.(6.2), reached through different coefficients.
Hence

```math
\begin{align*}
    F \left(\frac{I_0^0}{15} + \frac{2 I_2^0}{21} + \frac{I_4^0}{35}\right)
    &\underset{\Delta\chi_2\rightarrow 0^{+}}{\sim}
        \underbrace{f_1\left(3p^2-1\right)\Delta\chi_2^2}_{\mathcal{O}(\Delta\chi_2^2)}
        \cdot \underbrace{\frac{\sigma_0}{15}}_{\mathcal{O}(1)} \\[10pt]
    &= \frac{f_1\left(3p^2-1\right)\sigma_0}{15}\,\Delta\chi_2^2 \\[10pt]
    &= \mathcal{O}(\Delta\chi_2^2)
    \; \xrightarrow[\Delta\chi_2\rightarrow 0^{+}]{} \; 0 \; .
    \quad \quad (7.3)
\end{align*}
```

### Term 2: ``J_{20} I_0^2``

With ``J_{20}`` of Eq.(D.17):

```math
\begin{align*}
    J_{20}^{\delta \int\!\phi} I_0^2 &= -\Delta\chi_2^2 (3 b_1+f_1) \cdot I_0^2(\Delta\chi_2) \\[10pt]
    (2.1\mathrm{a}) \; \rightarrow \quad
        &\underset{\Delta\chi_2\rightarrow 0^{+}}{\sim}
            -\cancel{\Delta\chi_2^2}(3 b_1+f_1)\,\frac{\sigma_2}{\cancel{\Delta\chi_2^2}} \\[10pt]
        &= -(3 b_1+f_1)\,\sigma_2 \; .
    \quad \quad (7.4)
\end{align*}
```

### The sum

```math
\begin{align*}
(7.3) : \quad J_{\beta}^{\delta \int\!\phi}\left(\tfrac{I_0^0}{15} + \tfrac{2 I_2^0}{21} + \tfrac{I_4^0}{35}\right)
    &\underset{\Delta\chi_2\rightarrow 0^{+}}{\sim} \mathcal{O}(\Delta\chi_2^2)
    &&\rightarrow 0 \\[10pt]
(7.4) : \quad J_{20}^{\delta \int\!\phi} I_0^2
    &\underset{\Delta\chi_2\rightarrow 0^{+}}{\sim} \mathcal{O}(1)
    &&\rightarrow -(3 b_1 + f_1)\,\sigma_2 \\[10pt]
\end{align*}
```

```math
    \boxed{\;
    \lim_{\Delta\chi_2 \rightarrow 0} \left[
        J_{\beta}^{\delta \int\!\phi} \left(\frac{I_0^0}{15} + \frac{2 I_2^0}{21} + \frac{I_4^0}{35}\right)
        + J_{20}^{\delta \int\!\phi} I_0^2
    \right] = -(3 b_1 + f_1) \, \sigma_2 \; . }
    \quad \quad (7.5)
```
