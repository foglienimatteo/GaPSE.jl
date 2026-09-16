# Family 7: the ``\Delta\chi^2 \times`` (vanishing factor) terms

## Recap: everything this page needs

This page is self-contained. All the results quoted here are derived in
[The ``\Delta\chi \rightarrow 0`` limits](DeltaChiLimits.md), and the ``I_\ell^n``
asymptotics in [The ``I_\ell^n`` integrals](IlnIntegrals.md).

**The separation and its singular point.** For this family the two competing distances are
``s_1`` and ``\chi_2``, so the relevant separation is

```math
    \Delta\chi_2 := \sqrt{s_1^2 + \chi_2^2 - 2 \, s_1 \, \chi_2 \, y} \quad \quad (1.3)
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
    \chi_2 := s_1 + p \, \Delta\chi_2  \quad \quad (2.2) \\[10pt]
    |s_1 - \chi_2| \leq \Delta\chi_2 \quad \Rightarrow \quad |p| \leq 1
```

**and its two consequences**, obtained by inverting (1.3):

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
at the singular point, and only then substitute (2.2), (2.3a), (2.3b).

## The integrand

Concerned functions: `integrand_ξ_GNC_Doppler_IntegratedGP`,
`integrand_ξ_GNCxLD_IntegratedGP_Doppler`, `integrand_ξ_GNCxLD_Doppler_IntegratedGP`,
`integrand_ξ_LD_Doppler_IntegratedGP`.

These carry an overall ``\Delta\chi_2^2`` together with a geometric factor
``(\chi_2 y - s_1)`` (or ``(s_1 - \chi_2 y)``, with the opposite sign, which changes nothing
since the limit is zero), multiplying the whole parenthesis:

```math
    \Delta\chi_2^2 \, \underbrace{(\chi_2 y - s_1)}_{=: \; G} \left(
        \frac{I_0^0}{15} + \frac{2 I_2^0}{21} + \frac{I_4^0}{35} + I_0^2 \right) \; .
```

### Term 1: the geometric factor ``G``

``G`` vanishes at the singular point, ``y \rightarrow 1`` and ``\chi_2 \rightarrow s_1``
giving ``s_1 - s_1 = 0``, so it has to be decomposed — but this time the decomposition is a
single line, and the ``(\chi_2 - s_1)`` factor appears **linearly**:

```math
\begin{align*}
    G &= \chi_2 y - s_1 \\[10pt]
    y = 1 + (y-1) \; \rightarrow \quad
        &= \chi_2 + \chi_2 (y-1) - s_1 \\[10pt]
        &= \underbrace{(\chi_2 - s_1)}_{\mathcal{O}(\Delta\chi_2)}
            + \underbrace{\chi_2 (y-1)}_{\mathcal{O}(\Delta\chi_2^2)} \; .
    \quad \quad (9.1)
\end{align*}
```

By (2.4) the first term dominates — compare with the ``B_{00}`` of Family 1,

```math
    B_{00} = \underbrace{8(\chi_1-\chi_2)^2}_{\mathcal{O}(\Delta\chi^2)}
        + \underbrace{8(y-1)(\chi_1^2+\chi_2^2) - 9\chi_1\chi_2(y^2-1)}_{\mathcal{O}(\Delta\chi^2)} \; ,
    \quad \quad (3.2)
```

where the ``(\chi_1-\chi_2)`` factor was *quadratic* and therefore of the **same** order as
the ``(y-1)`` terms. Here it is linear, so it wins:

```math
\begin{align*}
    G &\underset{\Delta\chi_2\rightarrow 0^{+}}{\sim}
        p \, \Delta\chi_2 + \cancel{s_1}\left[- \frac{(1-p^2)}{2 s_1^{\cancel{2}}}
            \Delta\chi_2^2\right] \\[10pt]
    &= p \, \Delta\chi_2 - \frac{(1-p^2)}{2 s_1}\Delta\chi_2^2 \\[10pt]
    &\underset{\Delta\chi_2\rightarrow 0^{+}}{\sim} \; p \, \Delta\chi_2 \; .
    \quad \quad (9.2)
\end{align*}
```

### Term 2: against the ``I_\ell^n``

Of the four ``I_\ell^n`` in the parenthesis, ``I_0^0 \rightarrow \sigma_0`` is finite and
``I_2^0, I_4^0 \rightarrow 0``, so the only one that matters is the divergent
``I_0^2 \sim \sigma_2 \Delta\chi_2^{-2}``:

```math
\begin{align*}
    \Delta\chi_2^2 \, G \left(
        \frac{I_0^0}{15} + \frac{2 I_2^0}{21} + \frac{I_4^0}{35} + I_0^2 \right)
    &\underset{\Delta\chi_2\rightarrow 0^{+}}{\sim}
        \cancel{\Delta\chi_2^2} \cdot p \, \Delta\chi_2 \cdot
        \frac{\sigma_2}{\cancel{\Delta\chi_2^2}} \\[10pt]
    &= p \, \sigma_2 \, \Delta\chi_2
    \; \xrightarrow[\Delta\chi_2\rightarrow 0^{+}]{} \; 0 \; .
\end{align*}
```

### The result

```math
    \boxed{\; \lim_{\Delta\chi_2 \rightarrow 0}
    \left[\Delta\chi_2^2 (\chi_2 y - s_1) \left(
        \frac{I_0^0}{15} + \frac{2 I_2^0}{21} + \frac{I_4^0}{35} + I_0^2 \right)\right]
    = 0 \; . }
    \quad \quad (9.3)
```

The rate of approach is ``\mathcal{O}(\Delta\chi_2)`` and carries a ``p``, but since the limit
itself is ``0`` the direction-independence is automatic.
