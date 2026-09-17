# Family 6: the ``\Delta\chi^4 \, \tilde{I}_0^4`` terms

## Recap: everything this page needs

This page is self-contained. All the results quoted here are derived in
[The ``\Delta\chi \rightarrow 0`` limits](DeltaChiLimits.md), and the ``I_\ell^n``
asymptotics in [The ``I_\ell^n`` integrals](IlnIntegrals.md).

**The separation and its singular point.** For this family the two competing distances are
``\chi_1`` and ``\chi_2``, so the relevant separation is

```math
    \Delta\chi := \sqrt{\chi_1^2 + \chi_2^2 - 2 \, \chi_1 \, \chi_2 \, y} \quad \quad (D.1)
```

```math
    \Delta\chi = 0 \quad \iff \quad y = 1 \; \land \; \chi_1 = \chi_2 \quad \quad (1.7a)
```

i.e. the limit is a **joint** one: ``y \rightarrow 1`` and ``\chi_2 \rightarrow \chi_1`` together.

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
    \chi_2 := \chi_1 + p \, \Delta\chi  \quad \quad (D.4) \\[10pt]
    |\chi_1 - \chi_2| \leq \Delta\chi \quad \Rightarrow \quad |p| \leq 1
```

**and its two consequences**, obtained by inverting (D.1):

```math
\Rightarrow \quad 
    y-1  \underset{\Delta\chi\rightarrow 0^{+}}{\sim}  - \frac{(1-p^2)}{2\chi_1^2}\Delta\chi^2
    \quad \quad (2.3\mathrm{a})
    \quad , \quad \quad \quad 
    y^2-1 \underset{\Delta\chi \rightarrow 0^{+}}{\sim}- \frac{(1-p^2)}{\chi_1^2}\Delta\chi^2
    \quad \quad (2.3\mathrm{b}) \\[10pt]
```

**The trap to avoid.** The two vanishing quantities have **different orders**,

```math
    \chi_2 - \chi_1 = p \, \Delta\chi = \mathcal{O}(\Delta\chi)
    \qquad \qquad
    y - 1 = \mathcal{O}(\Delta\chi^2)
    \quad \quad (2.4)
```

so a term **quadratic** in ``(\chi_1 - \chi_2)`` contributes at the **same order** as a term
**linear** in ``(y-1)``. One must therefore never set ``\chi_1 = \chi_2`` first and expand in ``y``
afterwards: the safe recipe is to rewrite every vanishing bracket **exactly** as a
combination of ``(\chi_1-\chi_2)``, ``(y-1)`` and ``(y^2-1)``, with coefficients that are regular
at the singular point, and only then substitute (D.4), (2.3a), (2.3b).

## The integrand

Concerned functions: `integrand_ξ_GNC_IntegratedGP`, `integrand_ξ_LD_IntegratedGP`,
`integrand_ξ_GNCxLD_IntegratedGP_IntegratedGP`, `integrand_ξ_GNC_LocalGP_IntegratedGP`,
`integrand_ξ_GNCxLD_IntegratedGP_LocalGP`, `integrand_ξ_GNCxLD_LocalGP_IntegratedGP`,
`integrand_ξ_LD_LocalGP_IntegratedGP`.

All of them contain a single term,

```math
    \Delta\chi^4 \, \tilde{I}_0^4(\Delta\chi) \; ,
```

everything else being regular at the singular point. This is the simplest family of all: no
bracket vanishes, no cancellation is involved, and the direction parameter ``p`` never
enters. Only the explicit ``\Delta\chi^4`` against the ``\Delta\chi^{-2}`` of (2.1b) matters.

### Term 1: ``\Delta\chi^4 \, \tilde{I}_0^4``

The full series of ``\tilde{I}_0^4``, derived in
[The ``I_\ell^n`` integrals](IlnIntegrals.md), is

```math
    \tilde{I}_0^4(s) = \sum_{k=1}^{+\infty}
        \frac{(-1)^k \, \sigma_{4-2k}}{(2k+1)!} \, s^{\,2k-4}
    = -\frac{\sigma_2}{6\,s^2} + \frac{\sigma_0}{120} - \frac{\sigma_{-2}}{5040}\,s^2 + \dots
```

whose leading term is exactly Eq.(2.1b). Multiplying by ``\Delta\chi^4``:

```math
\begin{align*}
    \Delta\chi^4 \, \tilde{I}_0^4(\Delta\chi)
    &= \Delta\chi^4 \left(-\frac{\sigma_2}{6\,\Delta\chi^2} + \frac{\sigma_0}{120}
        - \frac{\sigma_{-2}}{5040}\Delta\chi^2 + \dots \right) \\[10pt]
    &= -\frac{\sigma_2}{6}\Delta\chi^2 + \frac{\sigma_0}{120}\Delta\chi^4
        - \frac{\sigma_{-2}}{5040}\Delta\chi^6 + \dots \\[10pt]
    &\underset{\Delta\chi\rightarrow 0^{+}}{\sim} -\frac{\sigma_2}{6}\Delta\chi^2
    \; \xrightarrow[\Delta\chi \rightarrow 0^{+}]{} \; 0 \; .
    \quad \quad (8.1)
\end{align*}
```

Note that the ``\Delta\chi^{4}`` overwhelms the ``\Delta\chi^{-2}`` divergence by two orders,
so the limit is not merely finite but **zero**, and it is approached as ``\Delta\chi^2``.

### The result

```math
    \boxed{\; \lim_{\Delta\chi \rightarrow 0}
        \left(\Delta\chi^4 \, \tilde{I}_0^4(\Delta\chi)\right) = 0 \; . }
    \quad \quad (8.2)
```
