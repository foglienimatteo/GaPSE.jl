# Family 8: the ``J_{22} I_2^2 + J_{31} I_1^3`` terms

## Recap: everything this page needs

This page is self-contained. All the results quoted here are derived in
[The ``\Delta\chi \rightarrow 0`` limits](DeltaChiLimits.md), and the ``I_\ell^n``
asymptotics in [The ``I_\ell^n`` integrals](IlnIntegrals.md).

**The separation and its singular point.** For this family the two competing distances are
``\chi_a`` and ``\chi_b``, so the relevant separation is

```math
    \Delta\chi := \sqrt{\chi_a^2 + \chi_b^2 - 2 \, \chi_a \, \chi_b \, y} \quad \quad (D.1)
```

```math
    \Delta\chi = 0 \quad \iff \quad y = 1 \; \land \; \chi_a = \chi_b \quad \quad (1.7a)
```

i.e. the limit is a **joint** one: ``y \rightarrow 1`` and ``\chi_b \rightarrow \chi_a`` together.

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
    \chi_b := \chi_a + p \, \Delta\chi  \quad \quad (D.4) \\[10pt]
    |\chi_a - \chi_b| \leq \Delta\chi \quad \Rightarrow \quad |p| \leq 1
```

**and its two consequences**, obtained by inverting (D.1):

```math
\Rightarrow \quad 
    y-1  \underset{\Delta\chi\rightarrow 0^{+}}{\sim}  - \frac{(1-p^2)}{2\chi_a^2}\Delta\chi^2
    \quad \quad (2.3\mathrm{a})
    \quad , \quad \quad \quad 
    y^2-1 \underset{\Delta\chi \rightarrow 0^{+}}{\sim}- \frac{(1-p^2)}{\chi_a^2}\Delta\chi^2
    \quad \quad (2.3\mathrm{b}) \\[10pt]
```

**The trap to avoid.** The two vanishing quantities have **different orders**,

```math
    \chi_b - \chi_a = p \, \Delta\chi = \mathcal{O}(\Delta\chi)
    \qquad \qquad
    y - 1 = \mathcal{O}(\Delta\chi^2)
    \quad \quad (2.4)
```

so a term **quadratic** in ``(\chi_a - \chi_b)`` contributes at the **same order** as a term
**linear** in ``(y-1)``. One must therefore never set ``\chi_a = \chi_b`` first and expand in ``y``
afterwards: the safe recipe is to rewrite every vanishing bracket **exactly** as a
combination of ``(\chi_a-\chi_b)``, ``(y-1)`` and ``(y^2-1)``, with coefficients that are regular
at the singular point, and only then substitute (D.4), (2.3a), (2.3b).

## The integrand

Concerned functions: `integrand_ξ_GNC_Lensing_IntegratedGP`,
`integrand_ξ_GNCxLD_IntegratedGP_Lensing`, `integrand_ξ_GNCxLD_Lensing_IntegratedGP`,
`integrand_ξ_GNCxLD_Lensing_LocalGP`, `integrand_ξ_LD_Lensing_IntegratedGP`,
`integrand_ξ_LD_Lensing_LocalGP`.
The GNC Lensing-IntegratedGP function is:

```math
\begin{split}
    \xi^{\kappa \int\!\phi} ( s_1 , s_2, y ) = 
    \int_0^{s_1}\mathrm{d} \chi_1 \int_0^{s_2}\mathrm{d} \chi_2 \;
    J_{\alpha}^{\kappa \int\!\phi} 
    \left[ 
        J_{31}^{\kappa \int\!\phi} I_1^3 ( \Delta \chi ) +
        J_{22}^{\kappa \int\!\phi} I_2^2 ( \Delta \chi ) 
     \right] \, ,
\end{split}
```

with

```math
\begin{align*}
     J_{\alpha}^{\kappa \int\!\phi} &=
    \frac{
        9 \chi_2 \ \mathcal{H}_0^4  \Omega_{\mathrm{M}0}^2 D(\chi_1) D(\chi_2) 
    }{
        a(\chi_1)  a(\chi_2) s_1  s_2
    }
    (\chi_1 - s_1) (5 s_{\mathrm{b}, 1} - 2) \times \\
    &\qquad\qquad\qquad\qquad
    \left[
        \ \mathcal{H}(\chi_2)  \mathcal{R}_2 s_1 (f(\chi_2) - 1) - 5 s_{\mathrm{b}, 1} + 2
    \right]
    \, , \\[6pt]
    J_{31}^{\kappa \int\!\phi} &=  y \, \Delta\chi^2
    \, , \\[6pt]
    J_{22}^{\kappa \int\!\phi} &= 
    \frac{1}{2} (y^2 - 1) \chi_1 \chi_2 
    \, .
\end{align*}
```

The six functions share this structure up to an overall constant ``A``, which we factor out
by writing (Eq.(D.19))

```math
    J_{22} := \frac{A}{2}\,\chi_a \chi_b \, (y^2-1) \; , \qquad
    J_{31} := A \, y \, \Delta\chi^2 \; ,
    \quad \quad (\mathrm{D}.19)
```

with

| ``A`` | functions |
|:-:|:--|
| ``1`` | `integrand_ξ_GNC_Lensing_IntegratedGP` |
| ``2`` | `integrand_ξ_GNCxLD_IntegratedGP_Lensing` , `integrand_ξ_GNCxLD_Lensing_IntegratedGP` |
| ``-2`` | `integrand_ξ_GNCxLD_Lensing_LocalGP` , `integrand_ξ_LD_Lensing_IntegratedGP` , `integrand_ξ_LD_Lensing_LocalGP` |

(in the last three the coefficients are written in the source as ``-2y\Delta\chi^2`` and
``\chi_a\chi_b(1-y^2)``, which is the same thing with ``A = -2``). The limit to be taken is:

```math
    \lim_{\Delta\chi\rightarrow 0^{+}}\left(J_{22} \, I_2^2 + J_{31} \, I_1^3\right) \; .
```

This analysis is valid for all the six functions listed above.

### Term 1: ``J_{22} I_2^2``

``J_{22}`` of Eq.(D.19) carries **no** negative power of ``\Delta\chi``, and its only
vanishing factor is ``(y^2-1)`` on its own — there is no bracket to decompose, because
nothing else cancels against it. Compare with the ``B_{00}`` of Family 1, Eq.(D.5) there,
where

```math
    B_{00} = 8(\chi_1-\chi_2)^2 + 8(y-1)(\chi_1^2+\chi_2^2) - 9\chi_1\chi_2(y^2-1) \; :
    \quad \quad (3.2)
```

there the ``(y^2-1)`` was one of three competing contributions, here it stands alone, so
(2.3b) can be applied directly:

```math
\begin{align*}
    J_{22}^{\kappa \int\!\phi} &= \frac{A}{2}\,\chi_a \chi_b \, (y^2-1) \\[10pt]
    (2.3\mathrm{b}) \; \rightarrow \quad
        &\underset{\Delta\chi\rightarrow 0^{+}}{\sim}
            \frac{A}{2}\,\chi_a \chi_b
            \left[- \frac{(1-p^2)}{\chi_a^2}\Delta\chi^2\right] \\[10pt]
    \chi_b \rightarrow \chi_a \; \rightarrow \quad
        &= \frac{A}{2}\cancel{\chi_a^2}
            \left[- \frac{(1-p^2)}{\cancel{\chi_a^2}}\Delta\chi^2\right] \\[10pt]
        &= -\frac{A}{2}(1-p^2)\,\Delta\chi^2 \; = \; \mathcal{O}(\Delta\chi^2) \; .
\end{align*}
```

Since ``I_2^2 \rightarrow \sigma_0/15`` is **finite**, the product vanishes:

```math
\begin{align*}
    J_{22}^{\kappa \int\!\phi} I_2^2
    &\underset{\Delta\chi\rightarrow 0^{+}}{\sim}
        \underbrace{\left[-\frac{A}{2}(1-p^2)\,\Delta\chi^2\right]}_{\mathcal{O}(\Delta\chi^2)}
        \cdot \underbrace{\frac{\sigma_0}{15}}_{\mathcal{O}(1)} \\[10pt]
    &= -\frac{A\,(1-p^2)\,\sigma_0}{30}\,\Delta\chi^2 \\[10pt]
    &= \mathcal{O}(\Delta\chi^2)
    \; \xrightarrow[\Delta\chi\rightarrow 0^{+}]{} \; 0 \; .
    \quad \quad (10.1)
\end{align*}
```

### Term 2: ``J_{31} I_1^3``

``J_{31}``, also defined in Eq.(D.19), is regular and its ``\Delta\chi^2`` cancels the
``\Delta\chi^{-2}`` of ``I_1^3``:

```math
\begin{align*}
    J_{31}^{\kappa \int\!\phi} I_1^3 &= A \, y \, \Delta\chi^2 \cdot I_1^3(\Delta\chi) \\[10pt]
    (2.1\mathrm{a}) \; \rightarrow \quad
        &\underset{\Delta\chi\rightarrow 0^{+}}{\sim}
            A \, y \, \Delta\chi^2 \cdot \frac{\sigma_2}{3\,\Delta\chi^2} \\[10pt]
    y \rightarrow 1 \; \rightarrow \quad
        &= A \underbrace{y}_{\rightarrow 1}\frac{\cancel{\Delta\chi^2}}{3\cancel{\Delta\chi^2}}\,\sigma_2 \\[10pt]
        &= \frac{A}{3}\,\sigma_2 \; .
    \quad \quad (10.2)
\end{align*}
```

### The sum

```math
\begin{align*}
(10.1) : \quad J_{22}^{\kappa \int\!\phi} I_2^2
    &\underset{\Delta\chi\rightarrow 0^{+}}{\sim} \mathcal{O}(\Delta\chi^2)
    &&\rightarrow 0 \\[10pt]
(10.2) : \quad J_{31}^{\kappa \int\!\phi} I_1^3
    &\underset{\Delta\chi\rightarrow 0^{+}}{\sim} \mathcal{O}(1)
    &&\rightarrow \frac{A}{3}\,\sigma_2 \\[10pt]
\end{align*}
```

```math
    \boxed{\; \lim_{\Delta\chi \rightarrow 0}
    \left(J_{22}^{\kappa \int\!\phi} I_2^2 + J_{31}^{\kappa \int\!\phi} I_1^3\right)
    = \frac{A}{3}\,\sigma_2 \; . }
    \quad \quad (10.3)
```

The only ``p`` appears in the term whose limit is ``0``, so nothing has to cancel.
