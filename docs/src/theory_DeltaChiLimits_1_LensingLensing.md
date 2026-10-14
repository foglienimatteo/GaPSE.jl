# Family 1: Lensing ``\times`` Lensing

## Recap: everything this page needs

This page is self-contained. All the results quoted here are derived in
[The ``\Delta\chi \rightarrow 0`` limits](theory_DeltaChiLimits.md), and the ``I_\ell^n``
asymptotics in [The ``I_\ell^n`` integrals](theory_IlnIntegrals.md).

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

Concerned functions: `integrand_ξ_GNC_Lensing`, `integrand_ξ_LD_Lensing`,
`integrand_ξ_GNCxLD_Lensing_Lensing`.

The analytical expression of the function `integrand_ξ_GNC_Lensing` is the following:

```math
f^{\kappa\kappa} (\chi_1, \chi_2, s_1, s_2, y) = 
J^{\kappa\kappa}_{\alpha}
\left[
    J^{\kappa\kappa}_{00} I_0^0(\Delta\chi) + 
    J^{\kappa\kappa}_{02} I_2^0(\Delta\chi) +
    J^{\kappa\kappa}_{31} I_1^3(\Delta\chi) +
    J^{\kappa\kappa}_{22} I_2^2(\Delta\chi)
\right]  \, , 
```

with

```math
\begin{align*}
    J^{\kappa\kappa}_{\alpha} & = 
    \frac{
        \mathcal{H}_0^4 \Omega_{\mathrm{M}0}^2 D(\chi_1) D(\chi_2) 
    }{
        s_1 s_2 a(\chi_1) a(\chi_2)}
    (\chi_1 - s_1)(\chi_2 - s_2)
    (5 s_{\mathrm{b}, 1} - 2)(5 s_{\mathrm{b}, 2} - 2) 
    \, , \\
    %%%%&%%%%%%%%%%%%%
    J^{\kappa\kappa}_{00} & = 
    -\frac{3}{4}\frac{\chi_1^2 \chi_2^2}{\Delta\chi^4} (y^2 - 1)
    \left[
        8 y (\chi_1^2 + \chi_2^2) - \chi_1\chi_2 (9y^2+7)
    \right] 
    \, , \\
    %%%%&%%%%%%%%%%%%%
    J^{\kappa\kappa}_{02} & = 
    -\frac{3}{2}\frac{\chi_1^2 \chi_2^2}{\Delta\chi^4}(y^2 - 1)
    \left[
        4 y (\chi_1^2 + \chi_2^2) - \chi_1 \chi_2 (3 y^2 + 5)
    \right] 
    \, , \\
    %%%%%%%%%%%%%%%%%%
    J^{\kappa\kappa}_{31} & = 9 y \Delta\chi^2 
    \, , \\
    %%%%%%%%%%%%%%%%%%
    J^{\kappa\kappa}_{22} & = 
    \frac{9}{4}\frac{\chi_1 \chi_2}{\Delta\chi^4}
    \left[
        2(\chi_1^4 + \chi_2^4)(7 y^2 - 3) - 
        16 y \chi_1 \chi_2 (\chi_1^2 + \chi_2^2)(y^2 + 1) + 
        \chi_1^2 \chi_2^2 (11y^4 + 14y^2 + 23) 
    \right] 
    \, .
\end{align*}
```



The limit to be taken is:

```math
    \lim_{\Delta\chi\rightarrow 0^{+}}
    \left(J_{00}I_0^0 + J_{02}I_2^0 + J_{31}I_1^3 + J_{22}I_2^2\right) \; .
```

This analysis is valid for Lensing-Lensing in all the 3 combinations (GNC, LD and GNCxLD).

We abbreviate the three square brackets as

```math
\begin{align*}
    B_{00} &:= 8 y (\chi_1^2 + \chi_2^2) - \chi_1\chi_2 (9y^2+7)
        &&\quad \quad (\mathrm{D}.5) \\[6pt]
    B_{02} &:= 4 y (\chi_1^2 + \chi_2^2) - \chi_1\chi_2 (3y^2+5)
        &&\quad \quad (\mathrm{D}.6) \\[6pt]
    B_{22} &:= 2(\chi_1^4 + \chi_2^4)(7 y^2 - 3) - 16 y \chi_1 \chi_2 (\chi_1^2 + \chi_2^2)(y^2 + 1)
        + \chi_1^2 \chi_2^2 (11y^4 + 14y^2 + 23)
        &&\quad \quad (\mathrm{D}.7)
\end{align*}
```

so that

```math
    J^{\kappa\kappa}_{00} = -\frac{3}{4}\frac{\chi_1^2\chi_2^2}{\Delta\chi^4}(y^2-1) B_{00}
    \; , \quad
    J^{\kappa\kappa}_{02} = -\frac{3}{2}\frac{\chi_1^2\chi_2^2}{\Delta\chi^4}(y^2-1) B_{02}
    \; , \quad
    J^{\kappa\kappa}_{22} = \frac{9}{4}\frac{\chi_1\chi_2}{\Delta\chi^4} B_{22} \; .
    \quad \quad (3.1)
```

All three brackets (D.5), (D.6), (D.7) vanish at the singular point ``y=1 \land \chi_1 = \chi_2 = \chi``:

```math
\begin{align*}
y=1 \; &\land \; \chi_1 = \chi_2 = \chi \; : \\[10pt]
    B_{00} &= 8(2\chi^2) - \chi^2 (9+7) = 16\chi^2 - 16\chi^2 = 0 \; , \\
    B_{02} &= 4(2\chi^2) - \chi^2 (3+5) = 8\chi^2 - 8\chi^2 = 0 \; , \\
    B_{22} &= 2(2\chi^4)(7-3) - 16\chi^2(2\chi^2)(1+1) + \chi^4(11+14+23)
        = 16\chi^4 - 64\chi^4 + 48\chi^4 = 0 \; ,
\end{align*}
```

so the leading order is not enough and the expansion must be pushed further.


### Term 1: ``J_{00} I_0^0``

We start by decomposing ``B_{00}``, Eq.(D.5), exactly:



```math
\begin{align*}
    B_{00} &= 8 y (\chi_1^2 + \chi_2^2) - \chi_1\chi_2 (9y^2+7) \\[10pt]
    8y = 8 + 8(y-1) \; \rightarrow \quad
        &= [8+8(y-1)](\chi_1^2+\chi_2^2) - \chi_1\chi_2 (9y^2+7) \\[10pt]
        &= 8(\chi_1^2+\chi_2^2) + 8(y-1)(\chi_1^2+\chi_2^2) - \chi_1\chi_2 (9y^2+7) \\[10pt]
    9y^2+7 = 16 + 9(y^2-1) \; \rightarrow \quad
        &= 8(\chi_1^2+\chi_2^2) + 8(y-1)(\chi_1^2+\chi_2^2) - \chi_1\chi_2[16+9(y^2-1)] \\[10pt]
        &= 8(\chi_1^2+\chi_2^2) - 16\chi_1\chi_2
            + 8(y-1)(\chi_1^2+\chi_2^2) - 9\chi_1\chi_2(y^2-1) \\[10pt]
        &= 8(\chi_1-\chi_2)^2
            + 8(y-1)(\chi_1^2+\chi_2^2) - 9\chi_1\chi_2(y^2-1) \; . \quad \quad (3.2)\\[10pt]
    (D.4),\; (2.3\mathrm{a}),\; &(2.3\mathrm{b})\; \quad\rightarrow \quad
        (\chi_2-\chi_1)=p\Delta\chi \; , \quad (\chi_1^2+\chi_2^2)=2\chi_1^2\\[10pt]
    &\underset{\Delta\chi\rightarrow 0^{+}}{\sim}
        8 \, p^2\Delta\chi^2
        + 8 \left[- \frac{(1-p^2)}{2\cancel{\chi_1^2}}\Delta\chi^2\right] 2\cancel{\chi_1^2}
        - 9 \cancel{\chi_1^2} \left[- \frac{(1-p^2)}{\cancel{\chi_1^2}}\Delta\chi^2\right] \\[10pt]
    &= \left[ 8p^2 - 8(1-p^2) + 9(1-p^2) \right] \Delta\chi^2 \\[10pt]
    &= \left[ 8p^2 + (1-p^2) \right] \Delta\chi^2 \\[10pt]
    &= \left( 1 + 7p^2 \right) \Delta\chi^2 \; . \quad \quad (3.3) \\[15pt]
\end{align*}
```


```math
\begin{align*}
    \Rightarrow \quad 
    J_{00} I_0^0 &= -\frac{3}{4}\frac{\chi_1^2\chi_2^2}{\Delta\chi^4}(y^2-1) \, B_{00} \, I_0^0 \\[10pt]
    (2.1\mathrm{a})\, ,\;(D.4)\, , \;  (2.3\mathrm{a}) \rightarrow \quad
    &\underset{\Delta\chi\rightarrow 0^{+}}{\sim}
        -\frac{3}{4}\frac{\chi_1^4}{\Delta\chi^4}
        \left[- \frac{(1-p^2)}{\chi_1^2}\Delta\chi^2\right]
        \left( 1 + 7p^2 \right) \Delta\chi^2 \; \sigma_0 \\[10pt]
    &= \frac{3}{4}\frac{\chi_1^4}{\cancel{\Delta\chi^4}}
        \frac{(1-p^2)}{\chi_1^2}\cancel{\Delta\chi^2}
        \left( 1 + 7p^2 \right) \cancel{\Delta\chi^2} \; \sigma_0 \\[10pt]
    &= \frac{3}{4}\,\chi_1^2 \, \sigma_0 \, (1-p^2)(1+7p^2) \\[10pt]
    &= \frac{3}{4}\,\chi_1^2 \, \sigma_0 \, \left(-7p^4 + 6p^2 + 1\right) \; . \quad \quad (3.4)
\end{align*}
```

### Term 2: ``J_{02} I_2^0``

We can do for ``B_{02}``, Eq.(D.6), exactly the same decomposition we did for ``B_{00}``,
Eq.(D.5), with ``8 \rightarrow 4`` and ``9 \rightarrow 3``:

```math
\begin{align*}
    B_{00} &:= 8 y (\chi_1^2 + \chi_2^2) - \chi_1\chi_2 (9y^2+7) \\[10pt]
    (3.2) \rightarrow \quad &= 8(\chi_1-\chi_2)^2
            + 8(y-1)(\chi_1^2+\chi_2^2) - 9\chi_1\chi_2(y^2-1) \\[15pt]
        &\underset{\Delta\chi\rightarrow 0^{+}}{\sim} \left[ 8p^2 - 8(1-p^2) + 9(1-p^2) \right] \Delta\chi^2 \\[10pt]
        &= \left( 1 + 7p^2 \right) \Delta\chi^2 \\[15pt]
    \Rightarrow \quad
    B_{02} &:= 4 y (\chi_1^2 + \chi_2^2) - \chi_1\chi_2 (3y^2+5) \\[10pt]
    &= 4(\chi_1-\chi_2)^2 + 4(y-1)(\chi_1^2+\chi_2^2) - 3\chi_1\chi_2(y^2-1) \\[10pt]
    
    &\underset{\Delta\chi\rightarrow 0^{+}}{\sim}
        \left[ 4p^2 - 4(1-p^2) + 3(1-p^2) \right] \Delta\chi^2 \\[10pt]
    &= \left( 5p^2 - 1 \right) \Delta\chi^2 \; . \quad \quad (3.5)
\end{align*}
```

Consequently, ``J_{02}`` is finite, ``\mathcal{O}(\Delta\chi^0)``, while ``I_2^0`` vanishes:

```math
\begin{align*}
    J_{02} I_2^0 &= -\frac{3}{2}\frac{\chi_1^2\chi_2^2}{\Delta\chi^4}(y^2-1) B_{02} I_2^0
    \; , \quad \\[10pt]
    &\underset{\Delta\chi\rightarrow 0^{+}}{\sim}
        -\frac{3}{2}\frac{\chi_1^4}{\cancel{\Delta\chi^4}} 
        \left[- \frac{(1-p^2)}{\chi_1^2}\cancel{\Delta\chi^2}\right]\left( 5p^2 - 1 \right) 
        \cancel{\Delta\chi^2} \cdot \frac{\sigma_{-2}}{15}\Delta\chi^2\\[10pt]
    &= \underbrace{
        \frac{3}{2}\,\chi_1^2 \,(1-p^2)(5p^2-1)
        }_{\mathcal{O}(\Delta\chi^0)}
        \cdot
        \underbrace{\frac{\sigma_{-2}}{15}\Delta\chi^2}_{\rightarrow \, 0} \\[25pt]
    &\xrightarrow[\Delta\chi \rightarrow 0^{+}]{} \; 0 \; . \quad \quad (3.6)
\end{align*}
```

### Term 3: ``J_{31} I_1^3``

This is the only regular ``J`` of the family, ``J_{31} = 9 y \Delta\chi^2``: nothing
vanishes, nothing has to be decomposed, and the ``\Delta\chi^2`` it carries is exactly what
cancels the ``\Delta\chi^{-2}`` of ``I_1^3``:

```math
\begin{align*}
    J_{31} I_1^3 &= 9 y \Delta\chi^2 \cdot I_1^3(\Delta\chi) \\[10pt]
    (2.1\mathrm{a}) \; \rightarrow \quad
        &\underset{\Delta\chi\rightarrow 0^{+}}{\sim}
            9 \, y \, \Delta\chi^2 \cdot \frac{\sigma_2}{3\,\Delta\chi^2} \\[10pt]
    y \rightarrow 1 \; \rightarrow \quad
        &= 9 \underbrace{y}_{\rightarrow 1} \frac{\cancel{\Delta\chi^2}}{3\cancel{\Delta\chi^2}}\,\sigma_2 \\[10pt]
        &= 3 \, \sigma_2 \; . \quad \quad (3.7)
\end{align*}
```

### Term 4: ``J_{22} I_2^2``

``I_2^2 \rightarrow \sigma_0/15`` is finite, while ``J_{22}`` carries an explicit
``\Delta\chi^{-4}``. We therefore need ``B_{22}``, Eq.(D.7), to order ``\Delta\chi^4`` — two orders
deeper than ``B_{00}`` and ``B_{02}``, which only needed ``\Delta\chi^2``.

The add-and-subtract trick used for ``B_{00}`` and ``B_{02}`` is impractical here, ``B_{22}``
being quartic in ``y``. But precisely because it is a **polynomial of degree 4**, its Taylor
expansion around ``y=1`` is finite and exact after five terms:

```math
\begin{align*}
    B_{22} &= 
         2(\chi_1^4 + \chi_2^4)(7 y^2 - 3) - 16 y \chi_1 \chi_2 (\chi_1^2 + \chi_2^2)(y^2 + 1)
        + \chi_1^2 \chi_2^2 (11y^4 + 14y^2 + 23) \\[10pt]
    &=\sum_{k=0}^{4} \frac{1}{k!}
        \frac{\partial^k B_{22}}{\partial y^k}\bigg|_{y=1} (y-1)^k \; .
\end{align*}
```

Before computing the five coefficients, one remark that makes all of them easy. ``B_{22}``
is **symmetric** under ``\chi_1 \leftrightarrow \chi_2``, so every coefficient can be written
in terms of the two elementary symmetric combinations

```math
    u := \chi_1^2 + \chi_2^2 \; , \qquad \qquad v := \chi_1 \chi_2 \; ,
    \quad \quad (\mathrm{D}.8)
```

for which

```math
    \chi_1^4 + \chi_2^4 = u^2 - 2v^2 \; , \qquad
    \chi_1\chi_2(\chi_1^2+\chi_2^2) = u\,v \; , \qquad
    \chi_1^2\chi_2^2 = v^2 \; ,
    \quad \quad (3.8\mathrm{a})
```

and, most importantly,

```math
    u - 2v = \chi_1^2 + \chi_2^2 - 2\chi_1\chi_2 = (\chi_1-\chi_2)^2 \; .
    \quad \quad (3.8\mathrm{b})
```

Every ``(u-2v)`` of the definitions (D.8) that appears is therefore a ``(\chi_1-\chi_2)^2``, i.e. an
``\mathcal{O}(\Delta\chi^2)``, and reading the order off a coefficient becomes a matter of
counting the powers of ``(u-2v)`` in it.

**The ``k=0`` coefficient.** Setting ``y=1`` in ``B_{22}``:

```math
\begin{align*}
    B_{22}\big|_{y=1} &= 2(\chi_1^4+\chi_2^4)(7-3) - 16\chi_1\chi_2(\chi_1^2+\chi_2^2)(1+1)
        + \chi_1^2\chi_2^2(11+14+23) \\[10pt]
        &= 8(\chi_1^4+\chi_2^4) - 32\chi_1\chi_2(\chi_1^2+\chi_2^2) + 48\chi_1^2\chi_2^2 \\[10pt]
    (\mathrm{D}.8),\;(3.8\mathrm{a}) \; \rightarrow \quad
        &= 8(u^2 - 2v^2) - 32uv + 48v^2 \\[10pt]
        &= 8u^2 - 16v^2 - 32uv + 48v^2 \\[10pt]
        &= 8\left(u^2 - 4uv + 4v^2\right) \\[10pt]
        &= 8(u-2v)^2 \\[10pt]
    (3.8\mathrm{b}) \; \rightarrow \quad
        &= 8(\chi_1-\chi_2)^4 \\[10pt]
        &\underset{\Delta\chi\rightarrow 0^{+}}{\rightarrow} \; \mathcal{O}(\Delta\chi^4) \; .
\end{align*}
```

(equivalently, expanding ``8[\chi_1^4 - 4\chi_1^3\chi_2 + 6\chi_1^2\chi_2^2 - 4\chi_1\chi_2^3 + \chi_2^4]`` and recognising the binomial coefficients of ``(\chi_1-\chi_2)^4``.)

**The ``k=1`` coefficient.** Differentiate ``B_{22}`` once with respect to ``y``, term by
term — only the ``y``-dependent factors move:

```math
\begin{align*}
    \frac{\partial B_{22}}{\partial y}
    &= 2(\chi_1^4+\chi_2^4)\frac{\partial (7y^2-3)}{\partial y}
        - 16\chi_1\chi_2(\chi_1^2+\chi_2^2)\frac{\partial \left[y(y^2+1)\right]}{\partial y}
        + \chi_1^2\chi_2^2 \frac{\partial (11y^4+14y^2+23)}{\partial y} \\[10pt]
    &= 2(\chi_1^4+\chi_2^4)(14y)
        - 16\chi_1\chi_2(\chi_1^2+\chi_2^2)(3y^2+1)
        + \chi_1^2\chi_2^2 (44y^3+28y) \\[10pt]
    &= 28y(\chi_1^4+\chi_2^4)
        - 16\chi_1\chi_2(\chi_1^2+\chi_2^2)(3y^2+1)
        + \chi_1^2\chi_2^2 (44y^3+28y) \; ,\\[15pt]
\end{align*}
```

```math
\begin{align*}
    \Rightarrow \quad  \frac{\partial B_{22}}{\partial y}\bigg|_{y=1}
    &= 28(\chi_1^4+\chi_2^4) - 16\chi_1\chi_2(\chi_1^2+\chi_2^2)(3+1)
        + \chi_1^2\chi_2^2 (44+28) \\[10pt]
    &= 28(\chi_1^4+\chi_2^4) - 64\chi_1\chi_2(\chi_1^2+\chi_2^2) + 72\chi_1^2\chi_2^2 \\[10pt]
    (3.8\mathrm{a}) \; \rightarrow \quad
    &= 28(u^2-2v^2) - 64uv + 72v^2 \\[10pt]
    &= 28u^2 - 56v^2 - 64uv + 72v^2 \\[10pt]
    &= 28u^2 - 64uv + 16v^2 \\[10pt]
    &= 4\left(7u^2 - 16uv + 4v^2\right) \; \\[10pt]
    &= 4(u-2v)(7u-2v) \\[10pt]
    (3.8\mathrm{b}) \; \rightarrow \quad
    &= 4(\chi_1-\chi_2)^2 \left(7\chi_1^2 - 2\chi_1\chi_2 + 7\chi_2^2\right)\\[10pt]
    &\underset{\Delta\chi\rightarrow 0^{+}}{\rightarrow} \; \mathcal{O}(\Delta\chi^2) \; .
\end{align*}
```

**The ``k=2`` coefficient.** Differentiating once more,

```math
\begin{align*}
    \frac{\partial^2 B_{22}}{\partial y^2}
    &= \frac{\partial}{\partial y}\left[\frac{\partial B_{22}}{\partial y}\right]\\[13pt]
    &= \frac{\partial}{\partial y}\left[
        28y(\chi_1^4+\chi_2^4)
        - 16\chi_1\chi_2(\chi_1^2+\chi_2^2)(3y^2+1)
        + \chi_1^2\chi_2^2 (44y^3+28y)
        \right]\\[13pt]
    &= 28(\chi_1^4+\chi_2^4) - 96 y \chi_1\chi_2(\chi_1^2+\chi_2^2)
        + \chi_1^2\chi_2^2 (132y^2+28) \\[10pt]
    \Rightarrow \quad
    \frac{1}{2}\frac{\partial^2 B_{22}}{\partial y^2}\bigg|_{y=1}
    &= 14(\chi_1^4+\chi_2^4) - 48\chi_1\chi_2(\chi_1^2+\chi_2^2) + 80\chi_1^2\chi_2^2 \\[10pt]
    (3.8\mathrm{a}) \; \rightarrow \quad
    &= 14(u^2-2v^2) - 48uv + 80v^2 \\[10pt]
    &= 14u^2 - 48uv + 52v^2 \; \\[10pt]
    &= \; 2\left(7u^2 - 24uv + 26v^2\right) \\[10pt]
    &\quad\quad \mathrm{This\; one \;does \; not \; contain \; a \; (u-2v) \; term!}\\[10pt]
    &\quad\quad \left[7u^2 - 24uv + 26v^2\right]\bigg|_{u = 2v} = 28v^2-48v^2+26v^2 = 6 v^2 \neq 0\\[10pt]
    (3.8\mathrm{b}) \; \rightarrow \quad
    &= 2\left(7\chi_1^4 - 24\chi_1^3\chi_2 + 40\chi_1^2\chi_2^2 - 24\chi_1\chi_2^3 + 7\chi_2^4\right)\\[10pt]
    &\underset{\Delta\chi\rightarrow 0^{+}}{\rightarrow} 2(7-24+40-24+7)\chi_1^4 \\[10pt]
    &= 12\chi_1^4 = \; \mathcal{O}(1) \; ,
\end{align*}
```


**The ``k=3`` and ``k=4`` coefficients.** Two more derivatives, the ``(\chi_1^4+\chi_2^4)``
term now being constant in ``y`` and dropping out:

```math
\begin{align*}
    \frac{\partial^3 B_{22}}{\partial y^3} 
        &= \frac{\partial}{\partial y}\left[\frac{\partial^2 B_{22}}{\partial y^2}\right]\\[15pt]
        &= \frac{\partial}{\partial y}\left[
            28(\chi_1^4+\chi_2^4) - 96 y \chi_1\chi_2(\chi_1^2+\chi_2^2)
            + \chi_1^2\chi_2^2 (132y^2+28) 
        \right]\\[15pt]
        &= - 96 \chi_1\chi_2(\chi_1^2+\chi_2^2) + 264 y \chi_1^2\chi_2^2\\[15pt]
    \Rightarrow \quad
    \frac{1}{6}\frac{\partial^3 B_{22}}{\partial y^3}\bigg|_{y=1}
        &= - 16 \chi_1\chi_2(\chi_1^2+\chi_2^2) + 44 \chi_1^2\chi_2^2 \\[16pt]
        &= -4\chi_1\chi_2\left(4\chi_1^2 - 11\chi_1\chi_2 + 4\chi_2^2\right) , \\[14pt]
        &\underset{\Delta\chi\rightarrow 0^{+}}{\rightarrow} -4(4-11+4)\chi_1^4 \\[10pt]
        &= 12\chi_1^4 = \; \mathcal{O}(1) \; \\[10pt]
\end{align*}
```

```math
\begin{align*}
    \frac{\partial^4 B_{22}}{\partial y^4}
        &= \frac{\partial}{\partial y}\left[\frac{\partial^3 B_{22}}{\partial y^3}\right]\\[15pt]
        &= \frac{\partial}{\partial y}\left[
            - 96 \chi_1\chi_2(\chi_1^2+\chi_2^2) + 264 y \chi_1^2\chi_2^2
        \right]\\[15pt]
        &= 264 \chi_1^2\chi_2^2\\[15pt]
    \quad \Rightarrow \quad
    \frac{1}{24}\frac{\partial^4 B_{22}}{\partial y^4}\bigg|_{y=1}
        &= 11 \chi_1^2\chi_2^2 , \\[14pt]
        &\underset{\Delta\chi\rightarrow 0^{+}}{\rightarrow} 11\chi_1^4 = \; \mathcal{O}(1) \; \\[10pt]
\end{align*}
```

both ``\mathcal{O}(1)``. Collecting the five coefficients:

```math
\begin{align*}
    B_{22}\big|_{y=1} &= 8(\chi_1-\chi_2)^4
        &&\rightarrow \; \mathcal{O}(\Delta\chi^4) \; , \\[6pt]
    \frac{\partial B_{22}}{\partial y}\bigg|_{y=1}
        &= 4(\chi_1-\chi_2)^2\left(7\chi_1^2 - 2\chi_1\chi_2 + 7\chi_2^2\right)
        &&\rightarrow \; \mathcal{O}(\Delta\chi^2) \; , \\[6pt]
    \frac{1}{2}\frac{\partial^2 B_{22}}{\partial y^2}\bigg|_{y=1}
        &= 2\left(7\chi_1^4 - 24\chi_1^3\chi_2 + 40\chi_1^2\chi_2^2 - 24\chi_1\chi_2^3
            + 7\chi_2^4\right)
        &&\rightarrow \; \mathcal{O}(1) \; , \\[6pt]
    \frac{1}{6}\frac{\partial^3 B_{22}}{\partial y^3}\bigg|_{y=1}
        &= -4\chi_1\chi_2\left(4\chi_1^2 - 11\chi_1\chi_2 + 4\chi_2^2\right)
        &&\rightarrow \; \mathcal{O}(1) \; , \\[6pt]
    \frac{1}{24}\frac{\partial^4 B_{22}}{\partial y^4}\bigg|_{y=1}
        &= 11 \chi_1^2\chi_2^2
        &&\rightarrow \; \mathcal{O}(1) \; ,
\end{align*}
```

so that, exactly,

```math
\begin{align*}
    B_{22} &=
        \sum_{k=0}^{4} \frac{1}{k!}
        \frac{\partial^k B_{22}}{\partial y^k}\bigg|_{y=1} (y-1)^k \; \\[15pt]
    &= 8(\chi_1-\chi_2)^4
        + 4(\chi_1-\chi_2)^2\left(7\chi_1^2 - 2\chi_1\chi_2 + 7\chi_2^2\right)(y-1) \\
        &\quad+ 2\left(7\chi_1^4 - 24\chi_1^3\chi_2 + 40\chi_1^2\chi_2^2 - 24\chi_1\chi_2^3
            + 7\chi_2^4\right)(y-1)^2 \\
        &\quad- 4\chi_1\chi_2\left(4\chi_1^2 - 11\chi_1\chi_2 + 4\chi_2^2\right)(y-1)^3
        + 11 \chi_1^2\chi_2^2 (y-1)^4 \; . \quad \quad (3.8)
\end{align*}
```

Now count the orders, remembering (D.4) and (2.4), ``(\chi_1-\chi_2) = -p\Delta\chi`` and
``(y-1) = \mathcal{O}(\Delta\chi^2)``:

```math
\begin{align*}
    \mathrm{1st \; term} \; &: \quad
        B_{22}\big|_{y=1} \cdot 1 &&=
        \mathcal{O}(\Delta\chi^4) \cdot 1 
        &&= \mathcal{O}(\Delta\chi^4)
        &&\Longrightarrow \; \mathrm{keep} \; , \\[6pt]
    \mathrm{2nd \; term} \; &: \quad
        \frac{\partial B_{22}}{\partial y}\bigg|_{y=1}(y-1) &&=
        \mathcal{O}(\Delta\chi^2) \cdot \mathcal{O}(\Delta\chi^2) 
        &&= \mathcal{O}(\Delta\chi^4)
        &&\Longrightarrow \; \mathrm{keep} \; , \\[6pt]
    \mathrm{3rd \; term} \; &: \quad
        \frac{1}{2}\frac{\partial^2 B_{22}}{\partial y^2}\bigg|_{y=1} (y-1)^2 &&= 
        \mathcal{O}(1) \cdot \mathcal{O}(\Delta\chi^4) 
        &&= \mathcal{O}(\Delta\chi^4)
        &&\Longrightarrow \; \mathrm{keep} \; , \\[6pt]
    \mathrm{4th \; term} \; &: \quad
        \frac{1}{6}\frac{\partial^3 B_{22}}{\partial y^3}\bigg|_{y=1}(y-1)^3 &&=
        \mathcal{O}(1) \cdot \mathcal{O}(\Delta\chi^6) 
        &&= \mathcal{O}(\Delta\chi^6)
        &&\Longrightarrow \; \mathrm{drop} \; , \\[6pt]
    \mathrm{5th \; term} \; &: \quad
        \frac{1}{24}\frac{\partial^4 B_{22}}{\partial y^4}\bigg|_{y=1}(y-1)^4 &&=
        \mathcal{O}(1) \cdot \mathcal{O}(\Delta\chi^8) 
        &&= \mathcal{O}(\Delta\chi^8)
        &&\Longrightarrow \; \mathrm{drop} \; .
\end{align*}
```

Three terms survive. Substituting (D.4) and (2.3a), and evaluating the regular coefficients
at ``\chi_1 = \chi_2`` (where ``7-2+7 = 12`` and ``7-24+40-24+7 = 6``):

```math
\begin{align*}
(D.4&) : \quad \chi_2 := \chi_1 + p \, \Delta\chi\\[10pt]
(2.3\mathrm{a}&): \quad y-1  \underset{\Delta\chi\rightarrow 0^{+}}{\sim}  - \frac{(1-p^2)}{2\chi_1^2}\Delta\chi^2\\[15pt]
    \mathrm{1st} \; : \left[B_{22}\big|_{y=1} \cdot 1 \right]\bigg|_{(D.4)} 
        &= \left[8(\chi_1-\chi_2)^4\right]\bigg|_{(D.4)}
        = \quad  8(-p\Delta\chi)^4  \\
        &= 8 \, p^4 \, \Delta\chi^4 \; , \\[8pt]
    \mathrm{2nd} \; : \left[\frac{\partial B_{22}}{\partial y}\bigg|_{y=1}(y-1)\right]\bigg|_{(D.4)} 
        &= \left[4(\chi_1-\chi_2)^2\left(7\chi_1^2 - 2\chi_1\chi_2 + 7\chi_2^2\right)(y-1)\right]\bigg|_{(D.4)}\\
        &= \quad  4 (-p\Delta\chi)^2 \cdot (7-2+7)\cancel{\chi_1^2} \cdot
        \left[- \frac{(1-p^2)}{2\cancel{\chi_1^2}}\Delta\chi^2\right]\\
        &= -24 \, p^2(1-p^2) \, \Delta\chi^4 \; , \\[8pt]
    \mathrm{3rd} \; :
        \left[\frac{1}{2}\frac{\partial^2 B_{22}}{\partial y^2}\bigg|_{y=1} (y-1)^2 \right]\bigg|_{(D.4)}
        &=\left[
            2\left(7\chi_1^4 - 24\chi_1^3\chi_2 + 40\chi_1^2\chi_2^2 - 24\chi_1\chi_2^3
            + 7\chi_2^4\right)(y-1)^2
        \right]\bigg|_{(D.4)} \\
        &= \quad  2 \, (7-24+40-24+7) \, \chi_1^4 \cdot
        \left[- \frac{(1-p^2)}{2\chi_1^2}\Delta\chi^2\right]^{2}\\
        &= 2 \cdot 6 \, \cancel{\chi_1^4} \, \frac{(1-p^2)^2}{4\cancel{\chi_1^4}}\Delta\chi^4\\
        &= 3 \, (1-p^2)^2 \, \Delta\chi^4 \; ,
\end{align*}
```

and summing:

```math
\begin{align*}
    B_{22} &= \sum_{k=0}^{4} \frac{1}{k!}
        \frac{\partial^k B_{22}}{\partial y^k}\bigg|_{y=1} (y-1)^k \; \\[15pt]
    &\underset{\Delta\chi\rightarrow 0^{+}}{\sim} 
        \left[B_{22}\big|_{y=1} \cdot 1 \right]\bigg|_{(D.4)} +
        \left[\frac{\partial B_{22}}{\partial y}\bigg|_{y=1}(y-1)\right]\bigg|_{(D.4)} +
        \left[\frac{1}{2}\frac{\partial^2 B_{22}}{\partial y^2}\bigg|_{y=1} (y-1)^2 \right]\bigg|_{(D.4)}
    \\[15pt]
    &= \left[8 \, p^4 \, \Delta\chi^4 \right] 
        + \left[- 24 \, p^2(1-p^2) \, \Delta\chi^4 \right]
        + \left[3 \, (1-p^2)^2 \, \Delta\chi^4 \right] \\[10pt]
    &=    \left[ 8p^4 - 24p^2(1-p^2) + 3(1-p^2)^2 \right] \Delta\chi^4 \\[10pt]
    &= \left[ 8p^4 - 24p^2 + 24p^4 + 3 - 6p^2 + 3p^4 \right] \Delta\chi^4 \\[10pt]
    &= \left( 35p^4 - 30p^2 + 3 \right) \Delta\chi^4 \\[10pt]
    &= \; 8 \, \mathcal{L}_4(p) \, \Delta\chi^4 \; , \quad \quad (3.9)
\end{align*}
```

where ``\mathcal{L}_4`` is  the fourth Legendre polynomial — not a coincidence, since the whole
construction is an expansion in ``y = \cos\theta``. Finally:

```math
\begin{align*}
    J_{22} I_2^2 &= \frac{9}{4}\frac{\chi_1\chi_2}{\Delta\chi^4} \, B_{22} \, I_2^2 \\[10pt]
    (2.1\mathrm{a})\, , \; (3.9) \; \rightarrow \quad
        &\underset{\Delta\chi\rightarrow 0^{+}}{\sim}
        \frac{9}{4}\frac{\chi_1^2}{\cancel{\Delta\chi^4}}
        \left( 35p^4 - 30p^2 + 3 \right) \cancel{\Delta\chi^4} \cdot \frac{\sigma_0}{15} \\[10pt]
    &= \frac{3}{20}\,\chi_1^2 \, \sigma_0 \left( 35p^4 - 30p^2 + 3 \right) \; . \quad \quad (3.10)
\end{align*}
```

### The sum


```math
\begin{align*}
(3.4)  : \quad J_{00} I_0^0 &\underset{\Delta\chi\rightarrow 0^{+}}{\sim} 
    \dfrac{3}{4}\chi_1^2\sigma_0 \left(-7p^4 + 6p^2 + 1\right) \\[13pt]
(3.6)  : \quad J_{02} I_2^0 &\underset{\Delta\chi\rightarrow 0^{+}}{\sim} 
    0 \\[13pt]
(3.7)  : \quad J_{31} I_1^3 &\underset{\Delta\chi\rightarrow 0^{+}}{\sim} 
    3\,\sigma_2 \\[13pt]
(3.10) : \quad J_{22} I_2^2 &\underset{\Delta\chi\rightarrow 0^{+}}{\sim} 
    \dfrac{3}{20}\chi_1^2\sigma_0 \left(35p^4 - 30p^2 + 3\right) \\[13pt]
\end{align*}
```

Only two terms carry a ``p``, and it must cancel between them:

```math
\begin{align*}
    J^{\kappa\kappa}_{00} I_0^0 + 
    J^{\kappa\kappa}_{02} I_2^0 +
    J^{\kappa\kappa}_{31} I_1^3 +
    J^{\kappa\kappa}_{22} I_2^2 
    &\underset{\Delta\chi\rightarrow 0^{+}}{\sim}
        \dfrac{3}{4}\chi_1^2\sigma_0 \left(-7p^4 + 6p^2 + 1\right) +
        0 + 3\,\sigma_2 + \dfrac{3}{20}\chi_1^2\sigma_0 \left(35p^4 - 30p^2 + 3\right)
    \\[10pt]
    &= 3\, \sigma_2+\frac{3}{20}\left[ 
            5\left(-7p^4+6p^2+1\right) + \left(35p^4-30p^2+3\right) 
        \right]\chi_1^2 \sigma_0 \\[10pt]
    &= 3\, \sigma_2+ \frac{3}{20}\left[ 
            \cancel{-35p^4}+\cancel{30p^2}+5 + \cancel{35p^4}-\cancel{30p^2}+3 
        \right]\chi_1^2 \sigma_0 \\[10pt]
    &= 3\, \sigma_2 + \frac{6}{5}\chi_1^2 \sigma_0  \; .
\end{align*}
```

The cancellation of ``p`` is the proof that the limit exists and does not depend on the
direction of approach.

```math
    \boxed{\;
    \lim_{\Delta\chi \rightarrow 0^{+}}
    \left(
        J_{00}^{\kappa\kappa}I_0^0 + J_{02}^{\kappa\kappa}I_2^0 + 
        J_{31}^{\kappa\kappa}I_1^3 + J_{22}^{\kappa\kappa}I_2^2
    \right)
    = 3\,\sigma_2 + \frac{6}{5}\,\chi_1^2\,\sigma_0 \; . }
    \quad \quad (3.11)
```



!!! warning "A common mistake"
    If in Term 1 one sets ``\chi_1 = \chi_2`` inside ``B_{00}`` before expanding, the
    ``8(\chi_1-\chi_2)^2 = 8p^2\Delta\chi^2`` term of Eq.(3.2) is lost. One is left with
    ``B_{00}\big|_{\chi_1=\chi_2} = \chi_1^2(16y-9y^2-7) = -\chi_1^2(y-1)(9y-7)``, hence
    ``B_{00} \sim 2\chi_1^2(1-y) \sim (1-p^2)\Delta\chi^2`` instead of Eq.(3.3)'s
    ``(1+7p^2)\Delta\chi^2``, and therefore

    ```math
        J_{00}I_0^0 = \frac{3}{4}(1-p^2)^2\chi_1^2\sigma_0
        \qquad \mathrm{instead \; of \; Eq.(3.4)} \; .
    ```

    That result is wrong: summed with Eq.(3.10) it leaves a residual ``p``-dependence,

    ```math
        \frac{3}{4}(1-p^2)^2 + \frac{3}{20}\left(35p^4-30p^2+3\right)
        = 6p^4 - 6p^2 + \frac{6}{5}
    ```

    instead of the constant ``6/5`` of Eq.(3.11), which would mean that the limit does not
    exist. Note that this wrong expression happens to give the right value ``6/5`` at both
    ``p = 0`` and ``p = \pm 1``, and is off by a factor ``-1/4`` at ``p^2 = 1/2``: spot
    checks at the endpoints would not catch it. Only the full ``p``-cancellation does.
