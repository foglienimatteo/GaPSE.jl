# Family 2: Lensing ``\times`` Doppler

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

Concerned functions: `integrand_ξ_GNC_Lensing_Doppler`,
`integrand_ξ_GNCxLD_Lensing_Doppler`, `integrand_ξ_GNCxLD_Doppler_Lensing`,
`integrand_ξ_LD_Lensing_Doppler`. The sum to be taken to the limit is

```math
    J_{00}I_0^0 + J_{02}I_2^0 + J_{04}I_4^0 + J_{20}I_0^2 \; ,
```

with

```math
\begin{align*}
    J_{00} &= \frac{1}{15}\left(\chi_1^2 y + \chi_1 s_2 (4y^2-3) - 2 y s_2^2\right) \; , \\[6pt]
    J_{02} &= \frac{1}{42\,\Delta\chi_1^2}\underbrace{\left(
        4\chi_1^4 y + 4\chi_1^3 (2y^2-3) s_2 + \chi_1^2 y (11 - 23y^2) s_2^2
        + \chi_1 (23y^2-3) s_2^3 - 8 y s_2^4 \right)}_{=: \; N_{02} \; , \quad (\mathrm{D}.9)} \; , \\[6pt]
    J_{04} &= \frac{1}{70\,\Delta\chi_1^2}\underbrace{\left(
        2\chi_1^4 y + 2\chi_1^3 (2y^2-3) s_2 - \chi_1^2 y (y^2+5) s_2^2
        + \chi_1 (y^2+9) s_2^3 - 4 y s_2^4 \right)}_{=: \; N_{04} \; , \quad (\mathrm{D}.10)} \; , \\[6pt]
    J_{20} &= y \, \Delta\chi_1^2 \; .
\end{align*}
```

### Term 1: ``J_{00} I_0^0``

``J_{00}`` is **regular** — it carries no negative power of ``\Delta\chi_1`` — so, since
``I_0^0 \rightarrow \sigma_0`` is finite as well, all we need is the value of ``J_{00}`` at
the singular point. It vanishes there, so we decompose it exactly, as in Family 1:

```math
\begin{align*}
    15 \, J_{00} &= \chi_1^2 y + \chi_1 s_2 (4y^2-3) - 2 y s_2^2 \\[10pt]
    y = 1 + (y-1) \; \rightarrow \quad
        &= \chi_1^2 + \chi_1^2 (y-1) + \chi_1 s_2 (4y^2-3) - 2 s_2^2 - 2 s_2^2 (y-1) \\[10pt]
    4y^2-3 = 1 + 4(y^2-1) \; \rightarrow \quad
        &= \chi_1^2 + \chi_1 s_2 - 2 s_2^2
            + (y-1)\left(\chi_1^2 - 2 s_2^2\right) + 4(y^2-1)\chi_1 s_2 \\[10pt]
    \chi_1^2 + \chi_1 s_2 - 2 s_2^2 = (\chi_1 - s_2)(\chi_1 + 2 s_2) \; \rightarrow \quad
        &= (\chi_1 - s_2)(\chi_1 + 2 s_2)
            + (y-1)\left(\chi_1^2 - 2 s_2^2\right) + 4(y^2-1)\chi_1 s_2 \; .
    \quad \quad (4.1)
\end{align*}
```

Now the crucial difference with Family 1: here the ``(\chi_1 - s_2)`` factor appears
**linearly**, not quadratically. By (2.4) the three terms are therefore

```math
    \underbrace{(\chi_1 - s_2)(\chi_1 + 2 s_2)}_{\mathcal{O}(\Delta\chi_1)}
    \; + \;
    \underbrace{(y-1)\left(\chi_1^2 - 2 s_2^2\right)}_{\mathcal{O}(\Delta\chi_1^2)}
    \; + \;
    \underbrace{4(y^2-1)\chi_1 s_2}_{\mathcal{O}(\Delta\chi_1^2)} \; ,
```

i.e. the first one **dominates** and the other two are subleading — the opposite of what
happened to ``B_{00}``, where the ``(\chi_1-\chi_2)^2`` was of the same order as the rest.
Hence, with ``(\chi_1 - s_2) = -p\Delta\chi_1`` from (D.4) and ``s_2 \rightarrow \chi_1``:

```math
\begin{align*}
    15 \, J_{00} &\underset{\Delta\chi_1\rightarrow 0^{+}}{\sim}
        (-p\Delta\chi_1)(\chi_1 + 2\chi_1) = -3 \, p \, \chi_1 \, \Delta\chi_1 \\[10pt]
    \Rightarrow \quad
    J_{00} I_0^0 &\underset{\Delta\chi_1\rightarrow 0^{+}}{\sim}
        -\frac{3 \, p \, \chi_1 \, \Delta\chi_1}{15}\,\sigma_0
        = -\frac{p \, \chi_1 \, \sigma_0}{5}\,\Delta\chi_1
    \; \xrightarrow[\Delta\chi_1 \rightarrow 0^{+}]{} \; 0 \; .
    \quad \quad (4.2)
\end{align*}
```

### Term 2: ``J_{02} I_2^0`` and Term 3: ``J_{04} I_4^0``

Both ``J`` carry an explicit ``\Delta\chi_1^{-2}``, so what matters is the order at which
their numerators ``N_{02}`` and ``N_{04}``, Eqs.(D.9) and (D.10), vanish. As for ``B_{22}`` in Family 1, they are
polynomials in ``y`` — of degree 2 and 3 respectively in each term, degree 3 overall — so
their Taylor expansion around ``y=1`` terminates:

```math
    N = \sum_{k} \frac{1}{k!}\frac{\partial^k N}{\partial y^k}\bigg|_{y=1} (y-1)^k \; .
```

**The ``k=0`` coefficients.** Setting ``y = 1``:

```math
\begin{align*}
    N_{02}\big|_{y=1} &= 4\chi_1^4 + 4\chi_1^3 (2-3) s_2 + \chi_1^2 (11-23) s_2^2
        + \chi_1 (23-3) s_2^3 - 8 s_2^4 \\[10pt]
        &= 4\chi_1^4 - 4\chi_1^3 s_2 - 12\chi_1^2 s_2^2 + 20\chi_1 s_2^3 - 8 s_2^4 \\[10pt]
        &= 4\left(\chi_1^4 - \chi_1^3 s_2 - 3\chi_1^2 s_2^2 + 5\chi_1 s_2^3 - 2 s_2^4\right)
        \\[10pt]
        &= 4(\chi_1-s_2)^3(\chi_1+2s_2) \; , \\[16pt]
    N_{04}\big|_{y=1} &= 2\chi_1^4 + 2\chi_1^3 (2-3) s_2 - \chi_1^2 (1+5) s_2^2
        + \chi_1 (1+9) s_2^3 - 4 s_2^4 \\[10pt]
        &= 2\chi_1^4 - 2\chi_1^3 s_2 - 6\chi_1^2 s_2^2 + 10\chi_1 s_2^3 - 4 s_2^4 \\[10pt]
        &= 2(\chi_1-s_2)^3(\chi_1+2s_2) \; = \; \frac{1}{2} N_{02}\big|_{y=1} \; ,
\end{align*}
```

where the factorisation is checked by expanding back,

```math
    (\chi_1-s_2)^3(\chi_1+2s_2)
    = \left(\chi_1^3 - 3\chi_1^2 s_2 + 3\chi_1 s_2^2 - s_2^3\right)(\chi_1+2s_2)
    = \chi_1^4 - \chi_1^3 s_2 - 3\chi_1^2 s_2^2 + 5\chi_1 s_2^3 - 2 s_2^4 \; .
```

Both therefore vanish as ``(\chi_1-s_2)^3 = \mathcal{O}(\Delta\chi_1^3)`` — **cubically**, not
just at the point.

**The ``k=1`` coefficients.** Differentiating once in ``y``,

```math
\begin{align*}
    \frac{\partial N_{02}}{\partial y} &= 4\chi_1^4 + 16\chi_1^3 y \, s_2
        + \chi_1^2 (11 - 69y^2) s_2^2 + 46\chi_1 y \, s_2^3 - 8 s_2^4 \; , \\[6pt]
    \frac{\partial N_{04}}{\partial y} &= 2\chi_1^4 + 8\chi_1^3 y \, s_2
        - \chi_1^2 (3y^2+5) s_2^2 + 2\chi_1 y \, s_2^3 - 4 s_2^4 \; ,
\end{align*}
```

and evaluating at ``y=1``,

```math
\begin{align*}
    \frac{\partial N_{02}}{\partial y}\bigg|_{y=1}
        &= 4\chi_1^4 + 16\chi_1^3 s_2 - 58\chi_1^2 s_2^2 + 46\chi_1 s_2^3 - 8 s_2^4 \\[10pt]
        &= 2\left(2\chi_1^4 + 8\chi_1^3 s_2 - 29\chi_1^2 s_2^2 + 23\chi_1 s_2^3 - 4 s_2^4\right) \\[10pt]
        &= 2(\chi_1-s_2)\left(2\chi_1^3 + 10\chi_1^2 s_2 - 19\chi_1 s_2^2 + 4 s_2^3\right)
        \; , \\[16pt]
    \frac{\partial N_{04}}{\partial y}\bigg|_{y=1}
        &= 2\chi_1^4 + 8\chi_1^3 s_2 - 8\chi_1^2 s_2^2 + 2\chi_1 s_2^3 - 4 s_2^4 \\[10pt]
        &= 2\left(\chi_1^4 + 4\chi_1^3 s_2 - 4\chi_1^2 s_2^2 + \chi_1 s_2^3 - 2 s_2^4\right) \\[10pt]
        &= 2(\chi_1-s_2)\left(\chi_1^3 + 5\chi_1^2 s_2 + \chi_1 s_2^2 + 2 s_2^3\right) \; ,
\end{align*}
```

both containing one power of ``(\chi_1-s_2)``, i.e. ``\mathcal{O}(\Delta\chi_1)``, so that —
multiplied by ``(y-1) = \mathcal{O}(\Delta\chi_1^2)`` — they too contribute at
``\mathcal{O}(\Delta\chi_1^3)``.

**The ``k \geq 2`` coefficients.** These are ``\mathcal{O}(1)``,

```math
    \frac{1}{2}\frac{\partial^2 N_{02}}{\partial y^2}\bigg|_{y=1}
        = \chi_1 s_2\left(8\chi_1^2 - 69\chi_1 s_2 + 23 s_2^2\right) \; , \qquad
    \frac{1}{6}\frac{\partial^3 N_{02}}{\partial y^3}\bigg|_{y=1} = -23\chi_1^2 s_2^2 \; ,
```

and multiply ``(y-1)^2 = \mathcal{O}(\Delta\chi_1^4)`` or higher, so they are **dropped**
(the same holds for ``N_{04}``).

**Putting the two surviving orders together.** With ``(\chi_1 - s_2) = -p\Delta\chi_1`` from
(D.4), ``(y-1)`` from (2.3a), and the cubics evaluated at ``s_2 = \chi_1``
(``2+10-19+4 = -3`` and ``1+5+1+2 = 9``):

```math
\begin{align*}
    N_{02} &\underset{\Delta\chi_1\rightarrow 0^{+}}{\sim}
        \underbrace{4(-p\Delta\chi_1)^3(3\chi_1)}_{k=0}
        + \underbrace{2(-p\Delta\chi_1)(-3\chi_1^3)
            \left[- \frac{(1-p^2)}{2\cancel{\chi_1^2}}\Delta\chi_1^2\right]}_{k=1} \\[10pt]
    &= -12\,\chi_1\, p^3 \Delta\chi_1^3 - 3\,\chi_1\, p \,(1-p^2)\,\Delta\chi_1^3 \\[10pt]
    &= -3\,\chi_1\, p \left[4p^2 + (1-p^2)\right]\Delta\chi_1^3 \\[10pt]
    &= -3\,\chi_1\, p \,(3p^2+1)\,\Delta\chi_1^3 \; , \\[20pt]
    N_{04} &\underset{\Delta\chi_1\rightarrow 0^{+}}{\sim}
        \underbrace{2(-p\Delta\chi_1)^3(3\chi_1)}_{k=0}
        + \underbrace{2(-p\Delta\chi_1)(9\chi_1^3)
            \left[- \frac{(1-p^2)}{2\cancel{\chi_1^2}}\Delta\chi_1^2\right]}_{k=1} \\[10pt]
    &= -6\,\chi_1\, p^3 \Delta\chi_1^3 + 9\,\chi_1\, p \,(1-p^2)\,\Delta\chi_1^3 \\[10pt]
    &= 3\,\chi_1\, p \left[-2p^2 + 3(1-p^2)\right]\Delta\chi_1^3 \\[10pt]
    &= 3\,\chi_1\, p \,(3-5p^2)\,\Delta\chi_1^3 \; .
\end{align*}
```

Dividing by the explicit ``\Delta\chi_1^{2}`` leaves both ``J`` of order
``\mathcal{O}(\Delta\chi_1)``, and the ``I_\ell^n`` they multiply vanish as well:

```math
\begin{align*}
    J_{02} I_2^0 &\underset{\Delta\chi_1\rightarrow 0^{+}}{\sim}
        \underbrace{-\frac{3\chi_1 p(3p^2+1)}{42}\Delta\chi_1}_{\mathcal{O}(\Delta\chi_1)}
        \cdot \underbrace{\frac{\sigma_{-2}}{15}\Delta\chi_1^2}_{\mathcal{O}(\Delta\chi_1^2)}
        = \mathcal{O}(\Delta\chi_1^3) \; \xrightarrow[\Delta\chi_1 \rightarrow 0^{+}]{} \; 0 \; , \\[10pt]
    J_{04} I_4^0 &\underset{\Delta\chi_1\rightarrow 0^{+}}{\sim}
        \underbrace{\frac{3\chi_1 p(3-5p^2)}{70}\Delta\chi_1}_{\mathcal{O}(\Delta\chi_1)}
        \cdot \underbrace{\frac{\sigma_{-4}}{945}\Delta\chi_1^4}_{\mathcal{O}(\Delta\chi_1^4)}
        = \mathcal{O}(\Delta\chi_1^5) \; \xrightarrow[\Delta\chi_1 \rightarrow 0^{+}]{} \; 0 \; .
    \quad \quad (4.3)
\end{align*}
```

### Term 4: ``J_{20} I_0^2``

The only surviving term, and the simplest: ``J_{20} = y\,\Delta\chi_1^2`` is regular and its
``\Delta\chi_1^2`` cancels the ``\Delta\chi_1^{-2}`` of ``I_0^2``.

```math
\begin{align*}
    J_{20} I_0^2 &= y \, \Delta\chi_1^2 \cdot I_0^2(\Delta\chi_1) \\[10pt]
    (2.1\mathrm{a}) \; \rightarrow \quad
        &\underset{\Delta\chi_1\rightarrow 0^{+}}{\sim}
        \underbrace{y}_{\rightarrow 1}\cancel{\Delta\chi_1^2} \,
        \frac{\sigma_2}{\cancel{\Delta\chi_1^2}} = \sigma_2 \; .
\end{align*}
```

### The sum

| term | order | limit |
|:--|:--|:--|
| ``J_{00} I_0^0`` , Eq.(4.2) | ``\mathcal{O}(\Delta\chi_1)`` | ``0`` |
| ``J_{02} I_2^0`` , Eq.(4.3) | ``\mathcal{O}(\Delta\chi_1^3)`` | ``0`` |
| ``J_{04} I_4^0`` , Eq.(4.3) | ``\mathcal{O}(\Delta\chi_1^5)`` | ``0`` |
| ``J_{20} I_0^2`` | ``\mathcal{O}(1)`` | ``\sigma_2`` |

```math
    \boxed{\;
    \lim_{\Delta\chi_1 \rightarrow 0}
    \left(J_{00}I_0^0 + J_{02}I_2^0 + J_{04}I_4^0 + J_{20}I_0^2\right) = \sigma_2 \; . }
    \quad \quad (4.4)
```

No ``p`` survives in any of the four limits, so there is nothing to cancel: here the
direction-independence is immediate.
