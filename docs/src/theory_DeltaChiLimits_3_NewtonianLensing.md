# Family 3: Newtonian ``\times`` Lensing

## Recap: everything this page needs

This page is self-contained. All the results quoted here are derived in
[The ``\Delta\chi \rightarrow 0`` limits](theory_DeltaChiLimits.md), and the ``I_\ell^n``
asymptotics in [The ``I_\ell^n`` integrals](theory_IlnIntegrals.md).

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

Concerned functions: `integrand_ξ_GNC_Newtonian_Lensing`,
`integrand_ξ_GNCxLD_Newtonian_Lensing`.
The GNC Newtonian-Lensing function is:

```math
\begin{split}
    \xi^{\delta \kappa} ( s_1 , s_2, y ) =
    D_1  \int_0^{s_2}\mathrm{d} \chi_2
    J^{\delta \kappa}_{\alpha}
    \left[ 
        J^{\delta \kappa}_{00} I_0^0 ( \Delta\chi_2 ) + 
        J^{\delta \kappa}_{02} I_2^0 ( \Delta \chi_2 ) + 
        J^{\delta \kappa}_{04} I_4^0 ( \Delta \chi_2 ) 
    \right] \, ,
\end{split}
```

with (``f_1 = f(s_1)`` the growth rate and ``b_1 = b(s_1)`` the bias)

```math
\begin{align*}
    J^{\delta \kappa}_{\alpha} &=
    \frac{
        \mathcal{H}_0 ^2 \Omega_{\mathrm{M}0} D (\chi_2)
    }{
        a(\chi_2 ) s_2
    } 
    (\chi_2 - s_2 ) (5s_{\mathrm{b}, 2} - 2) 
    \, , \\[6pt]
    J^{\delta \kappa}_{00} &=
        \frac{1}{5}\left[ f_1 \chi_2 (3 y^2-1) - 3 y s_1 f_1 - 5 y s_1 b_1 \right] \; , \\[6pt]
    J^{\delta \kappa}_{02} &=
        \frac{1}{14\,\Delta\chi_2^2}\left[
        7 s_1 b_1 \underbrace{\left(-2\chi_2^2 y + \chi_2 s_1 (y^2+3) - 2 y s_1^2\right)}_{=: \; N_{02}^{(b)} \; , \; (\mathrm{D}.11)}
        \right. \\
        &\phantom{= \frac{1}{14\,\Delta\chi_2^2}\Big[} \left.
        + f_1 \underbrace{\left(
            4\chi_2^3 (3y^2-1) - 2\chi_2^2 y s_1 (3y^2+8)
            + \chi_2 s_1^2 (9y^2+11) - 6 y s_1^3
        \right)}_{=: \; N_{02}^{(f)} \; , \; (\mathrm{D}.12)} \right] \; , \\[6pt]
    J^{\delta \kappa}_{04} &=
        \frac{3 \, f_1}{70\,\Delta\chi_2^4} \underbrace{\left[
        \chi_2^5 (6y^2-2) + 6\chi_2^4 y s_1 (y^2-3) - \chi_2^3 s_1^2 (y^4+12y^2-21)
        + 2\chi_2^2 y s_1^3 (y^2+3) - 12 \chi_2 s_1^4 + 4 y s_1^5 \right]}_{=: \; N_{04} \; , \; (\mathrm{D}.13)} \; .
\end{align*}
```

The limit to be taken is:

```math
    \lim_{\Delta\chi_2\rightarrow 0^{+}}
    \left(J_{00}I_0^0 + J_{02}I_2^0 + J_{04}I_4^0\right) \; .
```

This analysis is valid for Newtonian-Lensing in both the combinations (GNC and GNCxLD).

### Term 1: ``J_{00} I_0^0``

This family is the exception: ``J_{00}`` is regular **and does not vanish** at the singular
point, so it is the only one whose surviving contribution comes from an ``I_0^0`` rather than
from an ``I_0^2`` or an ``I_1^3``. Precisely because nothing cancels, here we *may* set
``y = 1`` and ``\chi_2 = s_1`` straight away — the trap of (2.4) only bites when a bracket
vanishes:

```math
\begin{align*}
    J_{00}^{\delta \kappa} &= \frac{1}{5}\left[ f_1 \chi_2 (3y^2-1) - 3 y s_1 f_1 - 5 y s_1 b_1 \right] \\[10pt]
    &\quad\quad y = 1 \; , \; \chi_2 = s_1 \\[10pt]
        &\xrightarrow[\Delta\chi_2\rightarrow 0^{+}]{}
            \frac{1}{5}\left[ f_1 s_1 (3-1) - 3 s_1 f_1 - 5 s_1 b_1 \right] \\[10pt]
        &= \frac{1}{5}\left[ 2 f_1 s_1 - 3 f_1 s_1 - 5 b_1 s_1 \right] \\[10pt]
        &= \frac{1}{5}\left[ (2-3) f_1 s_1 - 5 b_1 s_1 \right] \\[10pt]
        &= \frac{1}{5}\left[ - f_1 s_1 - 5 b_1 s_1 \right] \\[10pt]
        &= -\frac{s_1 \, (f_1 + 5 b_1)}{5} \; ,
    \quad \quad (5.1)
\end{align*}
```

and therefore, with ``I_0^0 \rightarrow \sigma_0`` from (2.1a),

```math
\begin{align*}
    J_{00}^{\delta \kappa} I_0^0
    &\xrightarrow[\Delta\chi_2\rightarrow 0^{+}]{}
        -\frac{s_1 \, (f_1 + 5 b_1) \, \sigma_0}{5} \; ,\\[10pt]
    &\underset{\Delta\chi_2 \rightarrow 0^{+}}{\sim}\mathcal{O}(1)
    \quad \quad (5.2)
\end{align*}
```

Being the value of a regular function at a point, this carries no ``p``: the
direction-independence is automatic for this term.

### Term 2: ``J_{02} I_2^0``

``J_{02}`` carries ``\Delta\chi_2^{-2}`` and ``I_2^0`` carries ``\Delta\chi_2^{+2}``, so the
two exactly compensate and the product is **finite**: what decides is the order at which the
two numerators ``N_{02}^{(b)}`` and ``N_{02}^{(f)}``, Eqs.(D.11) and (D.12), vanish. Both are
polynomials in ``y`` — of degree 2 and 3 — so their Taylor expansion around ``y=1``
terminates:

```math
    N = \sum_{k} \frac{1}{k!}\frac{\partial^k N}{\partial y^k}\bigg|_{y=1} (y-1)^k \; .
```

**The ``k=0`` coefficients.** Setting ``y = 1``:

```math
\begin{align*}
(\mathrm{D}.11) : \quad N_{02}^{(b)} &:= \left[-2\chi_2^2 y + \chi_2 s_1 (y^2+3) - 2 y s_1^2\right] \\[10pt]
    \Rightarrow \quad
    N_{02}^{(b)}\big|_{y=1} &= -2\chi_2^2 + \chi_2 s_1 (1+3) - 2 s_1^2 \\[10pt]
        &= -2\chi_2^2 + 4\chi_2 s_1 - 2 s_1^2 \\[10pt]
        &= -2\left(\chi_2^2 - 2\chi_2 s_1 + s_1^2\right) \\[10pt]
        &= -2(s_1-\chi_2)^2 \; , \\[10pt]
        &\underset{\Delta\chi_2 \rightarrow 0^{+}}{\sim}\mathcal{O}(\Delta\chi_2^2) \\[16pt]
(\mathrm{D}.12) : \quad N_{02}^{(f)} &:= \left[
        4\chi_2^3 (3y^2-1) - 2\chi_2^2 y s_1 (3y^2+8) + \chi_2 s_1^2 (9y^2+11) - 6 y s_1^3 \right] \\[10pt]
    \Rightarrow \quad
    N_{02}^{(f)}\big|_{y=1} &= 4\chi_2^3 (3-1) - 2\chi_2^2 s_1 (3+8) + \chi_2 s_1^2 (9+11) - 6 s_1^3 \\[10pt]
        &= 8\chi_2^3 - 22\chi_2^2 s_1 + 20\chi_2 s_1^2 - 6 s_1^3 \\[10pt]
        &= -2\left(3 s_1^3 - 10 s_1^2\chi_2 + 11 s_1\chi_2^2 - 4\chi_2^3\right) \\[10pt]
        &= -2\left[3 s_1^3 + (-6-4) s_1^2\chi_2 + (3+8) s_1\chi_2^2 - 4\chi_2^3\right] \\[10pt]
        &= -2\left[\left(3 s_1^3 - 6 s_1^2\chi_2 + 3 s_1\chi_2^2\right)
            - \left(4 s_1^2\chi_2 - 8 s_1\chi_2^2 + 4\chi_2^3\right)\right] \\[10pt]
        &= -2\left[3 s_1\left(s_1^2 - 2 s_1\chi_2 + \chi_2^2\right)
            - 4\chi_2\left(s_1^2 - 2 s_1\chi_2 + \chi_2^2\right)\right] \\[10pt]
        &= -2\left(3 s_1 - 4\chi_2\right)\left(s_1^2 - 2 s_1\chi_2 + \chi_2^2\right) \\[10pt]
        &= -2(s_1-\chi_2)^2\left(3 s_1 - 4\chi_2\right) \; , \\[10pt]
        &\underset{\Delta\chi_2 \rightarrow 0^{+}}{\sim}\mathcal{O}(\Delta\chi_2^2)
\end{align*}
```

**The ``k=1`` coefficients.** Differentiating once in ``y``, and evaluating at ``y=1``:

```math
\begin{align*}
\frac{\partial N_{02}^{(b)}}{\partial y}
    &= \frac{\partial}{\partial y}\left[-2\chi_2^2 y + \chi_2 s_1 (y^2+3) - 2 y s_1^2\right] \\[10pt]
    &= -2\chi_2^2 + 2\chi_2 s_1 y - 2 s_1^2 \; , \\[10pt]
\Rightarrow \quad
\frac{\partial N_{02}^{(b)}}{\partial y}\bigg|_{y=1}
    &= -2\chi_2^2 + 2\chi_2 s_1 - 2 s_1^2 \\[10pt]
    &= -2\left(s_1^2 - s_1\chi_2 + \chi_2^2\right) \; , \\[10pt]
    &\underset{\Delta\chi_2 \rightarrow 0^{+}}{\sim}\mathcal{O}(1) \\[16pt]
\frac{\partial N_{02}^{(f)}}{\partial y}
    &= \frac{\partial}{\partial y}\left[
        4\chi_2^3 (3y^2-1) - 2\chi_2^2 y s_1 (3y^2+8) + \chi_2 s_1^2 (9y^2+11) - 6 y s_1^3\right] \\[10pt]
    &= 24 y\chi_2^3 - 2\chi_2^2 s_1 (9y^2+8) + 18 y\chi_2 s_1^2 - 6 s_1^3 \; , \\[10pt]
\Rightarrow \quad
\frac{\partial N_{02}^{(f)}}{\partial y}\bigg|_{y=1}
    &= 24\chi_2^3 - 2\chi_2^2 s_1 (9+8) + 18\chi_2 s_1^2 - 6 s_1^3 \\[10pt]
    &= 24\chi_2^3 - 34\chi_2^2 s_1 + 18\chi_2 s_1^2 - 6 s_1^3 \\[10pt]
    &= -2\left(3 s_1^3 - 9 s_1^2\chi_2 + 17 s_1\chi_2^2 - 12\chi_2^3\right) \; , \\[10pt]
    &\underset{\Delta\chi_2 \rightarrow 0^{+}}{\sim}\mathcal{O}(1)
\end{align*}
```

The ``k=2`` and ``k=3`` coefficients are also ``\mathcal{O}(1)``, so that, counting the orders:

```math
\begin{align*}
    &N_{02}\big|_{y=1} \cdot 1 &&=
        \mathcal{O}(\Delta\chi_2^2) \cdot 1
        &&= \mathcal{O}(\Delta\chi_2^2)
        &&\Longrightarrow \; \mathrm{keep} \; , \\[6pt]
    &\frac{\partial N_{02}}{\partial y}\bigg|_{y=1}(y-1) &&=
        \mathcal{O}(1) \cdot \mathcal{O}(\Delta\chi_2^2)
        &&= \mathcal{O}(\Delta\chi_2^2)
        &&\Longrightarrow \; \mathrm{keep} \; , \\[6pt]
    \frac{1}{2}&\frac{\partial^2 N_{02}}{\partial y^2}\bigg|_{y=1} (y-1)^2 &&=
        \mathcal{O}(1) \cdot \mathcal{O}(\Delta\chi_2^4)
        &&= \mathcal{O}(\Delta\chi_2^4)
        &&\Longrightarrow \; \mathrm{drop} \; , \\[6pt]
    \frac{1}{6}&\frac{\partial^3 N_{02}}{\partial y^3}\bigg|_{y=1}(y-1)^3 &&=
        \mathcal{O}(1) \cdot \mathcal{O}(\Delta\chi_2^6)
        &&= \mathcal{O}(\Delta\chi_2^6)
        &&\Longrightarrow \; \mathrm{drop} \; .
\end{align*}
```

**Putting the two surviving orders together.** With ``(s_1 - \chi_2) = -p\Delta\chi_2`` from
(D.4), ``(y-1)`` from (2.3a), and the regular factors evaluated at ``\chi_2 = s_1``
(``3-4 = -1`` and ``3-9+17-12 = -1``):

```math
\begin{align*}
N_{02}^{(b)} &= \sum_{k} \frac{1}{k!}\frac{\partial^k N_{02}^{(b)}}{\partial y^k}\bigg|_{y=1} (y-1)^k \\[15pt]
    &\underset{\Delta\chi_2\rightarrow 0^{+}}{\sim}
        N_{02}^{(b)}\big|_{y=1} \cdot 1 + \frac{\partial N_{02}^{(b)}}{\partial y}\bigg|_{y=1}(y-1) \\[15pt]
    &= -2(s_1-\chi_2)^2 \cdot 1 - 2\left(s_1^2 - s_1\chi_2 + \chi_2^2\right)\cdot(y-1) \\[15pt]
    &\underset{\Delta\chi_2\rightarrow 0^{+}}{\sim}
        \underbrace{-2(-p\Delta\chi_2)^2}_{k=0}
        + \underbrace{\left(-2 s_1^2\right)\left[- \frac{(1-p^2)}{2 s_1^2}\Delta\chi_2^2\right]}_{k=1} \\[10pt]
    &= -2 p^2\Delta\chi_2^2 + (1-p^2)\Delta\chi_2^2 \\[10pt]
    &= \left[-2p^2 + (1-p^2)\right]\Delta\chi_2^2 \\[10pt]
    &= \left(1 - 3p^2\right)\Delta\chi_2^2 \; , \\[20pt]
N_{02}^{(f)} &= \sum_{k} \frac{1}{k!}\frac{\partial^k N_{02}^{(f)}}{\partial y^k}\bigg|_{y=1} (y-1)^k \\[15pt]
    &\underset{\Delta\chi_2\rightarrow 0^{+}}{\sim}
        N_{02}^{(f)}\big|_{y=1} \cdot 1 + \frac{\partial N_{02}^{(f)}}{\partial y}\bigg|_{y=1}(y-1) \\[15pt]
    &= -2(s_1-\chi_2)^2\left(3 s_1 - 4\chi_2\right) \cdot 1
        - 2\left(3 s_1^3 - 9 s_1^2\chi_2 + 17 s_1\chi_2^2 - 12\chi_2^3\right)\cdot(y-1) \\[15pt]
    &\underset{\Delta\chi_2\rightarrow 0^{+}}{\sim}
        \underbrace{-2(-p\Delta\chi_2)^2 \left(-s_1\right)}_{k=0}
        + \underbrace{\left(-2\right)\left(- s_1^3\right)
            \left[- \frac{(1-p^2)}{2 s_1^2}\Delta\chi_2^2\right]}_{k=1} \\[10pt]
    &= 2 s_1 p^2\Delta\chi_2^2 - s_1(1-p^2)\Delta\chi_2^2 \\[10pt]
    &= s_1\left[2p^2 - (1-p^2)\right]\Delta\chi_2^2 \\[10pt]
    &= s_1\left(3p^2 - 1\right)\Delta\chi_2^2 \; .
\end{align*}
```

Combining them, the two ``\mathcal{O}(\Delta\chi_2^2)`` coefficients add up to

```math
\begin{align*}
    7 s_1 b_1 \, N_{02}^{(b)} + f_1 \, N_{02}^{(f)}
        &\underset{\Delta\chi_2\rightarrow 0^{+}}{\sim}
            7 s_1 b_1 \left(1-3p^2\right)\Delta\chi_2^2
            + f_1 s_1 \left(3p^2-1\right)\Delta\chi_2^2 \\[10pt]
        &= s_1\left(3p^2-1\right)\left[-7 b_1 + f_1\right]\Delta\chi_2^2 \\[10pt]
        &= s_1\left(3p^2-1\right)\left(f_1 - 7 b_1\right)\Delta\chi_2^2 \; ,
\end{align*}
```

so that the explicit ``\Delta\chi_2^{-2}`` is exactly cancelled and ``J_{02}`` is **finite**,
but it multiplies a vanishing ``I_2^0``:

```math
\begin{align*}
    J_{02}^{\delta \kappa} I_2^0
    &= \frac{1}{14\,\Delta\chi_2^2}\left[7 s_1 b_1 N_{02}^{(b)} + f_1 N_{02}^{(f)}\right] I_2^0 \\[12pt]
    &\underset{\Delta\chi_2\rightarrow 0^{+}}{\sim}
        \frac{1}{14\,\cancel{\Delta\chi_2^2}}
        \left[s_1\left(3p^2-1\right)\left(f_1 - 7 b_1\right)\cancel{\Delta\chi_2^2}\right]
        \cdot \frac{\sigma_{-2}}{15}\Delta\chi_2^2 \\[12pt]
    &= \frac{s_1\left(3p^2-1\right)\left(f_1 - 7 b_1\right)\sigma_{-2}}{210}\,\Delta\chi_2^2 \\[12pt]
    &= \mathcal{O}(\Delta\chi_2^2) \; \xrightarrow[\Delta\chi_2\rightarrow 0^{+}]{} \; 0 \; .
    \quad \quad (5.3)
\end{align*}
```

### Term 3: ``J_{04} I_4^0``

This is the delicate one: ``J_{04} \propto \Delta\chi_2^{-4}`` against
``I_4^0 \propto \Delta\chi_2^{4}``, so ``N_{04}``, Eq.(D.13), must be expanded to **fourth**
order, and checking that it vanishes at the point is nowhere near enough.

**The ``k=0`` coefficient.** Setting ``y=1``:

```math
\begin{align*}
(\mathrm{D}.13) : \quad N_{04} &:= \left[
        \chi_2^5 (6y^2-2) + 6\chi_2^4 y s_1 (y^2-3) - \chi_2^3 s_1^2 (y^4+12y^2-21)
        + 2\chi_2^2 y s_1^3 (y^2+3) - 12 \chi_2 s_1^4 + 4 y s_1^5 \right] \\[10pt]
    \Rightarrow \quad
    N_{04}\big|_{y=1} &= \chi_2^5 (6-2) + 6\chi_2^4 s_1 (1-3) - \chi_2^3 s_1^2 (1+12-21)
        + 2\chi_2^2 s_1^3 (1+3) - 12\chi_2 s_1^4 + 4 s_1^5 \\[10pt]
        &= 4\chi_2^5 - 12\chi_2^4 s_1 + 8\chi_2^3 s_1^2 + 8\chi_2^2 s_1^3 - 12\chi_2 s_1^4 + 4 s_1^5 \\[10pt]
        &= 4\left(s_1^5 - 3 s_1^4\chi_2 + 2 s_1^3\chi_2^2 + 2 s_1^2\chi_2^3 - 3 s_1\chi_2^4 + \chi_2^5\right) \\[10pt]
        &= 4\left[s_1^5 + (-4+1) s_1^4\chi_2 + (6-4) s_1^3\chi_2^2
            + (-4+6) s_1^2\chi_2^3 + (1-4) s_1\chi_2^4 + \chi_2^5\right] \\[10pt]
        &= 4\left[\left(s_1^5 - 4 s_1^4\chi_2 + 6 s_1^3\chi_2^2 - 4 s_1^2\chi_2^3 + s_1\chi_2^4\right)
            \right. \\
        &\phantom{= 4[} \left.
            + \left(s_1^4\chi_2 - 4 s_1^3\chi_2^2 + 6 s_1^2\chi_2^3 - 4 s_1\chi_2^4 + \chi_2^5\right)\right] \\[10pt]
        &= 4\left[s_1\left(s_1^4 - 4 s_1^3\chi_2 + 6 s_1^2\chi_2^2 - 4 s_1\chi_2^3 + \chi_2^4\right)
            \right. \\
        &\phantom{= 4[} \left.
            + \chi_2\left(s_1^4 - 4 s_1^3\chi_2 + 6 s_1^2\chi_2^2 - 4 s_1\chi_2^3 + \chi_2^4\right)\right] \\[10pt]
        &= 4\left(s_1 + \chi_2\right)\left(s_1^4 - 4 s_1^3\chi_2 + 6 s_1^2\chi_2^2 - 4 s_1\chi_2^3 + \chi_2^4\right) \\[10pt]
        &= 4(s_1-\chi_2)^4\left(s_1+\chi_2\right) \; , \\[10pt]
        &\underset{\Delta\chi_2 \rightarrow 0^{+}}{\sim}\mathcal{O}(\Delta\chi_2^4)
\end{align*}
```

**The ``k=1`` and ``k=2`` coefficients.** Differentiating,

```math
\begin{align*}
\frac{\partial N_{04}}{\partial y}
    &= \frac{\partial}{\partial y}\left[
        \chi_2^5 (6y^2-2) + 6\chi_2^4 y s_1 (y^2-3) - \chi_2^3 s_1^2 (y^4+12y^2-21)
        + 2\chi_2^2 y s_1^3 (y^2+3) - 12 \chi_2 s_1^4 + 4 y s_1^5\right] \\[10pt]
    &= 12 y\chi_2^5 + 6\chi_2^4 s_1 (3y^2-3) - \chi_2^3 s_1^2 (4y^3+24y)
        + 2\chi_2^2 s_1^3 (3y^2+3) + 4 s_1^5 \\[10pt]
\Rightarrow \quad
\frac{\partial N_{04}}{\partial y}\bigg|_{y=1}
    &= 12\chi_2^5 + 0 - 28\chi_2^3 s_1^2 + 12\chi_2^2 s_1^3 + 4 s_1^5 \\[10pt]
    &= 4\left(3\chi_2^5 - 7\chi_2^3 s_1^2 + 3\chi_2^2 s_1^3 + s_1^5\right) \; ,
\end{align*}
```

which must contain ``(s_1-\chi_2)^2``, since it vanishes twice there; dividing it out,

```math
    3\chi_2^5 - 7\chi_2^3 s_1^2 + 3\chi_2^2 s_1^3 + s_1^5
    = (s_1-\chi_2)^2\left(3\chi_2^3 + 6\chi_2^2 s_1 + 2\chi_2 s_1^2 + s_1^3\right) \; ,
```

(check, by expanding the right-hand side:
``3\chi_2^5 + 6\chi_2^4 s_1 + 2\chi_2^3 s_1^2 + \chi_2^2 s_1^3
- 6\chi_2^4 s_1 - 12\chi_2^3 s_1^2 - 4\chi_2^2 s_1^3 - 2\chi_2 s_1^4
+ 3\chi_2^3 s_1^2 + 6\chi_2^2 s_1^3 + 2\chi_2 s_1^4 + s_1^5
= 3\chi_2^5 - 7\chi_2^3 s_1^2 + 3\chi_2^2 s_1^3 + s_1^5`` ✓), so that

```math
\begin{align*}
    \frac{\partial N_{04}}{\partial y}\bigg|_{y=1}
    &= 4(s_1-\chi_2)^2\left(3\chi_2^3 + 6\chi_2^2 s_1 + 2\chi_2 s_1^2 + s_1^3\right) \; , \\[10pt]
    &\underset{\Delta\chi_2 \rightarrow 0^{+}}{\sim}\mathcal{O}(\Delta\chi_2^2) \\[16pt]
    \frac{\partial^2 N_{04}}{\partial y^2}
    &= \frac{\partial}{\partial y}\left[
        12 y\chi_2^5 + 6\chi_2^4 s_1 (3y^2-3) - \chi_2^3 s_1^2 (4y^3+24y)
        + 2\chi_2^2 s_1^3 (3y^2+3) + 4 s_1^5\right] \\[10pt]
    &= 12\chi_2^5 + 36 y\chi_2^4 s_1 - \chi_2^3 s_1^2 (12y^2+24) + 12 y\chi_2^2 s_1^3 \\[10pt]
    \Rightarrow \quad
    \frac{1}{2}\frac{\partial^2 N_{04}}{\partial y^2}\bigg|_{y=1}
    &= \frac{1}{2}\left(12\chi_2^5 + 36\chi_2^4 s_1 - 36\chi_2^3 s_1^2 + 12\chi_2^2 s_1^3\right) \\[10pt]
    &= 6\chi_2^5 + 18\chi_2^4 s_1 - 18\chi_2^3 s_1^2 + 6\chi_2^2 s_1^3 \\[10pt]
    &= 6\chi_2^2\left(\chi_2^3 + 3\chi_2^2 s_1 - 3\chi_2 s_1^2 + s_1^3\right) \; , \\[10pt]
    &\underset{\Delta\chi_2 \rightarrow 0^{+}}{\sim}\mathcal{O}(1)
\end{align*}
```

The ``k=3`` and ``k=4`` coefficients are ``\mathcal{O}(1)`` as well. Counting the orders:

```math
\begin{align*}
    &N_{04}\big|_{y=1} \cdot 1 &&=
        \mathcal{O}(\Delta\chi_2^4) \cdot 1
        &&= \mathcal{O}(\Delta\chi_2^4)
        &&\Longrightarrow \; \mathrm{keep} \; , \\[6pt]
    &\frac{\partial N_{04}}{\partial y}\bigg|_{y=1}(y-1) &&=
        \mathcal{O}(\Delta\chi_2^2) \cdot \mathcal{O}(\Delta\chi_2^2)
        &&= \mathcal{O}(\Delta\chi_2^4)
        &&\Longrightarrow \; \mathrm{keep} \; , \\[6pt]
    \frac{1}{2}&\frac{\partial^2 N_{04}}{\partial y^2}\bigg|_{y=1} (y-1)^2 &&=
        \mathcal{O}(1) \cdot \mathcal{O}(\Delta\chi_2^4)
        &&= \mathcal{O}(\Delta\chi_2^4)
        &&\Longrightarrow \; \mathrm{keep} \; , \\[6pt]
    \frac{1}{6}&\frac{\partial^3 N_{04}}{\partial y^3}\bigg|_{y=1}(y-1)^3 &&=
        \mathcal{O}(1) \cdot \mathcal{O}(\Delta\chi_2^6)
        &&= \mathcal{O}(\Delta\chi_2^6)
        &&\Longrightarrow \; \mathrm{drop} \; , \\[6pt]
    \frac{1}{24}&\frac{\partial^4 N_{04}}{\partial y^4}\bigg|_{y=1}(y-1)^4 &&=
        \mathcal{O}(1) \cdot \mathcal{O}(\Delta\chi_2^8)
        &&= \mathcal{O}(\Delta\chi_2^8)
        &&\Longrightarrow \; \mathrm{drop} \; .
\end{align*}
```

**Putting the three surviving orders together.** With the regular factors at
``\chi_2 = s_1`` (``1+1 = 2`` , ``3+6+2+1 = 12`` and ``1+3-3+1 = 2``):

```math
\begin{align*}
N_{04} &= \sum_{k} \frac{1}{k!}\frac{\partial^k N_{04}}{\partial y^k}\bigg|_{y=1} (y-1)^k \\[15pt]
    &\underset{\Delta\chi_2\rightarrow 0^{+}}{\sim}
        N_{04}\big|_{y=1} \cdot 1 + \frac{\partial N_{04}}{\partial y}\bigg|_{y=1}(y-1)
        + \frac{1}{2}\frac{\partial^2 N_{04}}{\partial y^2}\bigg|_{y=1}(y-1)^2 \\[15pt]
    &= 4(s_1-\chi_2)^4\left(s_1+\chi_2\right) \cdot 1
        + 4(s_1-\chi_2)^2\left(3\chi_2^3 + 6\chi_2^2 s_1 + 2\chi_2 s_1^2 + s_1^3\right)(y-1) \\
    &\phantom{=} + 6\chi_2^2\left(\chi_2^3 + 3\chi_2^2 s_1 - 3\chi_2 s_1^2 + s_1^3\right)(y-1)^2 \\[15pt]
    &\underset{\Delta\chi_2\rightarrow 0^{+}}{\sim}
        \underbrace{4(-p\Delta\chi_2)^4\left(2 s_1\right)}_{k=0}
        + \underbrace{4(-p\Delta\chi_2)^2\left(12 s_1^3\right)
            \left[- \frac{(1-p^2)}{2 s_1^2}\Delta\chi_2^2\right]}_{k=1} \\
    &\phantom{\underset{\Delta\chi_2\rightarrow 0^{+}}{\sim}}
        + \underbrace{6 s_1^2\left(2 s_1^3\right)
            \left[- \frac{(1-p^2)}{2 s_1^2}\Delta\chi_2^2\right]^{2}}_{k=2} \\[10pt]
    &= 8 s_1 p^4 \Delta\chi_2^4 - 24 s_1 p^2(1-p^2)\Delta\chi_2^4
        + 12\cancel{s_1^5}\frac{(1-p^2)^2}{4\cancel{s_1^4}}\Delta\chi_2^4 \\[10pt]
    &= s_1\left[8p^4 - 24p^2(1-p^2) + 3(1-p^2)^2\right]\Delta\chi_2^4 \\[10pt]
    &= s_1\left[8p^4 - 24p^2 + 24p^4 + 3 - 6p^2 + 3p^4\right]\Delta\chi_2^4 \\[10pt]
    &= s_1\left(35p^4 - 30p^2 + 3\right)\Delta\chi_2^4
    \; = \; 8 \, s_1 \, \mathcal{L}_4(p) \, \Delta\chi_2^4 \; ,
\end{align*}
```

the very same ``8\mathcal{L}_4(p)`` that ``B_{22}`` produced in Family 1, Eq.(3.9). The
explicit ``\Delta\chi_2^{-4}`` is therefore exactly cancelled and ``J_{04}`` is **finite**,
while ``I_4^0`` vanishes as ``\Delta\chi_2^4``:

```math
\begin{align*}
    J_{04}^{\delta \kappa} I_4^0
    &= \frac{3 f_1}{70\,\Delta\chi_2^4} N_{04} \, I_4^0 \\[12pt]
    &\underset{\Delta\chi_2\rightarrow 0^{+}}{\sim}
        \frac{3 f_1}{70\,\cancel{\Delta\chi_2^4}}
        \left[s_1\left(35p^4 - 30p^2 + 3\right)\cancel{\Delta\chi_2^4}\right]
        \cdot \frac{\sigma_{-4}}{945}\Delta\chi_2^4 \\[12pt]
    &= \frac{f_1 s_1 \left(35p^4-30p^2+3\right)\sigma_{-4}}{22050}\,\Delta\chi_2^4 \\[12pt]
    &= \mathcal{O}(\Delta\chi_2^4) \; \xrightarrow[\Delta\chi_2\rightarrow 0^{+}]{} \; 0 \; .
    \quad \quad (5.4)
\end{align*}
```

### The sum

```math
\begin{align*}
(5.2) : \quad J_{00}^{\delta \kappa} I_0^0 &\underset{\Delta\chi_2\rightarrow 0^{+}}{\sim} \mathcal{O}(1)
    &&\rightarrow -\frac{s_1 (f_1+5 b_1)\,\sigma_0}{5} \\[10pt]
(5.3) : \quad J_{02}^{\delta \kappa} I_2^0 &\underset{\Delta\chi_2\rightarrow 0^{+}}{\sim} \mathcal{O}(\Delta\chi_2^2)
    &&\rightarrow 0 \\[10pt]
(5.4) : \quad J_{04}^{\delta \kappa} I_4^0 &\underset{\Delta\chi_2\rightarrow 0^{+}}{\sim} \mathcal{O}(\Delta\chi_2^4)
    &&\rightarrow 0 \\[10pt]
\end{align*}
```

```math
    \boxed{\;
    \lim_{\Delta\chi_2 \rightarrow 0}
    \left(J_{00}^{\delta \kappa}I_0^0 + J_{02}^{\delta \kappa}I_2^0
    + J_{04}^{\delta \kappa}I_4^0\right)
    = -\frac{s_1 \, (f_1 + 5 b_1) \, \sigma_0}{5} \; . }
    \quad \quad (5.5)
```

The two terms that do carry a ``p`` are the ones whose limit is ``0``, so nothing has to
cancel.
