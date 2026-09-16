# Family 3: Newtonian ``\times`` Lensing

## Recap: everything this page needs

This page is self-contained. All the results quoted here are derived in
The ``\Delta\chi \rightarrow 0`` limits, and the ``I_\ell^n``
asymptotics in The ``I_\ell^n`` integrals.

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

Concerned functions: `integrand_ξ_GNC_Newtonian_Lensing`,
`integrand_ξ_GNCxLD_Newtonian_Lensing`. The sum to be taken to the limit is

```math
    J_{00}I_0^0 + J_{02}I_2^0 + J_{04}I_4^0 \; ,
```

with (``f`` the growth rate and ``b`` the bias, both evaluated at ``s_1``)

```math
\begin{align*}
    J_{00} &= \frac{1}{5}\left[ f \chi_2 (3y^2-1) - 3 y s_1 f - 5 y s_1 b \right] \; , \\[6pt]
    J_{02} &= \frac{1}{14\,\Delta\chi_2^2}\Big[
        7 s_1 b \underbrace{\left(-2\chi_2^2 y + \chi_2 s_1 (y^2+3) - 2 y s_1^2\right)}_{=: \; N_{02}^{(b)}} \\
        &\phantom{= \frac{1}{14\,\Delta\chi_2^2}\Big[} + f\underbrace{\left(
            4\chi_2^3 (3y^2-1) - 2\chi_2^2 y s_1 (3y^2+8) + \chi_2 s_1^2 (9y^2+11) - 6 y s_1^3
        \right)}_{=: \; N_{02}^{(f)}} \Big] \; , \\[6pt]
    J_{04} &= \frac{3 \, f}{70\,\Delta\chi_2^4} \underbrace{\Big[
        \chi_2^5 (6y^2-2) + 6\chi_2^4 y s_1 (y^2-3) - \chi_2^3 s_1^2 (y^4+12y^2-21)
        + 2\chi_2^2 y s_1^3 (y^2+3) - 12 \chi_2 s_1^4 + 4 y s_1^5 \Big]}_{=: \; N_{04}} \; .
\end{align*}
```

### Term 1: ``J_{00} I_0^0``

This family is the exception: ``J_{00}`` is regular **and does not vanish** at the singular
point, so it is the only one whose surviving contribution comes from an ``I_0^0`` rather than
from an ``I_0^2`` or an ``I_1^3``. Precisely because nothing cancels, here we *may* set
``y = 1`` and ``\chi_2 = s_1`` straight away — the trap of (2.4) only bites when a bracket
vanishes:

```math
\begin{align*}
    J_{00} &= \frac{1}{5}\left[ f \chi_2 (3y^2-1) - 3 y s_1 f - 5 y s_1 b \right] \\[10pt]
    y = 1 \; , \; \chi_2 = s_1 \; \rightarrow \quad
        &\xrightarrow[\Delta\chi_2\rightarrow 0^{+}]{}
            \frac{1}{5}\left[ f s_1 (3-1) - 3 s_1 f - 5 s_1 b \right] \\[10pt]
        &= \frac{1}{5}\left[ 2 f s_1 - 3 f s_1 - 5 b s_1 \right] \\[10pt]
        &= \frac{1}{5}\left[ - f s_1 - 5 b s_1 \right] \\[10pt]
        &= -\frac{s_1 \, (f + 5b)}{5} \; ,
    \quad \quad (5.1)
\end{align*}
```

and therefore, with ``I_0^0 \rightarrow \sigma_0`` from (2.1a),

```math
    J_{00} I_0^0 \; \xrightarrow[\Delta\chi_2\rightarrow 0^{+}]{} \;
        -\frac{s_1 \, (f + 5b) \, \sigma_0}{5} \; .
    \quad \quad (5.2)
```

Being the value of a regular function at a point, this carries no ``p``: the
direction-independence is automatic for this term.

### Term 2: ``J_{02} I_2^0``

``J_{02}`` carries ``\Delta\chi_2^{-2}`` and ``I_2^0`` carries ``\Delta\chi_2^{+2}``, so the
two exactly compensate and the product is **finite**: what decides is the order at which the
two numerators vanish. They both vanish at the singular point:

```math
\begin{align*}
    N_{02}^{(b)}\big|_{y=1,\,\chi_2=s_1}
        &= \left(-2 + (1+3) - 2\right)s_1^2 = (-2 + 4 - 2)\,s_1^2 = 0 \; , \\[6pt]
    N_{02}^{(f)}\big|_{y=1,\,\chi_2=s_1}
        &= \left(4(3-1) - 2(3+8) + (9+11) - 6\right)s_1^3 = (8 - 22 + 20 - 6)\,s_1^3 = 0 \; ,
\end{align*}
```

and expanding both along ``\chi_2 = s_1 + p\,\Delta\chi_2`` gives their common
``\mathcal{O}(\Delta\chi_2^2)`` coefficient:

```math
\begin{align*}
    N_{02}^{(b)} &\underset{\Delta\chi_2\rightarrow 0^{+}}{\sim} \left(1 - 3p^2\right)\Delta\chi_2^2 \; ,
    \qquad
    N_{02}^{(f)} \underset{\Delta\chi_2\rightarrow 0^{+}}{\sim} s_1\left(3p^2 - 1\right)\Delta\chi_2^2 \; , \\[10pt]
    \Rightarrow \quad
    7 s_1 b \, N_{02}^{(b)} + f \, N_{02}^{(f)}
        &\underset{\Delta\chi_2\rightarrow 0^{+}}{\sim}
            \left[7 s_1 b \left(1-3p^2\right) + f s_1 \left(3p^2-1\right)\right]\Delta\chi_2^2 \\[10pt]
        &= s_1\left(3p^2-1\right)\left(f - 7b\right)\Delta\chi_2^2 \; .
\end{align*}
```

The explicit ``\Delta\chi_2^{-2}`` is thus exactly cancelled and ``J_{02}`` is finite,

```math
    J_{02} \underset{\Delta\chi_2\rightarrow 0^{+}}{\sim}
        \frac{s_1\left(3p^2-1\right)\left(f - 7b\right)}{14}
        = \frac{s_1 \, \mathcal{L}_2(p) \left(f - 7b\right)}{7} \; ,
```

but it multiplies a **vanishing** ``I_2^0``:

```math
    J_{02} I_2^0 \underset{\Delta\chi_2\rightarrow 0^{+}}{\sim}
        \underbrace{\frac{s_1\left(3p^2-1\right)\left(f - 7b\right)}{14}}_{\mathcal{O}(1)} \cdot
        \underbrace{\frac{\sigma_{-2}}{15}\Delta\chi_2^2}_{\mathcal{O}(\Delta\chi_2^2)}
    \; \xrightarrow[\Delta\chi_2\rightarrow 0^{+}]{} \; 0 \; .
    \quad \quad (5.3)
```

### Term 3: ``J_{04} I_4^0``

This is the delicate one: ``J_{04} \propto \Delta\chi_2^{-4}`` against
``I_4^0 \propto \Delta\chi_2^{4}``, so ``N_{04}`` must be expanded to **fourth** order, and
checking that it vanishes at the point is nowhere near enough. It does vanish there,

```math
\begin{align*}
    N_{04}\big|_{y=1,\,\chi_2=s_1} &= \left[(6-2) + 6(1-3) - (1+12-21) + 2(1+3) - 12 + 4\right]s_1^5 \\[6pt]
        &= \left[4 - 12 + 8 + 8 - 12 + 4\right]s_1^5 = 0 \; ,
\end{align*}
```

and the expansion along ``\chi_2 = s_1 + p\,\Delta\chi_2`` gives

```math
    N_{04} \underset{\Delta\chi_2\rightarrow 0^{+}}{\sim}
        s_1 \left(35p^4 - 30p^2 + 3\right)\Delta\chi_2^4
    = 8 \, s_1 \, \mathcal{L}_4(p) \, \Delta\chi_2^4 \; ,
```

the very same ``8\mathcal{L}_4(p)`` that ``B_{22}`` produced in Family 1, Eq.(3.9). So here
too the explicit ``\Delta\chi_2^{-4}`` is exactly cancelled and ``J_{04}`` is **finite**,

```math
    J_{04} \underset{\Delta\chi_2\rightarrow 0^{+}}{\sim}
        \frac{3 \, f \, s_1 \left(35p^4 - 30p^2 + 3\right)}{70} \; ,
```

while ``I_4^0`` vanishes as ``\Delta\chi_2^4``:

```math
    J_{04} I_4^0 \underset{\Delta\chi_2\rightarrow 0^{+}}{\sim}
        \underbrace{\frac{3 \, f \, s_1 \left(35p^4 - 30p^2 + 3\right)}{70}}_{\mathcal{O}(1)} \cdot
        \underbrace{\frac{\sigma_{-4}}{945}\Delta\chi_2^4}_{\mathcal{O}(\Delta\chi_2^4)}
    \; \xrightarrow[\Delta\chi_2\rightarrow 0^{+}]{} \; 0 \; .
    \quad \quad (5.4)
```

### The sum

| term | limit |
|:--|:--|
| ``J_{00} I_0^0`` , Eq.(5.2) | ``-\dfrac{s_1 (f+5b)\,\sigma_0}{5}`` |
| ``J_{02} I_2^0`` , Eq.(5.3) | ``0`` |
| ``J_{04} I_4^0`` , Eq.(5.4) | ``0`` |

```math
    \boxed{\;
    \lim_{\Delta\chi_2 \rightarrow 0}
    \left(J_{00}I_0^0 + J_{02}I_2^0 + J_{04}I_4^0\right)
    = -\frac{s_1 \, (f + 5b) \, \sigma_0}{5} \; . }
    \quad \quad (5.5)
```

The two terms that do carry a ``p`` are the ones whose limit is ``0``, so again nothing has
to cancel.
