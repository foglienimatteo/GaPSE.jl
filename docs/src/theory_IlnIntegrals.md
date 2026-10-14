# The ``I_\ell^n`` integrals

```@contents
Pages = ["theory_IlnIntegrals.md"]
Depth = 3
```

Every Two-Point Correlation Function (TPCF) that GaPSE computes is, in the end, a sum of
terms of the form ``J(\chi, s, y) \, I_\ell^n(\Delta\chi)``. This page collects the
definition of these ``I_\ell^n``, proves their behaviour for small separations, and shows
what they look like.

The plots are produced by the notebook `theory/Iln_terms.ipynb`; see the end of this page.







## Definitions

The integrals are

```math
    I_\ell^n(s) := \int_{k_\mathrm{min}}^{k_\mathrm{max}} \frac{\mathrm{d}q}{2\pi^2} \, q^2 \, P(q) \, \frac{j_\ell(qs)}{(qs)^n} \;, \quad \quad (1)
```

with ``P(q)`` the matter Power Spectrum at ``z=0`` stored inside the `Cosmology``, and
``j_\ell`` the spherical Bessel function of order ``\ell``.

!!! note "The integration range is part of the definition"
    ``k_\mathrm{min}`` and ``k_\mathrm{max}`` are **not** a numerical detail to be sent to
    ``0`` and ``+\infty`` at the end: they are part of what ``I_\ell^n`` *means* here.
    `IPSTools` hands ``k_\mathrm{min}, k_\mathrm{max} = 10^{-5}, 10^{3}`` to `xicalc` for
    the ``I_\ell^n``, and every result on this page is a statement about *that* integral.
    The ``\sigma_i`` of (3) are cut elsewhere, for a reason given in
    [2. The ranges the code uses](@ref "2. The ranges the code uses"), and
    [Why the cut cannot be dropped](@ref "Why the cut cannot be dropped") shows what changes
    if one insists on ``\int_0^{+\infty}`` instead — it is not a harmless idealisation.

The eight combinations ``(\ell, n)`` that appear in the code are

```math
    (0,0) \, , \quad (2,0) \, , \quad (4,0) \, , \quad (0,2) \, , \quad
    (2,2) \, , \quad (3,1) \, , \quad (1,3) \, , \quad (1,1) \; ,
```

and they are stored in `IPSTools` as `I00`, `I20`, `I40`, `I02`, `I22`, `I31`, `I13` and
`I11` respectively, i.e. the name is `I` followed by ``\ell`` and then ``n``. Alongside
them there is the auxiliary integral

```math
\begin{align*}
    \tilde{I}_0^4(s) &:= \int_{k_\mathrm{min}}^{k_\mathrm{max}} \frac{\mathrm{d}q}{2\pi^2} \, q^2 \, P(q) \, \frac{j_0(qs) - 1 }{(qs)^4} \quad \quad (2) \\[10pt]
        &= \frac{1}{s^4} \int_{k_\mathrm{min}}^{k_\mathrm{max}} \frac{\mathrm{d}q}{2\pi^2} \, \frac{P(q)}{q^2} \, \left[ j_0(qs) - 1 \right] \; .
\end{align*}
```

stored in `IPSTools` as `I04_tilde` and computed by [`GaPSE.func_I04_tilde`](@ref).

We also need the moments of the Power Spectrum

```math
    \sigma_i := \int_{k_\mathrm{min}}^{k_\mathrm{max}} \frac{\mathrm{d}q}{2 \pi^2} \, q^{2-i} \, P(q) \; , \quad \quad (3)
```

over **the same** ``[k_\mathrm{min}, k_\mathrm{max}]`` as (1) and (2), of which `IPSTools`
stores ``\sigma_0``, ``\sigma_1``, ``\sigma_2``, ``\sigma_3`` and ``\sigma_4``. Note that
``\sigma_i`` with ``i<0`` also appear below; they are not stored, so the notebook
recomputes them when needed.

Every ``\sigma_i`` is a finite number, for every ``i``, because the range is finite. That
is worth stating explicitly, because some of these integrands do **not** decay fast enough
for the integral to exist over ``(0, +\infty)``: the moment is defined by its range, and
the range is the one in (1). See
[The input Power Spectrum](theory_InputPowerSpectrum.md) for which ones, and why that is a
property of ``P(q)`` rather than a defect of the code.



## TLDR; the easy-to-get small-``s`` behaviour


```math
\begin{align*}
j_\ell(x) \mathrm{ \; series \; expansion \; near \; 0}&: \quad \quad  
    j_\ell(x) = x^\ell \left( 
        \frac{\sqrt{\pi}}{2^{\ell+1}\, 
        \Gamma(\ell + 3/2)} + O(x^2)
    \right) \quad \quad \mathrm{(1.1)}\\[10pt]
\Rightarrow j_\ell(x) \; &\underset{x \rightarrow 0}{\sim}  x^{\ell}
\quad \quad \mathrm{(1.2)}
\end{align*}
```

So:

```math
\begin{align*}
    (1): \quad \quad I_\ell^n(s) 
    &:= \int_{k_\mathrm{min}}^{k_\mathrm{max}} \frac{\mathrm{d}q}{2\pi^2} \, q^2 \, P(q) \, \frac{j_\ell(qs)}{(qs)^n} \\[10pt]
    \mathrm{inserting\; (1.2) }  \; \rightarrow \;\; \quad
    &\underset{s \rightarrow 0}{\sim} \int_{k_\mathrm{min}}^{k_\mathrm{max}} \frac{\mathrm{d}q}{2\pi^2} \, q^2 \, P(q) \, \frac{(qs)^{\ell}}{(qs)^n}\\[10pt] 
    &=  s^{\ell-n}\int_{k_\mathrm{min}}^{k_\mathrm{max}} \frac{\mathrm{d}q}{2\pi^2} \, q^{2-(n-\ell)} \, P(q) \, \\[10pt]
    &= \sigma_{n-\ell} \; s^{\ell-n}\\[10pt]
    &\underset{s \rightarrow 0}{\sim} s^{\ell-n}
\end{align*}
```

The step that matters is the first one, and it is legitimate **only because the range is
finite**. Replacing ``j_\ell(qs)`` by its small-argument form requires ``qs \ll 1``; on
``[k_\mathrm{min}, k_\mathrm{max}]`` the largest argument is ``k_\mathrm{max} s``, so the
single condition

```math
    s \; \ll \; \frac{1}{k_\mathrm{max}}
    \quad \quad \mathrm{(1.3)}
```

makes ``qs \ll 1`` hold **uniformly** over the whole range, and the series may be inserted
and integrated term by term. Had the integral run to ``+\infty`` there would always be a
region ``q > 1/s`` in which the replacement is simply false, no matter how small ``s`` is,
and the manipulation above would be wrong. This is the content of
[Why the cut cannot be dropped](@ref "Why the cut cannot be dropped").


## The exact small-``s`` behaviour

```math
\boxed{
    I_\ell^n(s) \; \underset{s \rightarrow 0}{\sim}  \;
        \frac{\sigma_{n-\ell}}{(2\ell+1)!!} \, s^{\,\ell - n} 
}
\quad \quad (4a)
\qquad\qquad
\boxed{
    \tilde{I}_0^4(s) \; \underset{s \rightarrow 0}{\sim} \; - \frac{\sigma_2}{6 \, s^2}
}
\quad \quad (4b)
```

**Proof.** 

Acronym used for the sources:

- DLMF = Digital Library of Mathematical Functions
- NIST = National Institute of Standards and Technology



The Taylor series of the spherical Bessel function ``j_\ell(x)`` of order ``\ell`` is the following everywhere-convergent expansion (source: [U.S. NIST DLMF, Spherical Bessel Functions - Power Series, Section 10.53](https://dlmf.nist.gov/10.53) ):

```math
    j_\ell(x) = \sum_{k=0}^{+\infty} \frac{(-1)^k \, x^{\ell + 2k}}{2^k \, k! \, (2\ell + 2k + 1)!!} \; \quad \quad (2.1)\\[10pt]

 \quad n!! := \begin{cases} 
    2 \cdot 4 \cdot 6 \cdot ... \cdot n & n \mathrm{\; is \; even} \\
    1 \cdot 3 \cdot 5 \cdot ... \cdot n & n \mathrm{\; is \; odd}  \\
    1                                   & n=0,-1
\end{cases} 

```

We can set ``x = qs`` and divide both terms by ``(qs)^n``:

```math
    x:=qs \; \Rightarrow \;  \frac{(2.1)}{(qs)^n} : \quad \quad
    \frac{j_\ell(qs)}{(qs)^n} = \sum_{k=0}^{+\infty} \frac{(-1)^k \, (qs)^{\ell + 2k - n}}
        {2^k \, k! \, (2\ell + 2k + 1)!!} \; . \quad \quad (2.2)
```

We then insert this last Eq.(2.2) into the definition of ``I_\ell^n`` integrals Eq.(1), we bring the sum outside the integral (which is legitimate because the series converges uniformly on the compact integration range ``[k_\mathrm{min}, k_\mathrm{max}]``) and we insert the definition of ``\sigma_i`` Eq.(3):


```math
\begin{align*}
     (1): \quad \quad I_\ell^n(s) 
    &:= \int_0^{+\infty} \frac{\mathrm{d}q}{2\pi^2} \, q^2 \, P(q) \, \frac{j_\ell(qs)}{(qs)^n} \\[10pt]
    \mathrm{inserting\; } (2.2) \; \rightarrow \;\; \quad 
    &=  \int_0^{+\infty} \frac{\mathrm{d}q}{2\pi^2} \, q^2 \, P(q) \, 
        \sum_{k=0}^{+\infty} \frac{(-1)^k \, (qs)^{\ell + 2k - n}}
        {2^k \, k! \, (2\ell + 2k + 1)!!}  \\[10pt]
    &= \sum_{k=0}^{+\infty} \frac{(-1)^k \; s^{\ell + 2k - n}}
        {2^k \, k! \, (2\ell + 2k + 1)!!}
        \int \frac{\mathrm{d}q}{2\pi^2} \, q^{\,2 - (n - \ell - 2k)} \, P(q) \; \\[10pt]
    \mathrm{inserting\; } (3) \; \rightarrow \;\; \quad 
    &= \sum_{k=0}^{+\infty} \frac{(-1)^k \; \sigma_{n - \ell - 2k}}
        {2^k \, k! \, (2\ell + 2k + 1)!!} \; s^{\,\ell + 2k - n} \; . \quad \quad (2.3)
\end{align*}
```


Every term carries two more powers of ``s`` than the previous one, so for
``s \rightarrow 0`` the ``k=0`` term dominates, which is the claimed result:

```math
\begin{align*}
(2.3):\quad\quad I_\ell^n(s) &= \sum_{k=0}^{+\infty} \frac{(-1)^k \; \sigma_{n - \ell - 2k}}
        {2^k \, k! \, (2\ell + 2k + 1)!!} \; s^{\,\ell + 2k - n} \\[10pt]
    &= s^{\,\ell - n} \left[
        \frac{\sigma_{n-\ell}}{(2\ell+1)!!}  + 
        \sum_{k=1}^{+\infty} \frac{(-1)^k \; \sigma_{n - \ell - 2k}}
        {2^k \, k! \, (2\ell + 2k + 1)!!} \; s^{2k}\right]\\[10pt]
    &= \frac{\sigma_{n-\ell}}{(2\ell+1)!!} s^{\,\ell - n}
        \left[1 + \alpha_1 s^2 + \alpha_2 s^4 + ...\right]\\[10pt]
    &\underset{s \rightarrow 0}{\sim} \frac{\sigma_{n-\ell}}{(2\ell+1)!!} s^{\,\ell - n} \; . \quad\quad (2.4)
\end{align*}
```


For ``\tilde{I}_0^4`` the same argument applies to ``j_0(x) - 1``, whose series starts at
``k=1``:

```math
\begin{align*}
    (2.1) \mathrm{\;with\;}\ell=0 : \quad \quad 
    j_0(x) &= \sum_{k=0}^{+\infty} \frac{(-1)^k \, x^{2k}}{2^k \, k! \, (2k + 1)!!}\; \\[10pt]
        &= 1+\sum_{k=1}^{+\infty} \frac{(-1)^k \, x^{2k}}{2^k \, k! \, (2k + 1)!!} \\[10pt]
    2^k \, k! \, (2k+1)!! = (2k+1)! \; \rightarrow \quad \quad
        &= 1+\sum_{k=1}^{+\infty} \frac{(-1)^k \, x^{2k}}{(2k+1)!}
\end{align*}
```
```math
\Rightarrow \quad \quad j_0(x) -1 = \sum_{k=1}^{+\infty} \frac{(-1)^k \, x^{2k}}{(2k+1)!} \quad \quad (2.5)
```

Then

```math
\begin{align*}
    (2): \quad \quad \tilde{I}_0^4(s) 
        &= \frac{1}{s^4} \int_0^{+\infty} \frac{\mathrm{d}q}{2\pi^2} \, \frac{P(q)}{q^2} \, \left[ j_0(qs) - 1 \right] \; \\[10pt]
    \mathrm{inserting\; } (2.6) \; \rightarrow \;\; \quad
        &= \frac{1}{s^4} \int_0^{+\infty} \frac{\mathrm{d}q}{2\pi^2} \, \frac{P(q)}{q^2}
        \sum_{k=1}^{+\infty} \frac{(-1)^k \, (qs)^{2k}}{(2k+1)!}\\[10pt]
        &= \frac{1}{s^4} \sum_{k=1}^{+\infty} \frac{(-1)^k \, s^{2k}}{(2k+1)!}
        \int_0^{+\infty} \frac{\mathrm{d}q}{2\pi^2} \, q^{\,2k-2} \, P(q)\\[10pt]
        &= \sum_{k=1}^{+\infty} \frac{(-1)^k}{(2k+1)!}s^{2k-4}
        \int_0^{+\infty} \frac{\mathrm{d}q}{2\pi^2} \, q^{2-(4-2k)} \, P(q)\\[10pt]
    \mathrm{inserting\; } (3) \; \rightarrow \;\; \quad 
        &=\sum_{k=1}^{+\infty} \frac{(-1)^k \, \sigma_{4-2k}}{(2k+1)!} \, s^{\,2k-4} \; ,
        \quad\quad (2.7)
\end{align*}

```

And analogously, every term carries two more powers of ``s`` than the previous one, so for
``s \rightarrow 0`` the ``k=1`` term dominates, which is again the claimed result:

```math
\begin{align*}
(2.7):\quad\quad \tilde{I}_0^4(s) &= \sum_{k=1}^{+\infty} \frac{(-1)^k \, \sigma_{4-2k}}{(2k+1)!} \, s^{\,2k-4} \\[10pt]
    &= s^{-2} \left[
        \sum_{k=1}^{+\infty} \frac{(-1)^k \; \sigma_{4 - 2k}}
        {(2k + 1)!} \; s^{2k-2}\right]\\[10pt]
    &= s^{-2} \left[
        - \frac{\sigma_{2}}{3!} + \frac{\sigma_0}{5!}s^2 - \frac{\sigma_{-2}}{7!}s^4 +...
        \right]\\[10pt]
    &= - \frac{\sigma_{2}}{6} s^{-2} + \frac{\sigma_0}{120} - \frac{\sigma_{-2}}{5040} s^2 + ... \\[10pt]
    &\underset{s \rightarrow 0}{\sim} - \frac{\sigma_{2}}{6\,s^2} \; . \quad\quad (2.8)
\end{align*}
```







### Why the cut cannot be dropped

It is tempting to write (1) with ``\int_0^{+\infty}``, treat
``k_\mathrm{min}, k_\mathrm{max}`` as a numerical approximation to it, and let the
``\sigma_i`` inherit the same infinite range. That substitution is not legitimate, and it
is worth seeing exactly where it breaks: it produces a spurious divergence of
``I_0^0(s \rightarrow 0)`` that the integral of (1) does not have.

Split the infinite integral at the two cuts:

```math
    \int_0^{+\infty} = \underbrace{\int_0^{k_\mathrm{min}}}_{\mathrm{(A)}}
    + \underbrace{\int_{k_\mathrm{min}}^{k_\mathrm{max}}}_{\mathrm{(B)}}
    + \underbrace{\int_{k_\mathrm{max}}^{+\infty}}_{\mathrm{(C)}}
    \quad \quad \mathrm{(1.4)}
```

**(B) is the one we computed**, and the derivation above applies to it word for word:
``qs \leq k_\mathrm{max}s \ll 1`` uniformly, so

```math
    \mathrm{(B)} \; = \; \frac{\sigma_{n-\ell}}{(2\ell+1)!!} \, s^{\,\ell-n}
    \left[ 1 + \mathcal{O}\!\left( (k_\mathrm{max}s)^2 \right) \right] \; .
```

**(A) keeps the expansion but may lose the integral.** For ``q < k_\mathrm{min}`` the
argument ``qs`` is even smaller, so ``j_\ell`` may again be replaced by ``x^\ell``; what
can fail is the ``q`` integral itself. With ``P(q) \sim B\,q^{\,n_s}`` as
``q \rightarrow 0`` (see [The input Power Spectrum](theory_InputPowerSpectrum.md)),

```math
    \int_0 \mathrm{d}q \; q^{\,2-i+n_s}
    \quad \text{converges} \iff 2 - i + n_s > -1
    \iff i < 3 + n_s \simeq 3.96 \; .
```

So ``\sigma_0 \ldots \sigma_3`` would survive ``k_\mathrm{min} \rightarrow 0``, while
``\sigma_4`` would not. Nothing about the code changes: with ``k_\mathrm{min} > 0`` the
piece (A) is simply not part of the definition.

**(C) is where the expansion itself fails.** For ``q > 1/s`` the argument ``qs`` is large,
however small ``s`` is, so ``j_\ell(qs)`` cannot be replaced by ``(qs)^\ell``. The honest
way to evaluate it is to rescale. With ``P(q) \simeq A \, q^{-\alpha}`` at large ``q``
(``\alpha \simeq 2.64`` for the fitted tail, ``\alpha = 3`` for the asymptotic CDM one)
and ``x := qs``,

```math
\begin{align*}
    \mathrm{(C)} &= \frac{1}{2\pi^2} \int_{k_\mathrm{max}}^{+\infty} \mathrm{d}q \;
        q^2 \, A \, q^{-\alpha} \, \frac{j_\ell(qs)}{(qs)^n} \\[10pt]
    x := qs \; \rightarrow \; \quad
    &= \frac{A}{2\pi^2} \, s^{\,\alpha-3} \int_{k_\mathrm{max}s}^{+\infty} \mathrm{d}x \;
        x^{\,2-\alpha-n} \, j_\ell(x)
    \quad \quad \mathrm{(1.5)}
\end{align*}
```

The remaining integral tends to a finite number as ``s \rightarrow 0`` (its lower limit
goes to ``0`` and ``j_\ell(x) \sim \sin(x - \ell\pi/2)/x`` makes the upper end converge),
so **(C) scales as ``s^{\,\alpha-3}``** — a power that has nothing to do with
``s^{\,\ell-n}``. For the asymptotic CDM tail ``\alpha = 3`` it is a constant; for the
fitted ``\alpha \simeq 2.64`` it grows as ``s^{-0.36}``.

The conclusion is the important part:

!!! warning "``\sigma_{n-\ell} s^{\ell-n}`` is the limit of (1), not of the infinite integral"
    With the range of (1) the answer is ``\sigma_{n-\ell}s^{\ell-n}/(2\ell+1)!!``, every
    ``\sigma_i`` is finite, and the numerical check below confirms it to six digits.

    With ``\int_0^{+\infty}`` instead, the piece (C) adds a term ``\propto s^{\alpha-3}``
    that the naive manipulation never produces, because that manipulation assumed an
    expansion valid over the whole range. For ``I_0^0`` that term dominates and one
    concludes that the limit "diverges" — which is a statement about the infinite-range
    integral, not about anything GaPSE computes.

Concretely, for ``I_0^0`` over the range actually used, ``[10^{-5}, 10^{3}]``:

| ``s`` | ``I_0^0(s)`` | ``\sigma_0`` | ratio |
|:--|--:|--:|--:|
| ``10^{-1}`` | ``23.730`` | ``143.285`` | ``0.166`` |
| ``10^{-2}`` | ``68.766`` | ``143.285`` | ``0.480`` |
| ``10^{-3}`` | ``139.478`` | ``143.285`` | ``0.973`` |
| ``10^{-4}`` | ``143.246`` | ``143.285`` | ``0.99973`` |
| ``10^{-5}`` | ``143.285`` | ``143.285`` | ``0.999997`` |
| ``10^{-6}`` | ``143.285`` | ``143.285`` | ``1.000000`` |

The approach to the limit is exactly what (1.3) predicts. At ``s = 10^{-1}`` only the part
of the range with ``q \lesssim 10`` satisfies ``qs \ll 1``, and indeed
``I_0^0(0.1) \simeq \sigma_0(<10) = 18.6``; as ``s`` decreases more of the range comes
inside the condition, and at ``s \lesssim 1/k_\mathrm{max} = 10^{-3}`` all of it does and
the ratio saturates at ``1``. It is a monotone approach to a finite number, not a
divergence.

### The three regimes

The sign of ``\ell - n`` decides everything:

```math
I_\ell^n \xrightarrow[s \rightarrow 0]{} 
\begin{cases}
\propto s^{\ell-n} \rightarrow 0                    & \ell > n \\[10pt]
\frac{\sigma_0}{(2\ell+1)!!}                        & \ell = n  \quad \quad (5) \\[10pt]
\propto \frac{1}{s^{n-\ell}} \rightarrow +\infty    & \ell < n
\end{cases} 
```

So, written explicitly:


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
\end{align*}
```

```math
\tilde{I}_0^4   \sim -\frac{\sigma_2}{6}\, s^{-2}     \rightarrow +\infty\\[10pt]
```


The diverging ones are never a problem in practice: inside a TPCF they always come
multiplied by a ``J`` carrying the matching positive power of ``\Delta\chi``, so that the
product stays finite. How the cancellation works, term by term and for every TPCF, is the
subject of "The ``\Delta\chi \rightarrow 0`` limits" page.

## A warning before looking at the plots

Everything above is exact. What `IPSTools` gives you, however, is **not** the integral
everywhere, and a naive plot of it down to ``s = 10^{-4}`` looks nothing like (4a).
Three distinct things have to be kept in mind.

### 1. An `IntegralIPS` is a spline only between `left` and `right`

Each ``I_\ell^n`` is stored as a [`GaPSE.IntegralIPS`](@ref), which evaluates as

```math
I_\ell^n(s) =
\begin{cases}
a_\mathrm{L} + b_\mathrm{L} \, s^{\,\nu_\mathrm{L}} \; ,  & s < \mathrm{left} \\[6pt]
\mathrm{spline}(s) \; ,  & \mathrm{left} \leq s \leq \mathrm{right} \\[6pt]
a_\mathrm{R} + b_\mathrm{R} \, s^{\,\nu_\mathrm{R}} \; ,  & s > \mathrm{right}
\end{cases}
```

with ``\mathrm{left} = \mathrm{fit\_min} = 0.05 \, h_0^{-1}\mathrm{Mpc}`` for all the
``I_\ell^n`` (and ``0.1`` for ``\tilde{I}_0^4``). The coefficients
``a_\mathrm{L}, b_\mathrm{L}, \nu_\mathrm{L}`` are fitted on
``[\mathrm{fit\_min}, \mathrm{fit\_max}] = [0.05, 0.5]``, i.e. on a region that is still
very far from the asymptotic one, and the fit is seeded with a *negative* exponent.
The consequence is that, with this input Power Spectrum, below ``s = 0.05``
**every** ``I_\ell^n`` comes out with a negative fitted exponent and grows without bound,
whatever its true behaviour. The extrapolation is there to keep the object callable
outside the tabulated range, not to reproduce the limits: GaPSE never evaluates the
``I_\ell^n`` below ``\Delta\chi_\mathrm{min}``, where the analytic branch takes over. The
practical consequence for a reader is that the region ``s < 0.05`` of any plot of an
`IntegralIPS` carries no information about the limits derived above; in the figures below
it is shaded in grey.

### 2. The ranges the code uses

`IPSTools` builds its two families of objects over two different ``k`` ranges:

- the ``I_\ell^n`` come from a `xicalc` call with
  ``k_\mathrm{min}, k_\mathrm{max} = 10^{-5}, 10^{3}``, fixed in the constructor. A Hankel
  transform needs a wide, densely sampled ``k`` grid, and this one covers the whole ``s``
  interval over which the resulting splines are then evaluated;
- the ``\sigma_i`` come from a direct quadrature over the `k_min`/`k_max` keywords of
  `IPSTools`, whose defaults are ``10^{-6}`` and ``10``.

For ``\sigma_1``, ``\sigma_2`` and ``\sigma_3`` the distinction is immaterial: they
converge at both ends, so any wide enough range returns the same number. For ``\sigma_0``
and for the negative-index moments it is not, because those are dominated by their upper
extreme — with ``P(q) \propto q^{-2.64}`` at large ``q``, the integrand of ``\sigma_{-2}``
grows as ``q^{1.36}`` and that of ``\sigma_{-4}`` as ``q^{3.36}``. On
`data/WideA_ZA_pk.dat`:

| ``i``  | over ``[10^{-6}, 10]`` | over ``[10^{-5}, 10^{3}]`` |                 ratio |
| :----: | ---------------------: | -------------------------: | --------------------: |
| ``0``  |              ``18.58`` |                  ``143.3`` |               ``7.7`` |
| ``2``  |             ``101.06`` |                 ``101.13`` |             ``1.001`` |
| ``-2`` |              ``437.8`` |     ``2.35 \times 10^{7}`` | ``5.4 \times 10^{4}`` |
| ``-4`` | ``2.41 \times 10^{4}`` |    ``1.27 \times 10^{13}`` | ``5.3 \times 10^{8}`` |

Which range is the right one depends on what the moment is for, and inside GaPSE
``\sigma_0`` has exactly one job: it is the value the ``\Delta\chi \rightarrow 0`` branch
uses in place of ``I_0^0(\Delta\chi)`` below ``\Delta\chi_\mathrm{min}``. The best stand-in
for that number is the moment cut where the Bessel function stops contributing, i.e. at
``k_\mathrm{max} \simeq 1 / \Delta\chi_\mathrm{min}``, and with the default
``\Delta\chi_\mathrm{min} = 0.1`` that is ``k_\mathrm{max} = 10`` — the `IPSTools`
default. Measured:

| quantity | value |
| :-- | --: |
| ``\sigma_0`` over ``[10^{-6}, 10]`` (the default) | ``18.58`` |
| ``I_0^0(\Delta\chi_\mathrm{min}) = I_0^0(0.1)``   | ``23.73`` |
| ``\sigma_0`` over ``[10^{-5}, 10^{3}]`` (the `xicalc` range) | ``143.3`` |

so the default is ``22\%`` below the number it replaces, while the `xicalc` range would be
a factor ``6`` above it. Should ``\Delta\chi_\mathrm{min}`` be changed, `k_max` should
follow it.

The one place where the two ranges *must* be the same is a comparison between an
``I_\ell^n`` and the asymptote (4a) built out of ``\sigma_{n-\ell}``: these are two
expressions for the same integral, so they only agree if they integrate the same thing.
The tables and figures below therefore use the `xicalc` range for both. Only
``\sigma_2`` is insensitive enough for the distinction not to show, which is why
``I_0^2``, ``I_1^3`` and ``\tilde{I}_0^4`` — the three whose limits depend on
``\sigma_2`` alone — are the only ones that appear to obey (4a) even with mismatched
extremes.

### 3. The asymptotic regime starts below ``1/k_\mathrm{max}``

Truncating the series at ``k = 0`` requires ``(qs)^2 \ll 1`` for every ``q`` that
carries weight in the integral, i.e.

```math
s \; \ll \; \frac{1}{k_\mathrm{max}} = 10^{-3} \; h_0^{-1}\mathrm{Mpc} \; .
```

This is 50 times *smaller* than ``\mathrm{fit\_min} = 0.05``, so the window where the
limits hold and the window where the spline is valid **do not overlap**. The limits are
therefore not observable through `IPSTools`: to see them one has to compute the integrals
directly, which is what the `I_direct` function of `theory/Iln_terms.ipynb` does.

Doing so confirms (4a) to four digits. Calling ``R_\ell^n(s)`` the ratio between the
directly-computed integral and its asymptote:

| ``s``       | ``R_0^0`` | ``R_2^0`` | ``R_4^0`` | ``R_0^2`` | ``R_2^2`` | ``R_3^1`` | ``R_1^3`` | ``R_1^1`` |
| :---------- | --------: | --------: | --------: | --------: | --------: | --------: | --------: | --------: |
| ``10^{-5}`` |    1.0000 |    1.0000 |    1.0000 |    1.0000 |    1.0000 |    1.0000 |    1.0000 |    1.0000 |
| ``10^{-4}`` |    0.9997 |    0.9996 |    0.9997 |    1.0000 |    0.9999 |    0.9997 |    1.0000 |    0.9998 |
| ``10^{-3}`` |    0.9734 |    0.9621 |    0.9693 |    1.0000 |    0.9885 |    0.9704 |    1.0000 |    0.9839 |
| ``10^{-2}`` |    0.4799 |    0.0661 |    0.0416 |    1.0000 |    0.5997 |    0.1018 |    1.0000 |    0.5521 |
| ``0.05``    |    0.2337 |    0.0015 |    0.0000 |    0.9998 |    0.3033 |    0.0023 |    0.9999 |    0.2760 |

![Convergence to the limits](assets/Iln_terms/ratios.png)

The three ``\sigma_2``-only columns sit at 1 over the whole range, for the reason given
in the previous point: their integrand ``q^2 P(q) \, j_\ell(qs)/(qs)^n`` is weighted so
that only ``q \ll 1/s`` contributes, so the expansion is never actually stressed. All
the others peel off from 1 exactly where ``s`` approaches ``1/k_\mathrm{max}``.

## The plots

All the curves show ``|I_\ell^n(s)|`` rather than ``I_\ell^n(s)``: these integrals
oscillate and change sign at large ``s``, where a logarithmic vertical axis would not be
defined. At small ``s``, the regime this page is about, they keep a constant sign, so the
absolute value is immaterial there and the asymptote (dashed, black) can be read off
directly.

Each of the individual figures carries four things: the `IntegralIPS` stored in
`IPSTools` (solid), the direct quadrature for ``s \leq 10^{-2}`` (dotted), the analytic
asymptote (dashed, black), and two grey bands marking where the `IntegralIPS` is a
power-law extrapolation rather than the integral. The combined figure below is instead
restricted to ``[\mathrm{left}, \mathrm{right}]``, where all of them are genuine splines.

```math
\boxed{
    I_\ell^n(s) \; \underset{s \rightarrow 0}{\sim}  \;
        \frac{\sigma_{n-\ell}}{(2\ell+1)!!} \, s^{\,\ell - n} 
}
\quad \quad (4a)
\qquad\qquad
\boxed{
    \tilde{I}_0^4(s) \; \underset{s \rightarrow 0}{\sim} \; - \frac{\sigma_2}{6 \, s^2}
}
\quad \quad (4b)
```

![All the I_l^n](assets/Iln_terms/all_Iln.png)



The three regimes are visible at a glance: the curves that flatten out
(``\ell = n``), the ones that fall off as a power law (``\ell > n``) and the ones that
blow up (``\ell < n``).

### One by one

|                                              |                                  |
| :------------------------------------------: | :------------------------------: |
|       ![I00](assets/Iln_terms/I00.png)       | ![I20](assets/Iln_terms/I20.png) |
|       ![I40](assets/Iln_terms/I40.png)       | ![I02](assets/Iln_terms/I02.png) |
|       ![I22](assets/Iln_terms/I22.png)       | ![I31](assets/Iln_terms/I31.png) |
|       ![I13](assets/Iln_terms/I13.png)       | ![I11](assets/Iln_terms/I11.png) |
| ![I04_tilde](assets/Iln_terms/I04_tilde.png) |                                  |


## Reproducing the figures

The figures are not built by the documentation: they are committed under
`docs/src/assets/Iln_terms/`, and regenerated on demand by the notebook in the `theory/`
directory. From `theory/`, after the one-time setup described in its `README.md`:

```bash
$ jupyter lab Iln_terms.ipynb
```

Running it top to bottom writes the plots and the tables of numerical values
(`Iln_values.txt`, `Iln_direct_values.txt` and `Iln_large_s_values.txt`) in
`theory/Iln_terms/`, and, since `SAVE_TO_DOCS = true`, also refreshes the copies used by
this page.

See also: [`GaPSE.IPSTools`](@ref), [`GaPSE.IntegralIPS`](@ref),
[`GaPSE.func_I04_tilde`](@ref).
