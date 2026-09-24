## TO-BE-CHECKED-STILL:  The `I_l^n` large-``s`` behaviour

The opposite end has a pleasant surprise: **the exponent is the same for every ``\ell``
and every ``n``**. Getting there, however, needs more care than the ``s \rightarrow 0``
side, because the obvious route does not work.

### The integral this section works with

Everything below is derived for the *idealised* integral

```math
    \bar{I}_\ell^n(s) := \int_0^{+\infty} \frac{\mathrm{d}q}{2\pi^2} \, q^2 \, P(q)
        \, \frac{j_\ell(qs)}{(qs)^n} \; ,
    \quad \quad (3.1)
```

i.e. Eq.(1) with the two cuts sent to ``0`` and ``+\infty``. Unlike the
``s \rightarrow 0`` side — where the cuts *are* the answer, as
[Why the cut cannot be dropped](#why-the-cut-cannot-be-dropped) shows — here the
replacement is harmless, and for a reason that is easy to state: the substitution
``x := qs`` makes it explicit that the weight of the integral sits around
``q \sim 1/s``. For large ``s`` that region is well inside ``[k_\mathrm{min},
k_\mathrm{max}]``, and both the piece below ``k_\mathrm{min}`` (integrand
``\propto q^{\,2+n_P}``, i.e. vanishing) and the piece above ``k_\mathrm{max}``
(argument ``qs \gg 1``, i.e. an oscillation with an amplitude already decaying as
``q^{-\alpha-n-1}``) contribute a vanishing fraction. So
``I_\ell^n(s) \rightarrow \bar{I}_\ell^n(s)`` as ``s`` grows, and the ``\bar{\cdot}``
is dropped from here on.

The two boundaries of that statement are exactly the two listed in
[Where it stops holding](#where-it-stops-holding) below, and they are what the last
section measures.

### The small-``k`` slope of the matter Power Spectrum

The proof that

```math
    P(k) \; \underset{k \rightarrow 0^{+}}{\sim}  \, k^{\,n_P} \; ,
    \quad \quad n_P = n_s \simeq 0.96
```

— together with the large-``k`` tail, the power laws `InputPS` extrapolates with, and
which ``\sigma_i`` survive the removal of the cuts — now lives in its own page,
[The input Power Spectrum](theory_InputPowerSpectrum.md), since it is needed by the
``\Delta\chi \rightarrow 0`` limits as much as by what follows here. Only the result
is used below: the small-``k`` slope of the *matter* Power Spectrum is ``n_P = n_s``,
**not** ``n_s - 1``, which is the slope of the dimensionless primordial curvature
spectrum ``\Delta^2_{\mathcal{R}}``.

### Why the limit cannot be taken inside the integral

The tempting move is to substitute ``q`` with ``x := q s`` in Eq.(3.1),

```math
x:=qs \quad \Rightarrow \quad q = \frac{x}{s} \quad \Rightarrow \quad \mathrm{d}q = \frac{\mathrm{d}x}{s} \quad \quad (3.2)
```

so that

```math
\begin{align*}
     (3.1): \quad \quad I_\ell^n(s)
    &:= \int_0^{+\infty} \frac{\mathrm{d}q}{2\pi^2} \, q^2 \, P(q) \, \frac{j_\ell(qs)}{(qs)^n} \\[10pt]
    \mathrm{inserting \; (3.2)} \; \rightarrow \;\; \quad
    &=  \int_0^{+\infty} \frac{\mathrm{d}x}{2 \pi^2 s} \, \frac{x^2}{s^2} \,
        P\!\left(\frac{x}{s}\right) \frac{j_\ell(x)}{x^n} \\[10pt]
    &= \frac{s^{-3-n}}{2\pi^2} \int_0^{+\infty} \mathrm{d}x \; x^{2-n} \,
        P\!\left(\frac{x}{s}\right) j_\ell(x) \; , \quad \quad (3.3)
\end{align*}
```

and then to replace ``P(x/s)`` by ``A (x/s)^{n_P}`` because ``x/s \rightarrow 0``.
**This is not legitimate.** The integration runs up to ``x = +\infty``, so ``x/s`` is an
indeterminate ``\infty/\infty``: there is no ``s`` beyond which ``x/s`` is small for
*every* ``x``, and no integrable dominating function, so neither dominated nor monotone
convergence applies.

The failure is not academic. Doing it anyway leaves

```math
    \int_0^{+\infty} \mathrm{d}x \; x^{\,\mu-1} \, j_\ell(x) \; ,
    \qquad \mu := 3 + n_P - n \; ,
    \quad \quad (3.4)
```

which converges only for ``-\ell < \mu < 2``. For ``I_0^0``, ``I_2^0``, ``I_4^0``
(``\mu = 3.96``), ``I_3^1`` and ``I_1^1`` (``\mu = 2.96``) it **diverges**: at large ``x``
the integrand behaves as ``x^{\mu-2}\sin(x - \ell\pi/2)``, whose primitive oscillates with
an unboundedly growing amplitude. Five of the eight cases would produce no answer at all.

What follows instead is an argument that never forms that object.

### The proof, by Mellin transform


Write Eq.(3.1) as a Mellin convolution, isolating the ``s^{-n}``:
```math
    (3.1): \quad \quad I_\ell^n(s) := \int_0^{+\infty} \frac{\mathrm{d}q}{2\pi^2} \, q^2 \, P(q) \, \frac{j_\ell(qs)}{(qs)^n} \\[14pt]
    \Rightarrow \quad I_\ell^n(s) = s^{-n} \int_0^{+\infty} \mathrm{d}q \; f(q) \, j_\ell(q s) \; ,
    \qquad f(q) := \frac{q^{\,2-n} \, P(q)}{2\pi^2} \; .
    \quad \quad (3.5)
```

Introduce the two Mellin transforms

```math
    \tilde{f}(z) := \int_0^{+\infty} \mathrm{d}q \; q^{\,z-1} f(q) \; ,
    \qquad
    \mathcal{M}_\ell(z) := \int_0^{+\infty} \mathrm{d}x \; x^{\,z-1} j_\ell(x)
    = \sqrt{\frac{\pi}{2}} \; 2^{\,z-3/2} \;
    \frac{\Gamma\!\left(\frac{\ell+z}{2}\right)}{\Gamma\!\left(\frac{\ell-z+3}{2}\right)} \; ,
    \quad (3.6)
```

the second one being the general case of the known integrals of the
[Spherical Bessel Functions](theory_SphericalBesselFunctions.md) page (``z = 1`` gives back
``\int_0^\infty j_\ell = \sqrt{\pi} \, \Gamma\!\left(\frac{\ell+1}{2}\right) /
2\Gamma\!\left(1+\frac{\ell}{2}\right)``). Inserting the inverse transform of ``j_\ell``
into Eq.(3.5) and exchanging the two integrals — legitimate here, the ``z``-contour being
a vertical line on which everything is absolutely convergent — gives the
Mellin-Parseval representation

```math
\begin{align*}
    I_\ell^n(s) &= s^{-n} \int_0^{+\infty} \mathrm{d}q \; f(q) \,
        \frac{1}{2\pi i}\int_{c - i\infty}^{c + i\infty} \mathrm{d}z \;
        \mathcal{M}_\ell(z) \, (q s)^{-z} \\[10pt]
    &= \frac{s^{-n}}{2\pi i} \int_{c - i\infty}^{c + i\infty} \mathrm{d}z \;
        \mathcal{M}_\ell(z) \, s^{-z} \int_0^{+\infty} \mathrm{d}q \; q^{-z} f(q) \\[10pt]
    &= \frac{s^{-n}}{2\pi i} \int_{c - i\infty}^{c + i\infty} \mathrm{d}z \;
        \mathcal{M}_\ell(z) \, \tilde{f}(1-z) \, s^{-z} \; .
    \quad \quad (3.7)
\end{align*}
```

This is exact: no limit has been taken yet. The large-``s`` behaviour is now read off the
**analytic structure in ``z``**, because ``s^{-z}`` decays faster the further right the
contour sits.

**Poles of ``\mathcal{M}_\ell``.** From Eq.(3.6),
``\Gamma\!\left(\frac{\ell+z}{2}\right)`` has simple poles at ``z = -\ell - 2k``,
``k \geq 0``, while ``1/\Gamma\!\left(\frac{\ell-z+3}{2}\right)`` is entire. So
``\mathcal{M}_\ell(z)`` is **analytic for ``\mathrm{Re}\,z > -\ell``**: in particular at
``z = \mu``, always, for all our cases.

**Poles of ``\tilde{f}(1-z)``.** Here

```math
    \tilde{f}(1-z) = \frac{1}{2\pi^2}\int_0^{+\infty}\mathrm{d}q \; q^{\,2-n-z} \, P(q) \;.
```

Splitting at any scale ``q_0`` inside the power-law region and using ``P(q) \sim A\,q^{\,n_P}`` below it (see [The input Power Spectrum](theory_InputPowerSpectrum.md)),

```math
    \frac{1}{2\pi^2}\int_0^{q_0}\mathrm{d}q \; A \, q^{\,2-n-z+n_P}
    = \frac{A}{2\pi^2} \, \frac{q_0^{\,\mu-z}}{\mu - z} \; ,
```

while ``\int_{q_0}^{+\infty}`` converges and is analytic for
``\mathrm{Re}\,z > 2 - n - \beta``, with ``-\beta`` the large-``q`` slope of ``P``. Hence
``\tilde{f}(1-z)`` has a **simple pole at ``z = \mu``**, of residue ``-A/2\pi^2``, and it
is the *rightmost* singularity: everything further right comes from the subleading
small-``q`` structure of ``P``, i.e. from ``T(k)`` departing from 1.

**Shift the contour.** Take ``c < \mu`` and push it to ``c' > \mu``. Moving a Mellin
contour to the right subtracts the residues crossed,

```math
    \frac{1}{2\pi i}\int_{c} = -\!\!\sum_{c < \mathrm{Re}\, z_k < c'}\!\!
        \mathrm{Res} \; + \; \frac{1}{2\pi i}\int_{c'} \; ,
```

and only ``z = \mu`` is crossed. Since
``\mathrm{Res}_{z=\mu}\left[\mathcal{M}_\ell(z)\tilde{f}(1-z)s^{-z}\right]
= -\frac{A}{2\pi^2}\mathcal{M}_\ell(\mu) \, s^{-\mu}`` and the leftover integral is
``\mathcal{O}(s^{-c'})``, Eq.(3.7) gives

```math
    I_\ell^n(s) = s^{-n}\left[
        \frac{A}{2\pi^2} \, \mathcal{M}_\ell(\mu) \, s^{-\mu}
        + \mathcal{O}\!\left(s^{-c'}\right)\right] \; ,
```

and with ``\mu + n = 3 + n_P``:

```math
\boxed{
    I_\ell^n(s) \; \underset{s \rightarrow +\infty}{\sim} \;
    \frac{A}{2\pi^2} \, \mathcal{M}_\ell(3 + n_P - n) \; s^{-(3+n_P)}
}
\quad \quad (3.8)
```

Three things are worth underlining.

- **The exponent depends on neither ``\ell`` nor ``n``.** The ``s^{-n}`` that the
  ``(qs)^{-n}`` factor pulls out of the integral is exactly compensated by the ``n`` that
  the same factor subtracts from ``\mu``. Only the amplitude ``\mathcal{M}_\ell(\mu)``
  knows about ``\ell`` and ``n``. The limit is always ``0``.
- **``\mathcal{M}_\ell(\mu)`` for ``\mu \geq 2`` is not a fudge.** It is never used as a
  convergent integral: the residue theorem evaluates the *function*
  ``\mathcal{M}_\ell(z)`` at ``z = \mu``, and that function is analytic there. The
  divergence of Eq.(3.4) is an artefact of the naive route, not a property of the answer.
- **The corrections are not a single clean power.** They are governed by the next
  singularities to the right, i.e. by how ``P`` leaves its ``A q^{n_P}`` behaviour — the
  transfer function and, on top of it, the BAO wiggles.

This is also why `IntegralIPS` seeds its right-hand power-law fit with
``p_0 = [-4.0, 1.0]``: with ``n_P = n_s \simeq 0.96``, ``3 + n_P \simeq 4``.

### The check

`theory/Iln_terms.ipynb` reads ``A`` and ``n_P`` out of the `InputPS` left fit (for
`data/WideA_ZA_pk.dat`, ``n_P = 0.960`` and ``A = 3.012 \times 10^6``) and compares
Eq.(3.8) with the stored ``I_\ell^n``. The ratio ``I_\ell^n(s) \, / \,`` Eq.(3.8):

| ``s``      | ``I_0^0`` | ``I_2^0`` | ``I_4^0`` | ``I_0^2`` | ``I_2^2`` | ``I_3^1`` | ``I_1^3`` | ``I_1^1`` |
| :--------- | --------: | --------: | --------: | --------: | --------: | --------: | --------: | --------: |
| ``\mu``    | 3.96      | 3.96      | 3.96      | 1.96      | 1.96      | 2.96      | 0.96      | 2.96      |
| ``300``    | 0.858     | 0.957     | 0.528     | 1.017     | 0.675     | 0.592     | 0.793     | 0.988     |
| ``10^{3}`` | 1.049     | 1.035     | 0.876     | 1.022     | 0.927     | 0.900     | 0.960     | 1.031     |
| ``3 \times 10^{3}`` | 1.024 | 1.015 | 0.978 | 1.006 | 0.990 | 0.984 | 0.995 | 1.011 |
| ``10^{4}`` | 1.005     | 1.003     | 1.000     | 0.993     | 1.001     | 1.000     | 0.998     | 1.003     |

![Convergence to the large-s limits](assets/Iln_terms/ratios_large_s.png)

Better than 0.7% for all eight at ``s = 10^4 \, h_0^{-1}\mathrm{Mpc}``, the five
``\mu > 2`` cases included — which is the numerical confirmation that the analytic
continuation is the right object. The approach is not monotonic and not a single power,
as anticipated above.

### Where it stops holding

Two boundaries, mirroring the small-``s`` ones:

- the pole at ``z = \mu`` is the rightmost singularity only as long as the region
  ``q \sim 1/s`` feeding the integral really is inside the ``P \propto q^{n_P}`` regime.
  It leaves it from below when ``1/s`` drops under ``k_\mathrm{min} = 10^{-5}``, i.e. for
  ``s \gtrsim 10^{5}``, and from above when ``1/s`` climbs past the turnover of ``P``
  around ``k_\mathrm{eq}``, which is why the ratios are still 10% off at ``s = 10^{3}``
  and only settle by ``10^{4}``;
- ``\mathrm{right} = 96466 \, h_0^{-1}\mathrm{Mpc}`` for the ``I_\ell^n``, beyond which an
  `IntegralIPS` is again a power-law extrapolation and not the integral. Note that
  ``\tilde{I}_0^4`` has instead ``\mathrm{right} = 9888``, which a realistic
  ``\Delta\chi \sim 2\chi_\mathrm{max}`` does reach.

### Why ``\tilde{I}_0^4`` is excluded

For ``\tilde{I}_0^4`` one has ``\ell = 0`` and ``n = 4``, hence
``\mu = n_P - 1 \simeq -0.04``, which sits essentially *on* the pole of
``\mathcal{M}_0(z)`` at ``z = 0``. That pole is not an accident: it is precisely the
``\sigma_4 / s^4`` divergence that the ``-1`` of the numerator subtracts away. What the
subtraction leaves behind is the finite part, whose two surviving powers —
``s^{-(3+n_P)}`` from Eq.(3.8) and ``s^{-4}`` from ``-\sigma_4/s^4`` — differ only by
``1 - n_P = 0.04``. They are degenerate for any practical purpose, so ``\tilde{I}_0^4``
never settles onto a clean power law: its measured local slope is still only ``-3.6`` at
``s = 10^{3}`` and ``-3.7`` at ``s = 9 \times 10^{3}``, where its spline already ends.
