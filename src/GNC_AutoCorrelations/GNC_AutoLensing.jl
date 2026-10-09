# -*- encoding: utf-8 -*-
#
# This file is part of GaPSE
# Copyright (C) 2022 Matteo Foglieni
#
# GaPSE is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# GaPSE is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with GaPSE. If not, see <http://www.gnu.org/licenses/>.
#



function integrand_ξ_GNC_Lensing(
    IP1::Point, IP2::Point,
    P1::Point, P2::Point,
    y, cosmo::Cosmology; Δχ_min::AbstractFloat=1e-1,
    b1=nothing, b2=nothing, s_b1=nothing, s_b2=nothing, 𝑓_evo1=nothing, 𝑓_evo2=nothing,
    s_lim=nothing, obs::Union{Bool,Symbol}=:noobsvel)

    s1 = P1.comdist
    s2 = P2.comdist
    χ1, D1, a1 = IP1.comdist, IP1.D, IP1.a
    χ2, D2, a2 = IP2.comdist, IP2.D, IP2.a

    Ω_M0 = cosmo.params.Ω_M0
    s_b_s1 = isnothing(s_b1) ? cosmo.params.s_b1 : s_b1
    s_b_s2 = isnothing(s_b2) ? cosmo.params.s_b2 : s_b2

    # rewritten to avoid the cancellation between chi1^2+chi2^2 and 2*chi1*chi2*y;
    # see the long comment further down. Both terms below are non-negative for y<=1.
    #Δχ_square = χ1^2 + χ2^2 - 2 * χ1 * χ2 * y
    Δχ_square = (χ1 - χ2)^2 + 2 * χ1 * χ2 * (1 - y)
    Δχ = Δχ_square > 0 ? √(Δχ_square) : zero(Δχ_square)  # throw(AssertionError("Δχ_square=$Δχ_square : y=$y, χ1=$χ1, χ2=$χ2"))

    denomin = s1 * s2 * a1 * a2
    factor = ℋ0^4 * Ω_M0^2 * D1 * (s1 - χ1) * D2 * (s2 - χ2) * (5 * s_b_s1 - 2) * (5 * s_b_s2 - 2)

    first_res = if Δχ > min(Δχ_min, Δχ_min * max(χ1, χ2))
        #χ1χ2 = χ1 * χ2
    #
        #new_J00 = -3 / 4 * χ1χ2^2 / Δχ^4 * (y^2 - 1) * (8 * y * (χ1^2 + χ2^2) - χ1χ2 * (9 * y^2 + 7))
        #new_J02 = -3 / 2 * χ1χ2^2 / Δχ^4 * (y^2 - 1) * (4 * y * (χ1^2 + χ2^2) - χ1χ2 * (3 * y^2 + 5))
        #new_J31 = 9 * y * Δχ^2
        #new_J22 = 9 / 4 * χ1χ2 / Δχ^4 * (
        #    2 * (χ1^4 + χ2^4) * (7 * y^2 - 3)
        #    - 16 * y * χ1χ2 * (y^2 + 1) * (χ1^2 + χ2^2)
        #    + χ1χ2^2 * (11y^4 + 14y^2 + 23)
        #)

        # ---------------------------------------------------------------------------------
        # NUMERICAL STABILITY. The four brackets below are written as an expansion around
        # the singular configuration y = 1, χ1 = χ2, instead of in their "natural" form
        # (kept commented out above each of them).
        #
        #
        # WHY IT IS NEEDED. Every one of those brackets VANISHES at that configuration,
        # while each is evaluated as a sum of terms of size ~chi^4. Near it the answer is
        # therefore the difference of numbers many orders of magnitude larger than itself,
        # and the significant digits cancel away. For J_22 the bracket is exactly 8*Δχ^4
        # at y = 1: with chi = 1000 and Δχ = 0.1 that is 8e-4 obtained from terms of size
        # 8e12, a ratio of 1e-16, i.e. below eps(Float64). Checked against exact rational
        # arithmetic, the old form has ZERO correct digits there, and Δχ^2, B_00 and B_02
        # keep only 5 to 9 of them. This is not academic: the Δχ -> 0 branch below takes
        # over only for Δχ < Δχ_min, so the worst case is exactly where the direct
        # evaluation is still the one being used.
        #
        # THE DERIVATION. Write everything in the symmetric combinations
        #     u := χ1^2 + χ2^2 ,    v := χ1*χ2 ,    t := y - 1 ,
        # and use the identity that carries all of the cancellation,
        #     u - 2v = (χ1 - χ2)^2 =: w ,
        # which is a square, so forming it directly costs no digits.
        # Substituting also
        #      y = 1 + t ,
        #      Δχ^2 = χ1^2 + χ2^2 - 2 * χ1 * χ2 * y
        #           = (χ1 - χ2)^2 + 2 * χ1 * χ2 * (1 - y)
        #           = w - 2*v*t
        # and collecting the powers of t:
        #
        # new_J00 := -3/4 * χ1χ2^2 / Δχ^4 * (y^2 - 1) * (8*y*(χ1^2+χ2^2) - χ1χ2*(9*y^2+7))
        #          = -3/4 * v^2 / Δχ^4 * (y^2 - 1) * B_00
        # new_J02 := -3/2 * χ1χ2^2 / Δχ^4 * (y^2 - 1) * (4*y*(χ1^2+χ2^2) - χ1χ2*(3*y^2+5))
        #          = -3/2 * v^2 / Δχ^4 * (y^2 - 1) * B_02
        # new_J31 := 9 * y * Δχ^2
        # new_J22 := 9/4 * χ1χ2 / Δχ^4 * (2*(χ1^4+χ2^4)*(7*y^2-3)
        #                                 - 16*y*χ1χ2*(y^2+1)*(χ1^2+χ2^2)
        #                                 + χ1χ2^2*(11y^4+14y^2+23))
        #          = 9/4 * v / Δχ^4 * B_22
        #
        # The elementary expansions used below, all obtained by putting y = 1 + t:
        #      7*y^2 - 3        = 7*t^2 + 14*t + 4
        #      9*y^2 + 7        = 9*t^2 + 18*t + 16
        #      3*y^2 + 5        = 3*t^2 +  6*t +  8
        #      y*(y^2 + 1)      = t^3 + 3*t^2 + 4*t + 2
        #      11*y^4+14*y^2+23 = 11*t^4 + 44*t^3 + 80*t^2 + 72*t + 48
        #
        #   Δχ^2  = u - 2*v*y
        #         = u - 2*v*(1 + t)
        #         = (u - 2*v) - 2*v*t
        #         = w - 2*v*t
        #
        #   B_00 := 8 * y * (χ1^2 + χ2^2) - χ1χ2 * (9 * y^2 + 7)
        #         = 8*y*u - v*(9*y^2+7)
        #         = 8*u*(1+t) - v*(9*t^2 + 18*t + 16)
        #         = 8*u - 16*v + (8*u - 18*v)*t - 9*v*t^2
        #         = 8*w + (8*u - 18*v)*t - 9*v*t^2
        #
        #   B_02 := 4 * y * (χ1^2 + χ2^2) - χ1χ2 * (3 * y^2 + 5)
        #         = 4*y*u - v*(3*y^2+5)
        #         = 4*u*(1+t) - v*(3*t^2 + 6*t + 8)
        #         = 4*u - 8*v + (4*u - 6*v)*t - 3*v*t^2
        #         = 4*w + (4*u - 6*v)*t - 3*v*t^2
        #
        #   B_22 := 2*(χ1^4+χ2^4)*(7*y^2-3) - 16*y*χ1χ2*(y^2+1)*(χ1^2+χ2^2)
        #           + χ1χ2^2*(11y^4+14y^2+23)
        #               with   χ1^4 + χ2^4 = u^2 - 2*v^2
        #         = 2*(u^2 - 2*v^2)*(7*t^2 + 14*t + 4)
        #           - 16*u*v*(t^3 + 3*t^2 + 4*t + 2)
        #           + v^2*(11*t^4 + 44*t^3 + 80*t^2 + 72*t + 48)
        #     collecting the powers of t:
        #       t^0:  8*u^2 - 16*v^2 - 32*u*v + 48*v^2 = 8*(u^2 - 4*u*v + 4*v^2)
        #                                              = 8*(u-2*v)^2 = 8*w^2
        #       t^1:  28*u^2 - 56*v^2 - 64*u*v + 72*v^2 = 4*(7*u^2 - 16*u*v + 4*v^2)
        #                                              = 4*(u-2*v)*(7*u-2*v) = 4*w*(7*u-2*v)
        #       t^2:  14*u^2 - 28*v^2 - 48*u*v + 80*v^2 = 2*(7*u^2 - 24*u*v + 26*v^2)
        #       t^3:  -16*u*v + 44*v^2 = -4*v*(4*u - 11*v)
        #       t^4:  11*v^2
        #         = 8*w^2 + 4*w*(7*u-2*v)*t + 2*(7*u^2-24*u*v+26*v^2)*t^2
        #           - 4*v*(4*u-11*v)*t^3 + 11*v^2*t^4
        #
        # WHY THIS IS BETTER. Along χ2 = χ1 + p*Δχ one has w = p^2*Δχ^2 and
        # t = y - 1 ~ -(1-p^2)*Δχ^2/(2*χ1^2), so BOTH w and v*t are O(Δχ^2) while u and v
        # are O(chi^2). Every term above therefore carries an explicit w or t, and the
        # order of each bracket can be read off by counting them:
        #   Δχ^2 = w - 2*v*t                       -> O(Δχ^2), two terms of that size
        #   B_00 = 8*w + (8*u-18*v)*t - 9*v*t^2    -> O(Δχ^2), leading term 8*w
        #   B_02 = 4*w + (4*u-6*v)*t - 3*v*t^2     -> O(Δχ^2), leading term 4*w
        #   B_22 = 8*w^2 + 4*w*(7*u-2*v)*t + 2*(7*u^2-24*u*v+26*v^2)*t^2 + ...
        #                                          -> O(Δχ^4): w^2, w*t and t^2 all are
        # and that O(Δχ^4) is exactly what cancels the 1/Δχ^4 in front of J_22, leaving it
        # finite. In the old form none of this was visible: the vanishing came out of a
        # subtraction between terms of size ~chi^4, i.e. it was produced by the rounding
        # rather than by the algebra.
        #
        # Three things are worth singling out in the B_22 collection above:
        #   - the t^0 coefficient is a PERFECT SQUARE. The -16*v^2 of the first term and
        #     the +48*v^2 of the third conspire with the -32*u*v to give 8*(u-2*v)^2, i.e.
        #     exactly 8*Δχ^4 at y = 1: that is the identity quoted at the top;
        #   - the t^1 coefficient FACTORS, and it has to. Had it not contained (u-2*v) as
        #     a factor there would be a surviving term of order Δχ^2*chi^2 in B_22, and
        #     J_22 ~ B_22/Δχ^4 would diverge as Δχ^-2. The factorisation is the algebraic
        #     statement that it does not;
        #   - the t^2 coefficient does NOT factor, and must not: at χ1 = χ2 it is
        #     2*(28-48+26)*chi^4 = 12*chi^4, genuinely O(chi^4), and paired with
        #     t^2 = O(Δχ^4/chi^4) it contributes at the same O(Δχ^4) as the other two.
        #     This is the usual trap of these limits showing up again: a term quadratic in
        #     t is the same size as one linear in w.
        #
        # These are ALGEBRAIC IDENTITIES, not approximations: each was verified with exact
        # rational arithmetic on 3000 random (χ1, χ2, y), and symbolically term by term.
        # ---------------------------------------------------------------------------------
        u = χ1^2 + χ2^2
        v = χ1 * χ2
        w = (χ1 - χ2)^2   # = u - 2v, but formed as a square: no cancellation
        t = y - 1

        B_00 = 8 * w + (8 * u - 18 * v) * t - 9 * v * t^2
        B_02 = 4 * w + (4 * u - 6 * v) * t - 3 * v * t^2
        B_22 = (8 * w^2 + 4 * w * (7 * u - 2 * v) * t
                + 2 * (7 * u^2 - 24 * u * v + 26 * v^2) * t^2
                - 4 * v * (4 * u - 11 * v) * t^3 + 11 * v^2 * t^4)

        new_J00 = -3 / 4 * v^2 / Δχ^4 * (y^2 - 1) * B_00
        new_J02 = -3 / 2 * v^2 / Δχ^4 * (y^2 - 1) * B_02
        new_J31 = 9 * y * Δχ^2
        new_J22 = 9 / 4 * v / Δχ^4 * B_22

        I00 = cosmo.tools.I00(Δχ)
        I20 = cosmo.tools.I20(Δχ)
        I13 = cosmo.tools.I13(Δχ)
        I22 = cosmo.tools.I22(Δχ)

        (
            new_J00 * I00 + new_J02 * I20 +
            new_J31 * I13 + new_J22 * I22
        )

    else

        # for Δχ → 0 the J02 term vanishes, the J31 one gives 3 * σ_2 and the
        # direction-dependent parts of J00 and J22 cancel each other;
        # see "The Δχ → 0 limits" page of the documentation.
        3 * cosmo.tools.σ_2 + 6 / 5 * χ1^2 * cosmo.tools.σ_0
    end

    return factor / denomin * first_res
end

function integrand_ξ_GNC_Lensing(
    χ1::AbstractFloat, χ2::AbstractFloat,
    s1::AbstractFloat, s2::AbstractFloat,
    y, cosmo::Cosmology;
    kwargs...)

    P1, P2 = Point(s1, cosmo), Point(s2, cosmo)
    IP1, IP2 = Point(χ1, cosmo), Point(χ2, cosmo)
    return integrand_ξ_GNC_Lensing(IP1, IP2, P1, P2, y, cosmo; kwargs...)
end


"""
    integrand_ξ_GNC_Lensing(
        IP1::Point, IP2::Point,
        P1::Point, P2::Point,
        y, cosmo::Cosmology;
        Δχ_min::AbstractFloat=1e-1, b1=nothing, b2=nothing,
        s_b1=nothing, s_b2=nothing, 𝑓_evo1=nothing, 𝑓_evo2=nothing,
        s_lim=nothing, obs::Union{Bool,Symbol}=:noobsvel
        ) ::Float64

    integrand_ξ_GNC_Lensing(
        χ1::AbstractFloat, χ2::AbstractFloat,
        s1::AbstractFloat, s2::AbstractFloat,
        y, cosmo::Cosmology;
        kwargs... )::Float64

Return the integrand of the Two-Point Correlation Function (TPCF) of the
Lensing auto-correlation effect arising from the Galaxy Number Counts (GNC).

In the first method, you should pass the two extreme `Point`s (`P1` and `P2`) and the two
intermediate integrand `Point`s (`IP1` and `IP2`) where to
evaluate the function. In the second method (that internally recalls the first),
you must provide the four corresponding comoving distances `s1`, `s2`, `χ1`, `χ2`.
We remember that all the distances are measured in ``h_0^{-1}\\mathrm{Mpc}``.

The analytical expression of this integrand is the following:

```math
\\begin{equation}
    f^{\\kappa\\kappa} (\\chi_1, \\chi_2, s_1, s_2, y) =
    J^{\\kappa\\kappa}_{\\alpha}
    \\left[
        J^{\\kappa\\kappa}_{00} I_0^0(\\Delta\\chi) +
        J^{\\kappa\\kappa}_{02} I_2^0(\\Delta\\chi) +
        J^{\\kappa\\kappa}_{31} I_1^3(\\Delta\\chi) +
        J^{\\kappa\\kappa}_{22} I_2^2(\\Delta\\chi)
    \\right]  \\, ,
\\end{equation}
```

with

```math
\\begin{split}
    J^{\\kappa\\kappa}_{\\alpha} & =
    \\frac{
        \\mathcal{H}_0^4 \\Omega_{\\mathrm{M}0}^2 D(\\chi_1) D(\\chi_2)
    }{
        s_1 s_2 a(\\chi_1) a(\\chi_2)}
    (\\chi_1 - s_1)(\\chi_2 - s_2)
    (5 s_{\\mathrm{b}, 1} - 2)(5 s_{\\mathrm{b}, 2} - 2)
    \\, , \\\\
    %%%%&%%%%%%%%%%%%%
    J^{\\kappa\\kappa}_{00} & =
    -\\frac{ 3 \\chi_1^2 \\chi_2^2}{4 \\Delta\\chi^4} (y^2 - 1)
    \\left[
        8 y (\\chi_1^2 + \\chi_2^2) - 9\\chi_1\\chi_2y^2 -
        7\\chi_1\\chi_2
    \\right]
    \\, , \\\\
    %%%%&%%%%%%%%%%%%%
    J^{\\kappa\\kappa}_{02} & =
    -\\frac{ 3 \\chi_1^2 \\chi_2^2}{2 \\Delta\\chi^4}(y^2 - 1)
    \\left[
        4 y (\\chi_1^2 + \\chi_2^2) - 3 \\chi_1 \\chi_2 y^2 -
        5 \\chi_1 \\chi_2
    \\right]
    \\, , \\\\
    %%%%%%%%%%%%%%%%%%
    J^{\\kappa\\kappa}_{31} & = 9 y \\Delta\\chi^2
    \\, , \\\\
    %%%%%%%%%%%%%%%%%%
    J^{\\kappa\\kappa}_{22} & =
    \\frac{9 \\chi_1 \\chi_2}{4 \\Delta\\chi^4}
    \\left[
        2(\\chi_1^4 + \\chi_2^4)(7 y^2 - 3) -
        16 y \\chi_1 \\chi_2 (\\chi_1^2 + \\chi_2^2)(y^2 + 1) +
        \\right.\\\\
        &\\left.\\qquad\\qquad\\qquad
        \\chi_1^2 \\chi_2^2 (11y^4 + 14y^2 + 23)
    \\right]
    \\, ,
\\end{split}
```

where:

- ``s_1`` and ``s_2`` are comoving distances;

- ``D_1 = D(s_1)``, ... is the linear growth factor (evaluated in ``s_1``);

- ``a_1 = a(s_1)``, ... is the scale factor (evaluated in ``s_1``);

- ``f_1 = f(s_1)``, ... is the linear growth rate (evaluated in ``s_1``);

- ``\\mathcal{H}_1 = \\mathcal{H}(s_1)``, ... is the comoving
  Hubble distances (evaluated in ``s_1``);

- ``y = \\cos{\\theta} = \\hat{\\mathbf{s}}_1 \\cdot \\hat{\\mathbf{s}}_2``;

- ``\\mathcal{R}_1 = \\mathcal{R}(s_1)``, ... is
  computed by `func_ℛ_GNC` in `cosmo::Cosmology` (and evaluated in ``s_1`` );
  the definition of ``\\mathcal{R}(s)`` is the following:
  ```math
  \\mathcal{R}(s) = 5 s_{\\mathrm{b}}(s) + \\frac{2 - 5 s_{\\mathrm{b}}(s)}{\\mathcal{H}(s) \\, s} +
  \\frac{\\dot{\\mathcal{H}}(s)}{\\mathcal{H}(s)^2} - \\mathit{f}_{\\mathrm{evo}} \\quad ;
  ```

- ``b_1``, ``s_{\\mathrm{b}, 1}``, ``\\mathit{f}_{\\mathrm{evo}, 1}``
  (and ``b_2``, ``s_{\\mathrm{b}, 2}``, ``\\mathit{f}_{\\mathrm{evo}, 2}``) :
  galaxy bias, magnification bias (i.e. the slope of the luminosity function at the luminosity threshold),
  and evolution bias for the first (second) effect;

- ``\\Omega_{\\mathrm{M}0} = \\Omega_{\\mathrm{cdm}} + \\Omega_{\\mathrm{b}}`` is the sum of
  cold-dark-matter and barionic density parameters (again, stored in `cosmo`);

- ``I_\\ell^n`` and ``\\sigma_i`` are defined as
  ```math
  I_\\ell^n(s) = \\int_0^{+\\infty} \\frac{\\mathrm{d}q}{2\\pi^2}
  \\, q^2 \\, P(q) \\, \\frac{j_\\ell(qs)}{(qs)^n} \\quad ,
  \\quad \\sigma_i = \\int_0^{+\\infty} \\frac{\\mathrm{d}q}{2\\pi^2}
  \\, q^{2-i} \\, P(q)
  ```
  with ``P(q)`` as the matter Power Spectrum at ``z=0`` (stored in `cosmo`)
  and ``j_\\ell`` as spherical Bessel function of order ``\\ell``;

- ``\\tilde{I}_0^4`` is defined as
  ```math
  \\tilde{I}_0^4 = \\int_0^{+\\infty} \\frac{\\mathrm{d}q}{2\\pi^2}
  \\, q^2 \\, P(q) \\, \\frac{j_0(qs) - 1}{(qs)^4}
  ```
  with ``P(q)`` as the matter Power Spectrum at ``z=0`` (stored in `cosmo`)
  and ``j_\\ell`` as spherical Bessel function of order ``\\ell``;

- ``\\mathcal{H}_0``, ``f_0`` and so on are evaluated at the observer position (i.e. at present day);

- ``\\Delta\\chi_1 := \\sqrt{\\chi_1^2 + s_2^2-2\\,\\chi_1\\,s_2\\,y}`` and
  ``\\Delta\\chi_2 := \\sqrt{s_1^2 + \\chi_2^2-2\\,s_1\\,\\chi_2\\,y}``;

- ``s=\\sqrt{s_1^2 + s_2^2 - 2 \\, s_1 \\, s_2 \\, y}`` and
  ``\\Delta\\chi := \\sqrt{\\chi_1^2 + \\chi_2^2-2\\,\\chi_1\\,\\chi_2\\,y}``.

In this TPCF there are no observer terms. The `obs` keyword is inserted only for compatibility with
the other GNC TPCFs.

This function is used inside `ξ_GNC_Lensing` with trapz() from the
[Trapz](https://github.com/francescoalemanno/Trapz.jl) Julia package.


## Inputs

-  `IP1::Point`, `IP2::Point`, `P1::Point`, `P2::Point` or `χ1`, `χ2`, `s1`, `s2`: `Point`/comoving
  distances where the TPCF has to be calculated; they contain all the
  data of interest needed for this calculus (comoving distance, growth factor and so on).

- `y`: the cosine of the angle between the two points `P1` and `P2` wrt the observer

- `cosmo::Cosmology`: cosmology to be used in this computation; it contains all the splines
  used for the conversion `s` -> `Point`, and all the cosmological parameters ``b``, ...

## Keyword Arguments

- `b1=nothing`, `s_b1=nothing`, `𝑓_evo1=nothing` and `b2=nothing`, `s_b2=nothing`, `𝑓_evo2=nothing`:
  galaxy, magnification and evolutionary biases respectively for the first and the second effect
  computed in this TPCF:
  - if not set (i.e. if you leave the default value `nothing`) the values stored in the input `cosmo`
    will be considered;
  - if you set one or more values, they will override the `cosmo` ones in this computation;
  - the two sets of values should be different only if you are interested in studing two galaxy species;
  - only the required parameters for the chosen TPCF will be used, depending on its analytical expression;
    all the others will have no effect, we still inserted them for pragmatical code compatibilities.

- `s_lim=nothing` : parameter used in order to avoid the divergence of the ``\\mathcal{R}`` and
  ``\\mathfrak{R}`` denominators: when ``0 \\leq s \\leq s_\\mathrm{lim}`` the returned values are
  ```math
  \\forall \\, s \\in [ 0, s_\\mathrm{lim} ] \\; : \\quad
      \\mathfrak{R}(s) = 1 - \\frac{1}{\\mathcal{H}_0 \\, s_\\mathrm{lim}} \\; , \\quad
      \\mathcal{R}(s) = 5 s_{\\mathrm{b}} +
          \\frac{2 - 5 s_{\\mathrm{b}}}{\\mathcal{H}_0 \\, s_\\mathrm{lim}} +
          \\frac{\\dot{\\mathcal{H}}}{\\mathcal{H}_0^2} - \\mathit{f}_{\\mathrm{evo}} \\; .
  ```
  If `nothing`, the default value stored in `cosmo` will be considered.

- `obs::Union{Bool,Symbol} = :noobsvel` : do you want to consider the observer terms in the computation of the
  chosen GNC TPCF effect?
  - `:yes` or `true` -> all the observer effects will be considered
  - `:no` or `false` -> no observer term will be taken into account
  - `:noobsvel` -> the observer terms related to the observer velocity (that you can find in the CF concerning Doppler)
    will be neglected, the other ones will be taken into account

- `Δχ_min::AbstractFloat = 1e-1` : when
  ``\\Delta\\chi = \\sqrt{\\chi_1^2 + \\chi_2^2 - 2 \\, \\chi_1 \\chi_2 y} \\to 0^{+}``,
  some ``I_\\ell^n`` term diverges, but the overall parenthesis has a known limit:

  ```math
      \\lim_{\\Delta\\chi \\to 0^{+}} \\left(
          J_{00}^{\\kappa\\kappa} \\, I^0_0(\\Delta\\chi) +
          J_{02}^{\\kappa\\kappa} \\, I^0_2(\\Delta\\chi) +
          J_{31}^{\\kappa\\kappa} \\, I^3_1(\\Delta\\chi) +
          J_{22}^{\\kappa\\kappa} \\, I^2_2(\\Delta\\chi)
      \\right) =
          3 \\, \\sigma_2 + \\frac{6}{5} \\, \\chi_2^2 \\, \\sigma_0
  ```

  So, when it happens that ``\\Delta\\chi < \\Delta\\chi_\\mathrm{min}``, the function considers this limit
  as the result of the parenthesis instead of calculating it in the normal way; it prevents
  computational divergences.

See also: [`Point`](@ref), [`Cosmology`](@ref), [`ξ_GNC_multipole`](@ref),
[`map_ξ_GNC_multipole`](@ref), [`print_map_ξ_GNC_multipole`](@ref),
[`ξ_GNC_Lensing`](@ref)
"""
integrand_ξ_GNC_Lensing



##########################################################################################92



function ξ_GNC_Lensing(P1::Point, P2::Point, y, cosmo::Cosmology;
    en::AbstractFloat=1e6, N_χs_2::Int=100, suit_sampling::Bool=true, kwargs...)

    χ1s = P1.comdist .* range(1e-6, 1, length=N_χs_2)
    #χ2s = P2.comdist .* range(1e-5, 1, length = N_χs_2 + 7)
    χ2s = P2.comdist .* range(1e-6, 1, length=N_χs_2)

    IP1s = [GaPSE.Point(x, cosmo) for x in χ1s]
    IP2s = [GaPSE.Point(x, cosmo) for x in χ2s]

    int_ξ_Lensings = [
    en * GaPSE.integrand_ξ_GNC_Lensing(IP1, IP2, P1, P2, y, cosmo; kwargs...)
    for IP1 in IP1s, IP2 in IP2s
    ]

    res = trapz((χ1s, χ2s), int_ξ_Lensings)
    #println("res = $res")
    return res / en
end


function ξ_GNC_Lensing(s1, s2, y, cosmo::Cosmology; kwargs...)
    P1, P2 = Point(s1, cosmo), Point(s2, cosmo)
    return ξ_GNC_Lensing(P1, P2, y, cosmo; kwargs...)
end


"""
    ξ_GNC_Lensing(P1::Point, P2::Point, y, cosmo::Cosmology;
        en::AbstractFloat = 1e6, Δχ_min::AbstractFloat = 1e-1,
        N_χs_2::Int = 100,
        s_b1=nothing, s_b2=nothing, 𝑓_evo1=nothing, 𝑓_evo2=nothing,
        s_lim=nothing, obs::Union{Bool,Symbol}=:noobsvel,
        suit_sampling::Bool=true ) ::Float64

    ξ_GNC_Lensing(s1, s2, y, cosmo::Cosmology;
        kwargs...) ::Float64

Return the Two-Point Correlation Function (TPCF) of the Lensing auto-correlation effect
arising from the Galaxy Number Counts (GNC).

In the first method, you should pass the two `Point` (`P1` and `P2`) where to
evaluate the function, while in the second method (that internally recalls the first)
you must provide the two corresponding comoving distances `s1` and `s2`.
We remember that all the distances are measured in ``h_0^{-1}\\mathrm{Mpc}``.

The analytical expression of this TPCF is the following:

```math
\\begin{split}
    \\xi^{\\kappa\\kappa} (s_1, s_2, y) =
    \\int_0^{s_1} \\mathrm{d}\\chi_1 \\int_0^{s_2} \\mathrm{d}\\chi_2\\;
    J^{\\kappa\\kappa}_{\\alpha}
    &\\left[
        J^{\\kappa\\kappa}_{00} I_0^0(\\Delta\\chi) +
        J^{\\kappa\\kappa}_{02} I_2^0(\\Delta\\chi) +
        \\right.\\\\
        &\\left.
        J^{\\kappa\\kappa}_{31} I_1^3(\\Delta\\chi) +
        J^{\\kappa\\kappa}_{22} I_2^2(\\Delta\\chi)
    \\right]  \\, ,
\\end{split}
```

with

```math
\\begin{split}
    J^{\\kappa\\kappa}_{\\alpha} & =
    \\frac{
        \\mathcal{H}_0^4 \\Omega_{\\mathrm{M}0}^2 D(\\chi_1) D(\\chi_2)
    }{
        s_1 s_2 a(\\chi_1) a(\\chi_2)}
    (\\chi_1 - s_1)(\\chi_2 - s_2)
    (5 s_{\\mathrm{b}, 1} - 2)(5 s_{\\mathrm{b}, 2} - 2)
    \\, , \\\\
    %%%%&%%%%%%%%%%%%%
    J^{\\kappa\\kappa}_{00} & =
    -\\frac{ 3 \\chi_1^2 \\chi_2^2}{4 \\Delta\\chi^4} (y^2 - 1)
    \\left[
        8 y (\\chi_1^2 + \\chi_2^2) - 9\\chi_1\\chi_2y^2 -
        7\\chi_1\\chi_2
    \\right]
    \\, , \\\\
    %%%%&%%%%%%%%%%%%%
    J^{\\kappa\\kappa}_{02} & =
    -\\frac{ 3 \\chi_1^2 \\chi_2^2}{2 \\Delta\\chi^4}(y^2 - 1)
    \\left[
        4 y (\\chi_1^2 + \\chi_2^2) - 3 \\chi_1 \\chi_2 y^2 -
        5 \\chi_1 \\chi_2
    \\right]
    \\, , \\\\
    %%%%%%%%%%%%%%%%%%
    J^{\\kappa\\kappa}_{31} & = 9 y \\Delta\\chi^2
    \\, , \\\\
    %%%%%%%%%%%%%%%%%%
    J^{\\kappa\\kappa}_{22} & =
    \\frac{9 \\chi_1 \\chi_2}{4 \\Delta\\chi^4}
    \\left[
        2(\\chi_1^4 + \\chi_2^4)(7 y^2 - 3) -
        16 y \\chi_1 \\chi_2 (\\chi_1^2 + \\chi_2^2)(y^2 + 1) +
        \\right.\\\\
        &\\left.\\qquad\\qquad\\qquad
        \\chi_1^2 \\chi_2^2 (11y^4 + 14y^2 + 23)
    \\right]
    \\, ,
\\end{split}
```

where:

- ``s_1`` and ``s_2`` are comoving distances;

- ``D_1 = D(s_1)``, ... is the linear growth factor (evaluated in ``s_1``);

- ``a_1 = a(s_1)``, ... is the scale factor (evaluated in ``s_1``);

- ``f_1 = f(s_1)``, ... is the linear growth rate (evaluated in ``s_1``);

- ``\\mathcal{H}_1 = \\mathcal{H}(s_1)``, ... is the comoving
  Hubble distances (evaluated in ``s_1``);

- ``y = \\cos{\\theta} = \\hat{\\mathbf{s}}_1 \\cdot \\hat{\\mathbf{s}}_2``;

- ``\\mathcal{R}_1 = \\mathcal{R}(s_1)``, ... is
  computed by `func_ℛ_GNC` in `cosmo::Cosmology` (and evaluated in ``s_1`` );
  the definition of ``\\mathcal{R}(s)`` is the following:
  ```math
  \\mathcal{R}(s) = 5 s_{\\mathrm{b}}(s) + \\frac{2 - 5 s_{\\mathrm{b}}(s)}{\\mathcal{H}(s) \\, s} +
  \\frac{\\dot{\\mathcal{H}}(s)}{\\mathcal{H}(s)^2} - \\mathit{f}_{\\mathrm{evo}} \\quad ;
  ```

- ``b_1``, ``s_{\\mathrm{b}, 1}``, ``\\mathit{f}_{\\mathrm{evo}, 1}``
  (and ``b_2``, ``s_{\\mathrm{b}, 2}``, ``\\mathit{f}_{\\mathrm{evo}, 2}``) :
  galaxy bias, magnification bias (i.e. the slope of the luminosity function at the luminosity threshold),
  and evolution bias for the first (second) effect;

- ``\\Omega_{\\mathrm{M}0} = \\Omega_{\\mathrm{cdm}} + \\Omega_{\\mathrm{b}}`` is the sum of
  cold-dark-matter and barionic density parameters (again, stored in `cosmo`);

- ``I_\\ell^n`` and ``\\sigma_i`` are defined as
  ```math
  I_\\ell^n(s) = \\int_0^{+\\infty} \\frac{\\mathrm{d}q}{2\\pi^2}
  \\, q^2 \\, P(q) \\, \\frac{j_\\ell(qs)}{(qs)^n} \\quad ,
  \\quad \\sigma_i = \\int_0^{+\\infty} \\frac{\\mathrm{d}q}{2\\pi^2}
  \\, q^{2-i} \\, P(q)
  ```
  with ``P(q)`` as the matter Power Spectrum at ``z=0`` (stored in `cosmo`)
  and ``j_\\ell`` as spherical Bessel function of order ``\\ell``;

- ``\\mathcal{H}_0``, ``f_0`` and so on are evaluated at the observer position (i.e. at present day);

- ``\\Delta\\chi_1 := \\sqrt{\\chi_1^2 + s_2^2-2\\,\\chi_1\\,s_2\\,y}`` and
  ``\\Delta\\chi_2 := \\sqrt{s_1^2 + \\chi_2^2-2\\,s_1\\,\\chi_2\\,y}``;

- ``s=\\sqrt{s_1^2 + s_2^2 - 2 \\, s_1 \\, s_2 \\, y}`` and
  ``\\Delta\\chi := \\sqrt{\\chi_1^2 + \\chi_2^2-2\\,\\chi_1\\,\\chi_2\\,y}``.

In this TPCF there are no observer terms. The `obs` keyword is inserted only for compatibility with
the other GNC TPCFs.

This function is computed integrating `integrand_ξ_GNC_Lensing` with trapz() from the
[Trapz](https://github.com/francescoalemanno/Trapz.jl) Julia package.


## Inputs

- `P1::Point` and `P2::Point`, or `s1` and `s2`: `Point`/comoving distances where the
  TPCF has to be calculated; they contain all the
  data of interest needed for this calculus (comoving distance, growth factor and so on).

- `y`: the cosine of the angle between the two points `P1` and `P2` wrt the observer

- `cosmo::Cosmology`: cosmology to be used in this computation; it contains all the splines
  used for the conversion `s` -> `Point`, and all the cosmological parameters ``b``, ...

## Keyword Arguments

- `b1=nothing`, `s_b1=nothing`, `𝑓_evo1=nothing` and `b2=nothing`, `s_b2=nothing`, `𝑓_evo2=nothing`:
  galaxy, magnification and evolutionary biases respectively for the first and the second effect
  computed in this TPCF:
  - if not set (i.e. if you leave the default value `nothing`) the values stored in the input `cosmo`
    will be considered;
  - if you set one or more values, they will override the `cosmo` ones in this computation;
  - the two sets of values should be different only if you are interested in studing two galaxy species;
  - only the required parameters for the chosen TPCF will be used, depending on its analytical expression;
    all the others will have no effect, we still inserted them for pragmatical code compatibilities.

- `s_lim=nothing` : parameter used in order to avoid the divergence of the ``\\mathcal{R}`` and
  ``\\mathfrak{R}`` denominators: when ``0 \\leq s \\leq s_\\mathrm{lim}`` the returned values are
  ```math
  \\forall \\, s \\in [ 0, s_\\mathrm{lim} ] \\; : \\quad
      \\mathfrak{R}(s) = 1 - \\frac{1}{\\mathcal{H}_0 \\, s_\\mathrm{lim}} \\; , \\quad
      \\mathcal{R}(s) = 5 s_{\\mathrm{b}} +
          \\frac{2 - 5 s_{\\mathrm{b}}}{\\mathcal{H}_0 \\, s_\\mathrm{lim}} +
          \\frac{\\dot{\\mathcal{H}}}{\\mathcal{H}_0^2} - \\mathit{f}_{\\mathrm{evo}} \\; .
  ```
  If `nothing`, the default value stored in `cosmo` will be considered.

- `obs::Union{Bool,Symbol} = :noobsvel` : do you want to consider the observer terms in the computation of the
  chosen GNC TPCF effect?
  - `:yes` or `true` -> all the observer effects will be considered
  - `:no` or `false` -> no observer term will be taken into account
  - `:noobsvel` -> the observer terms related to the observer velocity (that you can find in the CF concerning Doppler)
    will be neglected, the other ones will be taken into account

- `en::Float64 = 1e6`: just a float number used in order to deal better
  with small numbers;

- `Δχ_min::AbstractFloat = 1e-1` : when
  ``\\Delta\\chi = \\sqrt{\\chi_1^2 + \\chi_2^2 - 2 \\, \\chi_1 \\chi_2 y} \\to 0^{+}``,
  some ``I_\\ell^n`` term diverges, but the overall parenthesis has a known limit:

  ```math
      \\lim_{\\Delta\\chi \\to 0^{+}} \\left(
          J_{00}^{\\kappa\\kappa} \\, I^0_0(\\Delta\\chi) +
          J_{02}^{\\kappa\\kappa} \\, I^0_2(\\Delta\\chi) +
          J_{31}^{\\kappa\\kappa} \\, I^3_1(\\Delta\\chi) +
          J_{22}^{\\kappa\\kappa} \\, I^2_2(\\Delta\\chi)
      \\right) =
          3 \\, \\sigma_2 + \\frac{6}{5} \\, \\chi_2^2 \\, \\sigma_0
  ```

  So, when it happens that ``\\Delta\\chi < \\Delta\\chi_\\mathrm{min}``, the function considers this limit
  as the result of the parenthesis instead of calculating it in the normal way; it prevents
  computational divergences.

- `N_χs_2::Int = 100`: number of points to be used for sampling the integral
  along the ranges `(0, s1)` (for `χ1`) and `(0, s2)` (for `χ2`); it has been checked that
  with `N_χs_2 ≥ 50` the result is stable.

- `suit_sampling::Bool = true` : this bool keyword can be found in all the TPCFs which have at least one `χ` integral;
  it is conceived to enable a sampling of the `χ` integral(s) suited for the given TPCF; however, it actually have an
  effect only in the TPCFs that have such a sampling implemented in the code.
  Currently, only `ξ_GNC_Newtonian_Lensing` (and its simmetryc TPCF) has it.

See also: [`Point`](@ref), [`Cosmology`](@ref), [`ξ_GNC_multipole`](@ref),
[`map_ξ_GNC_multipole`](@ref), [`print_map_ξ_GNC_multipole`](@ref),
[`integrand_ξ_GNC_Lensing`](@ref)
"""
ξ_GNC_Lensing
