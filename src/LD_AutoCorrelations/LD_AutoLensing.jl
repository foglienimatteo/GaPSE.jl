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


function integrand_ξ_LD_Lensing(
    IP1::Point, IP2::Point,
    P1::Point, P2::Point,
    y, cosmo::Cosmology;
    Δχ_min::Float64 = 1e-4)

    s1 = P1.comdist
    s2 = P2.comdist
    χ1, D1, a_χ1 = IP1.comdist, IP1.D, IP1.a
    χ2, D2, a_χ2 = IP2.comdist, IP2.D, IP2.a
    Ω_M0 = cosmo.params.Ω_M0

    Δχ_square = χ1^2 + χ2^2 - 2 * χ1 * χ2 * y
    Δχ = Δχ_square > 0 ? √(Δχ_square) : 0.0
    
    denomin = s1 * s2 * a_χ1 * a_χ2
    factor = ℋ0^4 * Ω_M0^2 * D1 * abs(s1 - χ1) * D2 * abs(s2 - χ2)

    first_res = if Δχ > Δχ_min
        χ1χ2 = χ1 * χ2

        new_J00 = 0.75 * χ1χ2^2 / Δχ^4 * abs(1.0 - y^2) * (8 * y * (χ1^2 + χ2^2) - χ1χ2 * (9 * y^2 + 7))
        new_J02 = 1.5 * χ1χ2^2 / Δχ^4 * abs(1.0 - y^2) * (4 * y * (χ1^2 + χ2^2) - χ1χ2 * (3 * y^2 + 5))
        new_J31 = 9 * y * Δχ^2
        new_J22 = 2.25 * χ1χ2 / Δχ^4 * (
            2 * (χ1^4 + χ2^4) * (7 * y^2 - 3)
            - 16 * y * χ1χ2 * (y^2 + 1) * (χ1^2 + χ2^2)
            + χ1χ2^2 * (11y^4 + 14y^2 + 23)
        )

        I00 = cosmo.tools.I00(Δχ)
        I20 = cosmo.tools.I20(Δχ)
        I13 = cosmo.tools.I13(Δχ)
        I22 = cosmo.tools.I22(Δχ)

        resss = (
            new_J00 * I00 + new_J02 * I20 +
            new_J31 * I13 + new_J22 * I22
        )

        resss
    else

        lim = 4.0 / 15.0 * (5.0 * cosmo.tools.σ_2 + 6.0 * cosmo.tools.σ_0 * χ2^2)
        9.0 / 4.0 * lim
    end

    res = factor / denomin * first_res

    return res
end

function integrand_ξ_LD_Lensing(
    χ1::Float64, χ2::Float64,
    s1::Float64, s2::Float64,
    y, cosmo::Cosmology;
    kwargs...)

    P1, P2 = Point(s1, cosmo), Point(s2, cosmo)
    IP1, IP2 = Point(χ1, cosmo), Point(χ2, cosmo)
    return integrand_ξ_LD_Lensing(IP1, IP2, P1, P2, y, cosmo; kwargs...)
end


"""
    integrand_ξ_LD_Lensing(
        IP1::Point, IP2::Point,
        P1::Point, P2::Point,
        y, cosmo::Cosmology;
        Δχ_min::Float64 = 1e-4) ::Float64

    integrand_ξ_LD_Lensing(
        χ1::Float64, χ2::Float64,
        s1::Float64, s2::Float64,
        y, cosmo::Cosmology; kwargs...) ::Float64

Return the integrand of the Two-Point Correlation Function (TPCF) of the Lensing 
auto-correlation effect arising from the Luminosity Distance (LD) perturbations.

In the first method, you should pass the two extreme `Point`s (`P1` and `P2`) and the two 
intermediate integrand `Point`s (`IP1` and `IP2`) where to 
evaluate the function. In the second method (that internally recalls the first),
you must provide the four corresponding comoving distances `s1`, `s2`, `χ1`, `χ2`.
We remember that all the distances are measured in ``h_0^{-1}\\mathrm{Mpc}``.

The analytical expression of this term is the following:

```math
\\begin{split}
    f^{\\kappa\\kappa} (\\chi_1, \\chi_2, s_1, s_2, y) =  
    \\mathcal{J}^{\\kappa\\kappa}_{\\alpha} \\left[
        \\mathcal{J}^{\\kappa\\kappa}_{00}I^0_0(\\Delta\\chi) + 
        \\mathcal{J}^{\\kappa\\kappa}_{02} I^0_2(\\Delta\\chi) +
        \\mathcal{J}^{\\kappa\\kappa}_{31}I^3_1(\\Delta\\chi) +
        \\mathcal{J}^{\\kappa\\kappa}_{22}I^2_2(\\Delta\\chi)
    \\right] \\nonumber \\, ,
\\end{split}
```
with

```math
\\begin{split}
    \\mathcal{J}^{\\kappa\\kappa}_{\\alpha} & = 
    \\frac{
        \\mathcal{H}_0^4 \\Omega_{\\mathrm{M}0}^2 D(\\chi_1) D(\\chi_2)
    }{
        s_1 s_2 a(\\chi_1) a(\\chi_2)
    }(\\chi_1 - s_1)(\\chi_2 - s_2)
    \\, , \\\\
    %%%%%%%%%%%%%%%%%%%%%%%%
    \\mathcal{J}^{\\kappa\\kappa}_{00} & = 
    -\\frac{ 3 \\chi_1^2 \\chi_2^2}{4 \\Delta\\chi^4} (y^2 - 1)
    \\left[
        8 y (\\chi_1^2 + \\chi_2^2) - 9\\chi_1\\chi_2y^2 - 
        7\\chi_1\\chi_2
    \\right] 
    \\, , \\\\
    %%%%%%%%%%%%%%%%%%%%%%%%
    \\mathcal{J}^{\\kappa\\kappa}_{02} & = 
    -\\frac{ 3 \\chi_1^2 \\chi_2^2}{2 \\Delta\\chi^4}(y^2 - 1)
    \\left[
        4 y (\\chi_1^2 + \\chi_2^2) - 3 \\chi_1 \\chi_2 y^2 -
        5 \\chi_1 \\chi_2
    \\right] 
    \\, , \\\\
    %%%%%%%%%%%%%%%%%%%%%%%%
    \\mathcal{J}^{\\kappa\\kappa}_{31} & = 9 y \\Delta\\chi^2
    \\, , \\\\
    %%%%%%%%%%%%%%%%%%%%%%%%
    \\mathcal{J}^{\\kappa\\kappa}_{22} & = 
    \\frac{9 \\chi_1 \\chi_2}{4 \\Delta\\chi^4}
    \\left[
        2(\\chi_1^4 + \\chi_2^4)(7 y^2 - 3) - 
        16 y \\chi_1 \\chi_2 (\\chi_1^2 + \\chi_2^2)(y^2 + 1) + 
        \\right.\\\\
        &\\left.\\qquad\\qquad\\qquad
        \\chi_1^2 \\chi_2^2 (11y^4 + 14y^2 + 23) 
    \\right] \\nonumber
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

- ``\\mathfrak{R}_1 = \\mathfrak{R}(s_1)``, ... is 
  computed by `func_ℛ_LD` in `cosmo::Cosmology` (and evaluated in ``s_1`` );
  the definition of ``\\mathcal{R}(s)`` is the following:
  ```math
  \\mathfrak{R}(s) = 1 - \\frac{1}{\\mathcal{H}(s) s} ;
  ```

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

This function is used inside `ξ_LD_Lensing` with trapz() from the 
[Trapz](https://github.com/francescoalemanno/Trapz.jl) Julia package.

## Inputs

- `IP1::Point`, `IP2::Point`, `P1::Point`, `P2::Point` or `χ1`, `χ2`, `s1`, `s2`: `Point`/comoving 
  distances where the TPCF has to be calculated; they contain all the 
  data of interest needed for this calculus (comoving distance, growth factor and so on).
  
- `y`: the cosine of the angle between the two points `P1` and `P2` wrt the observer

- `cosmo::Cosmology`: cosmology to be used in this computation; it contains all the splines
  used for the conversion `s` -> `Point`, and all the cosmological parameters ``\\Omega_{\\mathrm{M}0}``, ...

## Keyword arguments

- `Δχ_min::Float64 = 1e-4` : when ``\\Delta\\chi = \\sqrt{\\chi_1^2 + \\chi_2^2 - 2 \\, \\chi_1 \\chi_2 y} \\to 0^{+}``,
  some ``I_\\ell^n`` term diverges, but the overall parenthesis has a known limit:

  ```math
  \\lim_{\\chi\\to 0^{+}} \\left(J^{\\kappa\\kappa}_{00} \\, I^0_0(\\Delta\\chi) + 
        J^{\\kappa\\kappa}_{02} \\, I^0_2(\\Delta\\chi) + 
        J^{\\kappa\\kappa}_{31} \\, I^3_1(\\Delta\\chi) + J^{\\kappa\\kappa}_{22} \\, I^2_2(\\Delta\\chi)
        \\right) = 
        \\frac{4}{15} \\, \\left(5 \\, \\sigma_2 + \\frac{2}{3} \\, σ_0 \\,s_1^2 \\, \\chi_2^2\\right)
  ```

See also: [`Point`](@ref), [`Cosmology`](@ref), [`ξ_LD_multipole`](@ref), 
[`map_ξ_LD_multipole`](@ref), [`print_map_ξ_LD_multipole`](@ref)
"""
integrand_ξ_LD_Lensing




##########################################################################################92





function ξ_LD_Lensing(P1::Point, P2::Point, y, cosmo::Cosmology;
    en::Float64 = 1e6, N_χs_2::Int = 100, Δχ_min::Float64 = 1e-4)

    χ1s = P1.comdist .* range(1e-6, 1.0, length = N_χs_2)
    χ2s = P2.comdist .* range(1e-6, 1.0, length = N_χs_2 + 7)

    IP1s = [GaPSE.Point(x, cosmo) for x in χ1s]
    IP2s = [GaPSE.Point(x, cosmo) for x in χ2s]

    int_ξ_Lensings = [
        en * GaPSE.integrand_ξ_LD_Lensing(IP1, IP2, P1, P2, y, cosmo; Δχ_min = Δχ_min)
        for IP1 in IP1s, IP2 in IP2s
    ]

    res = trapz((χ1s, χ2s), int_ξ_Lensings)
    #println("res = $res")
    return res / en
end


function ξ_LD_Lensing(s1, s2, y, cosmo::Cosmology; kwargs...)
    P1, P2 = Point(s1, cosmo), Point(s2, cosmo)
    return ξ_LD_Lensing(P1, P2, y, cosmo; kwargs...)
end



"""
    ξ_LD_Lensing(P1::Point, P2::Point, y, cosmo::Cosmology;
        en::Float64 = 1e6, Δχ_min::Float64 = 1e-3,
        N_χs_2::Int = 100) ::Float64

    ξ_LD_Lensing(s1, s2, y, cosmo::Cosmology; kwargs...) ::Float64

Return the Two-Point Correlation Function (TPCF) of the Lensing 
auto-correlation effect arising from the Luminosity Distance (LD) perturbations.

In the first method, you should pass the two `Point` (`P1` and `P2`) where to 
evaluate the function, while in the second method (that internally recalls the first) 
you must provide the two corresponding comoving distances `s1` and `s2`.
We remember that all the distances are measured in ``h_0^{-1}\\mathrm{Mpc}``.

The analytical expression of this term is the following:

```math
\\begin{split}
    \\xi^{\\kappa\\kappa} (s_1, s_2, y) = 
    \\int_0^{s_1} \\mathrm{d}\\chi_1 \\int_0^{s_2} \\mathrm{d}\\chi_2 \\;
    \\mathcal{J}^{\\kappa\\kappa}_{\\alpha} \\left[
        \\mathcal{J}^{\\kappa\\kappa}_{00}I^0_0(\\Delta\\chi) + 
        \\mathcal{J}^{\\kappa\\kappa}_{02} I^0_2(\\Delta\\chi) +
        \\mathcal{J}^{\\kappa\\kappa}_{31}I^3_1(\\Delta\\chi) +
        \\mathcal{J}^{\\kappa\\kappa}_{22}I^2_2(\\Delta\\chi)
    \\right] \\nonumber \\, ,
\\end{split}
```
with

```math
\\begin{split}
    \\mathcal{J}^{\\kappa\\kappa}_{\\alpha} & = 
    \\frac{
        \\mathcal{H}_0^4 \\Omega_{\\mathrm{M}0}^2 D(\\chi_1) D(\\chi_2)
    }{
        s_1 s_2 a(\\chi_1) a(\\chi_2)
    }(\\chi_1 - s_1)(\\chi_2 - s_2)
    \\, , \\\\
    %%%%%%%%%%%%%%%%%%%%%%%%
    \\mathcal{J}^{\\kappa\\kappa}_{00} & = 
    -\\frac{ 3 \\chi_1^2 \\chi_2^2}{4 \\Delta\\chi^4} (y^2 - 1)
    \\left[
        8 y (\\chi_1^2 + \\chi_2^2) - 9\\chi_1\\chi_2y^2 - 
        7\\chi_1\\chi_2
    \\right] 
    \\, , \\\\
    %%%%%%%%%%%%%%%%%%%%%%%%
    \\mathcal{J}^{\\kappa\\kappa}_{02} & = 
    -\\frac{ 3 \\chi_1^2 \\chi_2^2}{2 \\Delta\\chi^4}(y^2 - 1)
    \\left[
        4 y (\\chi_1^2 + \\chi_2^2) - 3 \\chi_1 \\chi_2 y^2 -
        5 \\chi_1 \\chi_2
    \\right] 
    \\, , \\\\
    %%%%%%%%%%%%%%%%%%%%%%%%
    \\mathcal{J}^{\\kappa\\kappa}_{31} & = 9 y \\Delta\\chi^2
    \\, , \\\\
    %%%%%%%%%%%%%%%%%%%%%%%%
    \\mathcal{J}^{\\kappa\\kappa}_{22} & = 
    \\frac{9 \\chi_1 \\chi_2}{4 \\Delta\\chi^4}
    \\left[
        2(\\chi_1^4 + \\chi_2^4)(7 y^2 - 3) - 
        16 y \\chi_1 \\chi_2 (\\chi_1^2 + \\chi_2^2)(y^2 + 1) + 
        \\right.\\\\
        &\\left.\\qquad\\qquad\\qquad
        \\chi_1^2 \\chi_2^2 (11y^4 + 14y^2 + 23) 
    \\right] \\nonumber
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

- ``\\mathfrak{R}_1 = \\mathfrak{R}(s_1)``, ... is 
  computed by `func_ℛ_LD` in `cosmo::Cosmology` (and evaluated in ``s_1`` );
  the definition of ``\\mathcal{R}(s)`` is the following:
  ```math
  \\mathfrak{R}(s) = 1 - \\frac{1}{\\mathcal{H}(s) s} ;
  ```

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

This function is computed integrating `integrand_ξ_LD_Lensing` with trapz() from the 
[Trapz](https://github.com/francescoalemanno/Trapz.jl) Julia package.

## Inputs

- `P1::Point` and `P2::Point`, or `s1` and `s2`: `Point`/comoving distances where the 
  TPCF has to be calculated; they contain all the 
  data of interest needed for this calculus (comoving distance, growth factor and so on).
  
- `y`: the cosine of the angle between the two points `P1` and `P2` wrt the observer

- `cosmo::Cosmology`: cosmology to be used in this computation; it contains all the splines
  used for the conversion `s` -> `Point`, and all the cosmological parameters ``\\Omega_{\\mathrm{M}0}``, ...

## Keyword arguments

- `en::Float64 = 1e6`: just a float number used in order to deal better 
  with small numbers;

- `N_χs_2::Int = 100`: number of points to be used for sampling the integral
  along the ranges `(0, s1)` (for `χ1`) and `(0, s2)` (for `χ2`); it has been checked that
  with `N_χs_2 ≥ 50` the result is stable.

- `Δχ_min::Float64 = 1e-4` : when ``\\Delta\\chi = \\sqrt{\\chi_1^2 + \\chi_2^2 - 2 \\, \\chi_1 \\chi_2 y} \\to 0^{+}``,
  some ``I_\\ell^n`` term diverges, but the overall parenthesis has a known limit:

  ```math
  \\lim_{\\Delta\\chi\\to 0^{+}} \\left(J^{\\kappa\\kappa}_{00} \\, I^0_0(\\Delta\\chi) + 
        J^{\\kappa\\kappa}_{02} \\, I^0_2(\\Delta\\chi) + 
        J^{\\kappa\\kappa}_{31} \\, I^3_1(\\Delta\\chi) + J^{\\kappa\\kappa}_{22} \\, I^2_2(\\Delta\\chi)
        \\right) = 
        \\frac{4}{15} \\, \\left(5 \\, \\sigma_2 + \\frac{2}{3} \\, σ_0 \\,s_1^2 \\, \\chi_2^2\\right)
  ```

See also: [`Point`](@ref), [`Cosmology`](@ref), [`ξ_LD_multipole`](@ref), 
[`map_ξ_LD_multipole`](@ref), [`print_map_ξ_LD_multipole`](@ref)
"""
ξ_LD_Lensing
