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




function integrand_ξ_GNC_Lensing_IntegratedGP(
     IP1::Point, IP2::Point,
     P1::Point, P2::Point,
     y, cosmo::Cosmology; obs::Union{Bool,Symbol}=:noobsvel)

     s1 = P1.comdist
     s2, ℛ_s2 = P2.comdist, P2.ℛ_GNC
     χ1, D1, a1 = IP1.comdist, IP1.D, IP1.a
     χ2, D2, a2, f2, ℋ2 = IP2.comdist, IP2.D, IP2.a, IP2.f, IP2.ℋ
     s_b_s2 = cosmo.params.s_b
     Ω_M0 = cosmo.params.Ω_M0

     Δχ_square = χ1^2 + χ2^2 - 2 * χ1 * χ2 * y
     Δχ = √(Δχ_square) > 0 ? √(Δχ_square) : 0

     denomin = a1 * a2 * s1 * s2
     common = 9 * χ2 * ℋ0^4 * Ω_M0^2 * D1 * (χ1 - s1) * D2 * (5 * s_b_s2 - 2)
     parenth = ℋ2 * ℛ_s2 * s2 * (f2 - 1) - 5 * s_b_s2 + 2

     new_J31 = y * Δχ^2
     new_J22 = χ1 * χ2 * (y^2 - 1) / 2

     I13 = cosmo.tools.I13(Δχ)
     I22 = cosmo.tools.I22(Δχ)

     return common * parenth * (new_J22 * I22 + new_J31 * I13) / denomin
end


function integrand_ξ_GNC_Lensing_IntegratedGP(
     χ1::Float64, χ2::Float64,
     s1::Float64, s2::Float64,
     y, cosmo::Cosmology;
     kwargs...)

     P1, P2 = Point(s1, cosmo), Point(s2, cosmo)
     IP1, IP2 = Point(χ1, cosmo), Point(χ2, cosmo)
     return integrand_ξ_GNC_Lensing_IntegratedGP(IP1, IP2, P1, P2, y, cosmo; kwargs...)
end




"""
     integrand_ξ_GNC_Lensing_IntegratedGP(
          IP1::Point, IP2::Point,
          P1::Point, P2::Point,
          y, cosmo::Cosmology; 
          obs::Union{Bool,Symbol}=:noobsvel
          ) ::Float64

     integrand_ξ_GNC_Lensing_IntegratedGP(
          χ1::Float64, χ2::Float64,
          s1::Float64, s2::Float64,
          y, cosmo::Cosmology;
          kwargs...) ::Float64

Return the integrand of the Two-Point Correlation Function (TPCF) given 
by the cross correlation between the Lensing
and the Integrated Gravitational Potential (GP) effects arising 
from the Galaxy Number Counts (GNC).

In the first method, you should pass the two extreme `Point`s (`P1` and `P2`) and the two 
intermediate integrand `Point`s (`IP1` and `IP2`) where to 
evaluate the function. In the second method (that internally recalls the first),
you must provide the four corresponding comoving distances `s1`, `s2`, `χ1`, `χ2`.
We remember that all the distances are measured in ``h_0^{-1}\\mathrm{Mpc}``.

The analytical expression of this integrand is the following:

```math
\\begin{split}
    f^{\\kappa \\int\\!\\phi} (\\chi_1, \\chi_2, s_1 , s_2, y ) = 
    J_{\\alpha}^{\\kappa \\int\\!\\phi} 
    \\left[ 
        J_{31}^{\\kappa \\int\\!\\phi} I_1^3 ( \\Delta \\chi ) +
        J_{22}^{\\kappa \\int\\!\\phi} I_2^2 ( \\Delta \\chi ) 
     \\right] \\, ,
\\end{split}
```

where

```math
\\begin{split}
     J_{\\alpha}^{\\kappa \\int\\!\\phi} &=
    \\frac{
        9 \\chi_2 \\ \\mathcal{H}_0^4  \\Omega_{\\mathrm{M}0}^2 D(\\chi_1) D(\\chi_2) 
    }{
        a(\\chi_1)  a(\\chi_2) s_1  s_2
    }
    (\\chi_1 - s_1) (5 s_{\\mathrm{b}, 1} - 2) \\times\\\\
    &\\qquad\\qquad\\qquad\\qquad\\qquad\\qquad
    \\left[
        \\ \\mathcal{H}(\\chi_2)  \\mathcal{R}_2 s_1 (f(\\chi_2) - 1) - 5 s_{\\mathrm{b}, 1} + 2
    \\right] \\nonumber
    \\, , \\\\
    %%%%%%%%%%%%%
    J_{31}^{\\kappa \\int\\!\\phi} &=  y \\Delta\\chi^2
    \\, , \\\\
    %%%%%%%%%%%%%
    J_{22}^{\\kappa \\int\\!\\phi} &= 
    \\frac{1}{2} (y^2 - 1) \\chi_1 \\chi_2 
    \\, .
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

- ``b_1 = b(s_1)``, ``s_{\\mathrm{b}, 1} = s_{\\mathrm{b}}(s_1)``, ``\\mathit{f}_{\\mathrm{evo}}``, ... : 
  galaxy bias, magnification bias (i.e. the slope of the luminosity function at the luminosity threshold), 
  and evolution bias (the first two evaluated in ``s_1``); they are
  all stored in `cosmo`;

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

This function is used inside `ξ_GNC_Lensing_IntegratedGP` with [`trapz`](@ref) from the 
[Trapz](https://github.com/francescoalemanno/Trapz.jl) Julia package.


## Inputs

- `IP1::Point`, `IP2::Point`, `P1::Point`, `P2::Point` or `χ1`, `χ2`, `s1`, `s2`: `Point`/comoving 
  distances where the TPCF has to be calculated; they contain all the 
  data of interest needed for this calculus (comoving distance, growth factor and so on).
  
- `y`: the cosine of the angle between the two points `P1` and `P2` wrt the observer

- `cosmo::Cosmology`: cosmology to be used in this computation; it contains all the splines
  used for the conversion `s` -> `Point`, and all the cosmological parameters ``b``, ...

## Keyword Arguments

- `obs::Union{Bool,Symbol} = :noobsvel` : do you want to consider the observer terms in the computation of the 
  chosen GNC TPCF effect?
  - `:yes` or `true` -> all the observer effects will be considered
  - `:no` or `false` -> no observer term will be taken into account
  - `:noobsvel` -> the observer terms related to the observer velocity (that you can find in the CF concerning IntegratedGP)
    will be neglected, the other ones will be taken into account


See also: [`Point`](@ref), [`Cosmology`](@ref), [`ξ_GNC_multipole`](@ref), 
[`map_ξ_GNC_multipole`](@ref), [`print_map_ξ_GNC_multipole`](@ref),
[`ξ_GNC_Lensing_IntegratedGP`](@ref)
"""
integrand_ξ_GNC_Lensing_IntegratedGP


#=
function func_Δχ_min(s1, s2, y; frac = 1e-4)
     Δs = s(s1, s2, y)
     return frac * Δs
end
=#


##########################################################################################92



function ξ_GNC_Lensing_IntegratedGP(P1::Point, P2::Point, y, cosmo::Cosmology;
     en::Float64=1e6, N_χs_2::Int=100, obs::Union{Bool,Symbol}=:noobsvel)

     χ1s = P1.comdist .* range(1e-6, 1, length=N_χs_2)
     χ2s = P2.comdist .* range(1e-6, 1, length=N_χs_2)

     IP1s = [GaPSE.Point(x, cosmo) for x in χ1s]
     IP2s = [GaPSE.Point(x, cosmo) for x in χ2s]

     int_ξs = [
          en * GaPSE.integrand_ξ_GNC_Lensing_IntegratedGP(IP1, IP2, P1, P2, y, cosmo; obs=obs)
          for IP1 in IP1s, IP2 in IP2s
     ]

     res = trapz((χ1s, χ2s), int_ξs)
     #println("res = $res")
     return res / en
end



function ξ_GNC_Lensing_IntegratedGP(s1, s2, y, cosmo::Cosmology; kwargs...)
     P1, P2 = Point(s1, cosmo), Point(s2, cosmo)
     return ξ_GNC_Lensing_IntegratedGP(P1, P2, y, cosmo; kwargs...)
end



"""
     ξ_GNC_Lensing_IntegratedGP(
          P1::Point, P2::Point, y, cosmo::Cosmology;
          en::Float64=1e6, N_χs_2::Int=100, 
          obs::Union{Bool,Symbol}=:noobsvel
          ) ::Float64

     ξ_GNC_Lensing_IntegratedGP(
          s1, s2, y, cosmo::Cosmology; 
          kwargs...) ::Float64

Return the Two-Point Correlation Function (TPCF) given 
by the cross correlation between the Lensing
and the Integrated Gravitational Potential (GP) effects arising 
from the Galaxy Number Counts (GNC).

In the first method, you should pass the two `Point` (`P1` and `P2`) where to 
evaluate the function, while in the second method (that internally recalls the first) 
you must provide the two corresponding comoving distances `s1` and `s2`.
We remember that all the distances are measured in ``h_0^{-1}\\mathrm{Mpc}``.

The analytical expression of this integrand is the following:

```math
\\begin{split}
    \\xi^{\\kappa \\int\\!\\phi} ( s_1 , s_2, y ) = 
    \\int_0^{s_1}\\mathrm{d} \\chi_1 \\int_0^{s_2}\\mathrm{d} \\chi_2 \\;
    J_{\\alpha}^{\\kappa \\int\\!\\phi} 
    \\left[ 
        J_{31}^{\\kappa \\int\\!\\phi} I_1^3 ( \\Delta \\chi ) +
        J_{22}^{\\kappa \\int\\!\\phi} I_2^2 ( \\Delta \\chi ) 
     \\right] \\, ,
\\end{split}
```

where

```math
\\begin{split}
     J_{\\alpha}^{\\kappa \\int\\!\\phi} &=
    \\frac{
        9 \\chi_2 \\ \\mathcal{H}_0^4  \\Omega_{\\mathrm{M}0}^2 D(\\chi_1) D(\\chi_2) 
    }{
        a(\\chi_1)  a(\\chi_2) s_1  s_2
    }
    (\\chi_1 - s_1) (5 s_{\\mathrm{b}, 1} - 2) \\times\\\\
    &\\qquad\\qquad\\qquad\\qquad\\qquad\\qquad
    \\left[
        \\ \\mathcal{H}(\\chi_2)  \\mathcal{R}_2 s_1 (f(\\chi_2) - 1) - 5 s_{\\mathrm{b}, 1} + 2
    \\right] \\nonumber
    \\, , \\\\
    %%%%%%%%%%%%%
    J_{31}^{\\kappa \\int\\!\\phi} &=  y \\Delta\\chi^2
    \\, , \\\\
    %%%%%%%%%%%%%
    J_{22}^{\\kappa \\int\\!\\phi} &= 
    \\frac{1}{2} (y^2 - 1) \\chi_1 \\chi_2 
    \\, .
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

- ``b_1 = b(s_1)``, ``s_{\\mathrm{b}, 1} = s_{\\mathrm{b}}(s_1)``, ``\\mathit{f}_{\\mathrm{evo}}``, ... : 
  galaxy bias, magnification bias (i.e. the slope of the luminosity function at the luminosity threshold), 
  and evolution bias (the first two evaluated in ``s_1``); they are
  all stored in `cosmo`;

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

This function is computed from `integrand_ξ_GNC_Lensing_IntegratedGP` with [`trapz`](@ref) from the 
[Trapz](https://github.com/francescoalemanno/Trapz.jl) Julia package.


## Inputs

-  `P1::Point`, `P2::Point` or `s1`,`s2`: `Point`/comoving 
  distances where the TPCF has to be calculated; they contain all the 
  data of interest needed for this calculus (comoving distance, growth factor and so on).
  
- `y`: the cosine of the angle between the two points `P1` and `P2` wrt the observer

- `cosmo::Cosmology`: cosmology to be used in this computation; it contains all the splines
  used for the conversion `s` -> `Point`, and all the cosmological parameters ``b``, ...

## Keyword Arguments

- `obs::Union{Bool,Symbol} = :noobsvel` : do you want to consider the observer terms in the computation of the 
  chosen GNC TPCF effect?
  - `:yes` or `true` -> all the observer effects will be considered
  - `:no` or `false` -> no observer term will be taken into account
  - `:noobsvel` -> the observer terms related to the observer velocity (that you can find in the CF concerning IntegratedGP)
    will be neglected, the other ones will be taken into account

- `en::Float64 = 1e6`: just a float number used in order to deal better 
  with small numbers;

- `N_χs_2::Int = 100`: number of points to be used for sampling the integral
  along the ranges `(0, s1)` (for `χ1`) and `(0, s2)` (for `χ2`); it has been checked that
  with `N_χs_2 ≥ 50` the result is stable.


See also: [`Point`](@ref), [`Cosmology`](@ref), [`ξ_GNC_multipole`](@ref), 
[`map_ξ_GNC_multipole`](@ref), [`print_map_ξ_GNC_multipole`](@ref),
[`integrand_ξ_GNC_Lensing_IntegratedGP`](@ref)
"""
ξ_GNC_Lensing_IntegratedGP


##########################################################################################92

##########################################################################################92

##########################################################################################92


"""
     ξ_GNC_IntegratedGP_Lensing(s1, s2, y, cosmo::Cosmology; kwargs...) = 
          ξ_GNC_Lensing_IntegratedGP(s2, s1, y, cosmo; kwargs...)

Return the Two-Point Correlation Function (TPCF) given by the cross correlation between the 
Integrated Gravitational Potential (GP) and the Lensing effects arising from the Galaxy Number Counts (GNC).

It's computed through the symmetric function `ξ_GNC_Lensing_IntegratedGP`; check its documentation for
more details about the analytical expression and the keyword arguments.
We remember that all the distances are measured in ``h_0^{-1}\\mathrm{Mpc}``.


## Inputs

- `s1` and `s2`: comoving distances where the TPCF has to be calculated;
  
- `y`: the cosine of the angle between the two points `P1` and `P2` wrt the observer

- `cosmo::Cosmology`: cosmology to be used in this computation; it contains all the splines
  used for the conversion `s` -> `Point`, and all the cosmological parameters ``b``, ...

## Keyword Arguments

- `kwargs...` : Keyword arguments to be passed to the symmetric TPCF

See also: [`Point`](@ref), [`Cosmology`](@ref), [`ξ_GNC_multipole`](@ref), 
[`map_ξ_GNC_multipole`](@ref), [`print_map_ξ_GNC_multipole`](@ref),
[`ξ_GNC_Lensing_IntegratedGP`](@ref)
"""
function ξ_GNC_IntegratedGP_Lensing(s1, s2, y, cosmo::Cosmology; kwargs...)
     ξ_GNC_Lensing_IntegratedGP(s2, s1, y, cosmo; kwargs...)
end

