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


function integrand_ξ_GNC_IntegratedGP(IP1::Point, IP2::Point,
  P1::Point, P2::Point,
  y, cosmo::Cosmology; obs::Union{Bool,Symbol}=:noobsvel)

  s1, ℛ_s1 = P1.comdist, P1.ℛ_GNC
  s2, ℛ_s2 = P2.comdist, P2.ℛ_GNC
  χ1, D1, a1, ℋ1, f1 = IP1.comdist, IP1.D, IP1.a, IP1.ℋ, IP1.f
  χ2, D2, a2, ℋ2, f2 = IP2.comdist, IP2.D, IP2.a, IP2.ℋ, IP2.f
  s_b_s1, s_b_s2 = cosmo.params.s_b, cosmo.params.s_b
  Ω_M0 = cosmo.params.Ω_M0

  Δχ_square = χ1^2 + χ2^2 - 2 * χ1 * χ2 * y
  Δχ = Δχ_square > 0 ? √(Δχ_square) : 0

  factor = 9 * Δχ^4 * ℋ0^4 * Ω_M0^2 * D1 * D2 / (s1 * s2 * a1 * a2)
  parenth_1 = s1 * ℋ1 * ℛ_s1 * (f1 - 1) - 5 * s_b_s1 + 2
  parenth_2 = s2 * ℋ2 * ℛ_s2 * (f2 - 1) - 5 * s_b_s2 + 2

  I04_tilde = cosmo.tools.I04_tilde(Δχ)

  return factor * parenth_1 * parenth_2 * I04_tilde
end

function integrand_ξ_GNC_IntegratedGP(
  χ1::Float64, χ2::Float64,
  s1::Float64, s2::Float64,
  y, cosmo::Cosmology;
  kwargs...)

  P1, P2 = Point(s1, cosmo), Point(s2, cosmo)
  IP1, IP2 = Point(χ1, cosmo), Point(χ2, cosmo)
  return integrand_ξ_GNC_IntegratedGP(IP1, IP2, P1, P2, y, cosmo; kwargs...)
end



"""
     integrand_ξ_GNC_IntegratedGP(
          IP1::Point, IP2::Point,
          P1::Point, P2::Point,
          y, cosmo::Cosmology; 
          obs::Union{Bool,Symbol}=:noobsvel
          ) ::Float64

     integrand_ξ_GNC_IntegratedGP(
          χ1::Float64, χ2::Float64,
          s1::Float64, s2::Float64,
          y, cosmo::Cosmology;
          kwargs...) ::Float64

Return the integrand of the Two-Point Correlation Function (TPCF) of the 
Integrated Gravitational Potential (GP) auto-correlation effect 
arising from the Galaxy Number Counts (GNC).

In the first method, you should pass the two extreme `Point`s (`P1` and `P2`) and the two 
intermediate integrand `Point`s (`IP1` and `IP2`) where to 
evaluate the function. In the second method (that internally recalls the first),
you must provide the four corresponding comoving distances `s1`, `s2`, `χ1`, `χ2`.
We remember that all the distances are measured in ``h_0^{-1}\\mathrm{Mpc}``.

The analytical expression of this integrand is the following:

```math
\\begin{equation}
    f^{\\int \\!\\phi\\int \\!\\phi} (\\chi_1, \\chi_2, s_1, s_2, y) = 
    J^{\\int \\!\\phi \\int \\!\\phi}_{40} 
    \\tilde{I}_0^4 ( \\Delta\\chi) \\, , 
\\end{equation}
```

where

```math
\\begin{split}
    J^{\\int \\!\\phi\\int \\!\\phi}_{40} =
    \\frac{
        9 \\Delta\\chi ^4 \\mathcal{H}_0^4 \\Omega_{\\mathrm{M}0}^2 D(\\chi_1) D(\\chi_2)
    }{
        a(\\chi_1) a(\\chi_2) s_1 s_2
    }
    &\\left[
        s_1 (f(\\chi_1) - 1) \\mathcal{H}(\\chi_1) \\mathcal{R}_1 - 5 s_{\\mathrm{b}, 1} + 2
    \\right] \\times
    \\nonumber \\\\
    &\\left[
        s_2 (f(\\chi_2) - 1) \\mathcal{H}(\\chi_2) \\mathcal{R}_2 - 5 s_{\\mathrm{b}, 2} + 2
    \\right]
    \\, .
\\end{split}
```

where:

- ``s_1`` and ``s_2`` are comoving distances;

- ``D_1 = D(s_1)``, ... is the linear growth factor (evaluated in ``s_1``);

- ``a_1 = a(s_1)``, ... is the scale factor (evaluated in ``s_1``);

- ``f_1 = f(s_1)``, ... is the linear growth rate (evaluated in ``s_1``);

- ``\\mathcal{H}_1 = \\mathcal{H}(s_1)``, ... is the comoving 
  Hubble parameter (evaluated in ``s_1``, ...);

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

This function is used inside `ξ_GNC_IntegratedGP` with [`trapz`](@ref) from the 
[Trapz](https://github.com/francescoalemanno/Trapz.jl) Julia package.


## Inputs

-  `IP1::Point`, `IP2::Point`, `P1::Point`, `P2::Point` or `χ1`, `χ2`, `s1`, `s2`: `Point`/comoving 
  distances where the TPCF has to be calculated; they contain all the 
  data of interest needed for this calculus (comoving distance, growth factor and so on).
  
- `y`: the cosine of the angle between the two points `P1` and `P2` wrt the observer

- `cosmo::Cosmology`: cosmology to be used in this computation; it contains all the splines
  used for the conversion `s` -> `Point`, and all the cosmological parameters ``b``, ...

## Keyword Arguments

- `obs::Union{Bool,Symbol} = :noobsvel` : do you want to consider the observer terms in the computation of the 
  chosen GNC TPCF effect?
  - `:yes` or `true` -> all the observer terms will be considered;
  - `:no` or `false` -> no observer term will be taken into account;
  - `:noobsvel` -> the observer terms related to the observer velocity (that you can find in the CF concerning Doppler)
    will be neglected, the other ones will be taken into account.


See also: [`Point`](@ref), [`Cosmology`](@ref), [`ξ_GNC_multipole`](@ref), 
[`map_ξ_GNC_multipole`](@ref), [`print_map_ξ_GNC_multipole`](@ref),
[`ξ_GNC_IntegratedGP`](@ref)
"""
integrand_ξ_GNC_IntegratedGP



##########################################################################################92




function ξ_GNC_IntegratedGP(P1::Point, P2::Point, y, cosmo::Cosmology;
  en::Float64=1e10, N_χs_2::Int=100, obs::Union{Bool,Symbol}=:noobsvel)

  #adim_χs = range(1e-12, 1, N_χs)
  #Δχ_min = func_Δχ_min(s1, s2, y; frac = frac_Δχ_min)

  χ1s = P1.comdist .* range(1e-6, 1, length=N_χs_2)
  χ2s = P2.comdist .* range(1e-6, 1, length=N_χs_2)

  IP1s = [GaPSE.Point(x, cosmo) for x in χ1s]
  IP2s = [GaPSE.Point(x, cosmo) for x in χ2s]

  int_ξ_igp = [
    en * GaPSE.integrand_ξ_GNC_IntegratedGP(IP1, IP2, P1, P2, y, cosmo; obs=obs)
    for IP1 in IP1s, IP2 in IP2s
  ]

  res = trapz((χ1s, χ2s), int_ξ_igp)
  #println("res = $res")

  #=
  χ1s = [x for x in range(0, P1.comdist, length = N_χs)[begin+1:end]]
  l = Int(floor(N_χs/2))
  matrix_χ2s = [begin
       a = [x for x in range(0, P2.comdist, length=l)[begin+1:end]];
       b = [x for x in range(x1-focus, x1+focus, length=l)];
       vcat(a[a.<x1-focus], b, a[a.>x1+focus])
       end for x1 in χ1s]

  IP1s = [GaPSE.Point(x, cosmo) for x in χ1s]
  matrix_IP2s = [[GaPSE.Point(x, cosmo) for x in y] for y in matrix_χ2s]
  matrix_int_ξs = [
       [en * GaPSE.integrand_ξ_GNC_IntegratedGP(IP1, IP2, P1, P2, y, cosmo) 
       for IP2 in matrix_IP2s[i]]
       for (i,IP1) in enumerate(IP1s)]

  vec_trapz = [trapz(χ2s,int_ξs) for (χ2s,int_ξs) in zip(matrix_χ2s, matrix_int_ξs)]
  res = trapz(χ1s, vec_trapz)
  =#
  return res / en
end


function ξ_GNC_IntegratedGP(s1, s2, y, cosmo::Cosmology; kwargs...)
  P1, P2 = Point(s1, cosmo), Point(s2, cosmo)
  return ξ_GNC_IntegratedGP(P1, P2, y, cosmo; kwargs...)
end



"""
     ξ_GNC_IntegratedGP(
          P1::Point, P2::Point, y, cosmo::Cosmology;
          en::Float64=1e10, N_χs_2::Int=100, 
          obs::Union{Bool,Symbol}=:noobsvel
          ) ::Float64

     ξ_GNC_IntegratedGP(s1, s2, y, cosmo::Cosmology; 
          kwargs...) ::Float64

Returns the Two-Point Correlation Function (TPCF) of the Integrated
Gravitational Potential (GP) auto-correlation effect
arising from the Galaxy Number Counts (GNC).

In the first method, you should pass the two `Point` (`P1` and `P2`) where to 
evaluate the function, while in the second method (that internally recalls the first) 
you must provide the two corresponding comoving distances `s1` and `s2`.
We remember that all the distances are measured in ``h_0^{-1}\\mathrm{Mpc}``.

The analytical expression of this term is the following:

```math
\\begin{split}
    \\xi^{\\int\\!\\phi \\int \\!\\phi }( s_1 , s_2, y ) = 
    \\int_0^{s_1}\\mathrm{d} \\chi_1  \\int_0^{s_2}\\mathrm{d} \\chi_2 \\;  
    J^{\\int \\!\\phi \\int \\!\\phi}_{40} 
    \\tilde{I}_0^4 ( \\Delta\\chi) \\, , 
\\end{split}
```

where

```math
\\begin{split}
    J^{\\int \\!\\phi\\int \\!\\phi}_{40} =
    \\frac{
        9 \\Delta\\chi ^4 \\mathcal{H}_0^4 \\Omega_{\\mathrm{M},  0}^2 D(\\chi_1) D(\\chi_2)
    }{
        a(\\chi_1) a(\\chi_2) s_1 s_2
    }
    &\\left[
        s_1 (f(\\chi_1) - 1) \\mathcal{H}(\\chi_1) \\mathcal{R}_1 - 5 s_{\\mathrm{b},  1} + 2
    \\right] \\times
    \\nonumber \\\\
    &\\left[
        s_2 (f(\\chi_2) - 1) \\mathcal{H}(\\chi_2) \\mathcal{R}_2 - 5 s_{\\mathrm{b},  2} + 2
    \\right]
    \\, .
\\end{split}
```

where:

- ``s_1`` and ``s_2`` are comoving distances;

- ``D_1 = D(s_1)``, ... is the linear growth factor (evaluated in ``s_1``);

- ``a_1 = a(s_1)``, ... is the scale factor (evaluated in ``s_1``);

- ``f_1 = f(s_1)``, ... is the linear growth rate (evaluated in ``s_1``);

- ``\\mathcal{H}_1 = \\mathcal{H}(s_1)``, ... is the comoving 
  Hubble parameter (evaluated in ``s_1``, ...);

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

- ``\\mathcal{H}_0``, ``f_0`` and so on are evaluated at the observer position (i.e. at present day);

- ``\\Delta\\chi_1 := \\sqrt{\\chi_1^2 + s_2^2-2\\,\\chi_1\\,s_2\\,y}`` and 
  ``\\Delta\\chi_2 := \\sqrt{s_1^2 + \\chi_2^2-2\\,s_1\\,\\chi_2\\,y}``;

- ``s=\\sqrt{s_1^2 + s_2^2 - 2 \\, s_1 \\, s_2 \\, y}`` and 
  ``\\Delta\\chi := \\sqrt{\\chi_1^2 + \\chi_2^2-2\\,\\chi_1\\,\\chi_2\\,y}``.

In this TPCF there are no observer terms. The `obs` keyword is inserted only for compatibility with 
the other GNC TPCFs.

This function is computed integrating `integrand_ξ_GNC_IntegratedGP` with [`trapz`](@ref) from the 
[Trapz](https://github.com/francescoalemanno/Trapz.jl) Julia package.


## Inputs

- `P1::Point` and `P2::Point`, or `s1` and `s2`: `Point`/comoving distances where the 
  TPCF has to be calculated; they contain all the 
  data of interest needed for this calculus (comoving distance, growth factor and so on);
  
- `y`: the cosine of the angle between the two points `P1` and `P2` wrt the observer;

- `cosmo::Cosmology`: cosmology to be used in this computation; it contains all the splines
  used for the conversion `s` -> `Point`, and all the cosmological parameters ``b``, ...

## Keyword Arguments

- `obs::Union{Bool,Symbol} = :noobsvel` : do you want to consider the observer terms in the computation of the 
  chosen GNC TPCF effect?
  - `:yes` or `true` -> all the observer terms will be considered;
  - `:no` or `false` -> no observer term will be taken into account;
  - `:noobsvel` -> the observer terms related to the observer velocity (that you can find in the CF concerning Doppler)
    will be neglected, the other ones will be taken into account.

- `en::Float64 = 1e6`: just a float number used in order to deal better 
  with small numbers;

- `N_χs_2::Int = 100`: number of points to be used for sampling the integral
  along the ranges `(0, s1)` (for `χ1`) and `(0, s2)` (for `χ2`); it has been checked that
  with `N_χs_2 ≥ 50` the result is stable.

See also: [`Point`](@ref), [`Cosmology`](@ref), [`ξ_GNC_multipole`](@ref), 
[`map_ξ_GNC_multipole`](@ref), [`print_map_ξ_GNC_multipole`](@ref),
[`integrand_ξ_GNC_IntegratedGP`](@ref)
"""
ξ_GNC_IntegratedGP
