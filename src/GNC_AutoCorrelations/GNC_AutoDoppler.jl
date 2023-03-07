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


function ξ_GNC_Doppler(P1::Point, P2::Point, y, cosmo::Cosmology; obs::Union{Bool,Symbol}=:noobsvel)

     s1, D1, f1, ℋ1, ℛ1 = P1.comdist, P1.D, P1.f, P1.ℋ, P1.ℛ_GNC
     s2, D2, f2, ℋ2, ℛ2 = P2.comdist, P2.D, P2.f, P2.ℋ, P2.ℛ_GNC
     s_b1, s_b2 = cosmo.params.s_b, cosmo.params.s_b

     Δs = s(P1.comdist, P2.comdist, y)
     common = D1 * D2 * f1 * f2 * ℛ1 * ℛ2 * ℋ1 * ℋ2
     factor = 3 * s1 * s2 - 2 * y * (s1^2 + s2^2) + s1 * s2 * y^2

     I00 = cosmo.tools.I00(Δs)
     I20 = cosmo.tools.I20(Δs)
     I40 = cosmo.tools.I40(Δs)
     I02 = cosmo.tools.I02(Δs)

     parenth = 1 / 45 * I00 + 2 / 63 * I20 + 1 / 105 * I40


     if obs == false || obs == :no || obs == :noobsvel
          return common * (factor * parenth + 1 / 3 * y * Δs^2 * I02)
     elseif obs == true || obs == :yes

          #### New observer terms #########

          I13_s1, I13_s2 = cosmo.tools.I13(s1), cosmo.tools.I13(s2)
          I31_s1, I31_s2 = cosmo.tools.I31(s1), cosmo.tools.I31(s2)
          I11_s1, I11_s2 = cosmo.tools.I11(s1), cosmo.tools.I11(s2)
          σ2 = cosmo.tools.σ_2

          obs_common_12 = y * f0 * ℋ0 * s1^2 * f1 * ℋ1 * ℛ1 * (ℛ2 - 5 * s_b2 + 2)
          obs_common_21 = y * f0 * ℋ0 * s2^2 * f2 * ℋ2 * ℛ2 * (ℛ1 - 5 * s_b1 + 2)
          J_σ2 = 1 / 3 * y * f0^2 * ℋ0^2 * (ℛ1 - 5 * s_b1 + 2) * (ℛ2 - 5 * s_b2 + 2)

          obs_terms_12 = D1 * obs_common_12 * (1 / 5 * (I11_s1 + I31_s1) - I13_s1)
          obs_terms_21 = D2 * obs_common_21 * (1 / 5 * (I11_s2 + I31_s2) - I13_s2)

          obs_terms = obs_terms_12 + obs_terms_21 + J_σ2 * σ2

          #################################

          return common * (factor * parenth + 1 / 3 * y * Δs^2 * I02) + obs_terms
     else
          throw(AssertionError(":$obs is not a valid Symbol for \"obs\"; they are: \n\t" *
                               "$(":".*string.(VALID_OBS_VALUES) .* vcat([" , " for i in 1:length(VALID_OBS_VALUES)-1], " .")... )"
          ))
     end

end


function ξ_GNC_Doppler(s1, s2, y, cosmo::Cosmology; obs::Union{Bool,Symbol}=:noobsvel)
     P1, P2 = Point(s1, cosmo), Point(s2, cosmo)
     return ξ_GNC_Doppler(P1, P2, y, cosmo; obs=obs)
end



"""
     ξ_GNC_Doppler(P1::Point, P2::Point, y, cosmo::Cosmology; 
        obs::Union{Bool, Symbol} = :noobsvel) ::Float64

     ξ_GNC_Doppler(s1, s2, y, cosmo::Cosmology; 
        kwargs...) ::Float64

Returns the Two-Point Correlation Function (TPCF) of the Doppler auto-correlation effect
arising from the Galaxy Number Counts (GNC).

In the first method, you should pass the two `Point` (`P1` and `P2`) where to 
evaluate the function, while in the second method (that internally recalls the first) 
you must provide the two corresponding comoving distances `s1` and `s2`.
We remember that all the distances are measured in ``h_0^{-1}\\mathrm{Mpc}``.

The analytical expression of this term is the following:

```math
\\begin{split}
    \\xi^{v_{\\parallel} v_{\\parallel}}( s_1, s_2 , y ) &=
    D_1 D_2 J_{\\alpha}^{v_{\\parallel} v_{\\parallel}}
    \\left[  
        J_{\\beta}^{v_{\\parallel} v_{\\parallel}}
        \\left( 
            \\frac{1}{45} I_0^0(s) +
            \\frac{2}{63} I_2^0 (s) + 
            \\frac{1}{105} I_4^0(s)
        \\right) +
        J^{v_{\\parallel} v_{\\parallel}}_{20} I_0^2 (s)
    \\right] 
     \\\\
    & + D_1 \\left[ 
        J^{v_{\\parallel} v_{\\parallel}}_{31} ( s_1 , s_2) I_1^3 (s_1) +   
        J^{v_{\\parallel} v_{\\parallel}}_{11} ( s_1 , s_2) I_1^1 (s_1) + 
        J^{v_{\\parallel} v_{\\parallel}}_{13} ( s_1 ,s_2) I_3^1 (s_1) 
    \\right] 
     \\\\
    & + D_2 \\left[
        J^{v_{\\parallel} v_{\\parallel}}_{31} ( s_2 , s_1) I_1^3 (s_2) +  
        J^{v_{\\parallel} v_{\\parallel}}_{11} ( s_2 , s_1) I_1^1 (s_2) +
        J^{v_{\\parallel} v_{\\parallel}}_{13} ( s_2 , s_1) I_3^1 (s_2) 
    \\right]
     \\\\
    & + J^{v_{\\parallel} v_{\\parallel}}_{\\sigma2} \\sigma_2 \\, ,
\\end{split}
```
with

```math
\\begin{split}
    J_{\\alpha}^{v_{\\parallel} v_{\\parallel}} &= 
    f_1 f_2\\mathcal{H}_1 \\mathcal{H}_2 \\mathcal{R}_1 \\mathcal{R}_2
    \\, , \\\\
    %%%%%%%%%%%%%%
    J_{\\beta}^{v_{\\parallel} v_{\\parallel}} &= 
    y^2 s_1 s_2 - 2 y (s_1^2 + s_2^2) + 3 s_1 s_2
    \\, , \\\\
    %%%%%%%%%%%%%%
    J_{20}^{v_{\\parallel} v_{\\parallel}}  &=
    \\frac{1}{3} y s^2
    \\, , \\\\
    %%%%%%%%%%%%%%
    J_{31}^{v_{\\parallel} v_{\\parallel}} (s_1, s_2)  &=
    -y f_0 \\mathcal{H}_0 s_1^2 f_1 \\mathcal{R}_1 (\\mathcal{R}_2 - 5 s_{\\mathrm{b}, 2} + 2)
    \\, , \\\\
    %%%%%%%%%%%%%%
    J_{11}^{v_{\\parallel} v_{\\parallel}} (s_1, s_2)  &=
    \\frac{1}{5} y f_0 \\mathcal{H}_0 s_1^2 f_1 \\mathcal{H}_1 \\mathcal{R}_1 
    (\\mathcal{R}_2 - 5 s_{\\mathrm{b}, 2} + 2)
    \\, , \\\\
    %%%%%%%%%%%%%%
    J_{13}^{v_{\\parallel} v_{\\parallel}} (s_1, s_2)  &=
    \\frac{1}{5} y f_0 \\mathcal{H}_0 s_1^2 f_1 \\mathcal{H}_1 \\mathcal{R}_1 
    (\\mathcal{R}_2 - 5 s_{\\mathrm{b}, 2} + 2)
    \\, , \\\\
    %%%%%%%%%%%%%%
    J_{\\sigma2}^{v_{\\parallel} v_{\\parallel}} &=
    \\frac{1}{3} y f_0^2 \\mathcal{H}_0^2 (\\mathcal{R}_1 - 5 s_{\\mathrm{b}, 1} + 2) 
    (\\mathcal{R}_2 - 5 s_{\\mathrm{b}, 2} + 2) 
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

In this TPCF, the only non-observer term is the first one (i.e. the bracket proportional to 
``D(s_1) \\, D(s_2)``). The other three ones are all observer terms related to the observer velocity, so
if you set `obs = :no`, `obs = false` or even `obs = :noobsvel` they will not be computed.

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

See also: [`Point`](@ref), [`Cosmology`](@ref), [`ξ_GNC_multipole`](@ref), 
[`map_ξ_GNC_multipole`](@ref), [`print_map_ξ_GNC_multipole`](@ref)
"""
ξ_GNC_Doppler
