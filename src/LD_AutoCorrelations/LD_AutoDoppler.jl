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


function ξ_LD_Doppler(P1::Point, P2::Point, y, cosmo::Cosmology)
    s1, D1, f1, ℋ1, ℛ1 = P1.comdist, P1.D, P1.f, P1.ℋ, P1.ℛ_LD
    s2, D2, f2, ℋ2, ℛ2 = P2.comdist, P2.D, P2.f, P2.ℋ, P2.ℛ_LD

    Δs = s(P1.comdist, P2.comdist, y)
    prefac = D1 * D2 * f1 * f2 * ℛ1 * ℛ2 * ℋ1 * ℋ2
    c1 = 3 * s1 * s2 - 2 * y * (s1^2 + s2^2) + s1 * s2 * y^2

    I00 = cosmo.tools.I00(Δs)
    I20 = cosmo.tools.I20(Δs)
    I40 = cosmo.tools.I40(Δs)
    I02 = cosmo.tools.I02(Δs)

    parenth = I00 / 45.0 + I20 / 31.5 + I40 / 105.0

    first = prefac * (c1 * parenth + I02 * y * Δs^2 / 3.0)

    return first
end


function ξ_LD_Doppler(s1, s2, y, cosmo::Cosmology)
    P1, P2 = Point(s1, cosmo), Point(s2, cosmo)
    return ξ_LD_Doppler(P1, P2, y, cosmo)
end



"""
    ξ_LD_Doppler(P1::Point, P2::Point, y, cosmo::Cosmology) ::Float64

    ξ_LD_Doppler(s1, s2, y, cosmo::Cosmology) ::Float64

Return the Two-Point Correlation Function (TPCF) of the Doppler auto-correlation effect
arising from the Luminosity Distance (LD) perturbations.

In the first method, you should pass the two `Point` (`P1` and `P2`) where to 
evaluate the function, while in the second method (that internally recalls the first) 
you must provide the two corresponding comoving distances `s1` and `s2`.
We remember that all the distances are measured in ``h_0^{-1}\\mathrm{Mpc}``.

The analytical expression of this term is the following:

```math
\\begin{split}
    \\xi^{v_{\\parallel}v_{\\parallel}} (s_1, s_2, y) = 
    \\mathcal{J}^{v_{\\parallel}v_{\\parallel}}_{\\alpha}
    \\left[
        \\mathcal{J}^{v_{\\parallel}v_{\\parallel}}_{00} I_0^0(s) + 
        \\mathcal{J}^{v_{\\parallel}v_{\\parallel}}_{02} I_2^0(s) +
        \\mathcal{J}^{v_{\\parallel}v_{\\parallel}}_{04} I_4^0(s) + 
        \\mathcal{J}^{v_{\\parallel}v_{\\parallel}}_{20} I_0^2(s)
    \\right] \\, ,
\\end{split}
```
with

```math
\\begin{split}
    \\mathcal{J}^{v_{\\parallel}v_{\\parallel}}_{\\alpha} & = 
    D_1 D_2 f_1 f_2 \\mathcal{H}_1 \\mathcal{H}_2  \\mathfrak{R}_1  \\mathfrak{R}_2
    \\, , \\\\
    %%%%%%%%%%%%%%%%%%%%%%%%
    \\mathcal{J}^{v_{\\parallel}v_{\\parallel}}_{00} & = 
    \\frac{1}{45} \\left[
        y^2 s_1 s_2 - 2y(s_1^2 + s_2^2) + 3s_1 s_2
    \\right]
    \\, , \\\\
    %%%%%%%%%%%%%%%%%%%%%%%%
    \\mathcal{J}^{v_{\\parallel}v_{\\parallel}}_{02}  & = 
    \\frac{2}{63} \\left[
        y^2 s_1 s_2 - 2y(s_1^2 + s_2^2) + 3s_1 s_2
    \\right]
    \\, , \\\\
    %%%%%%%%%%%%%%%%%%%%%%%%
    \\mathcal{J}^{v_{\\parallel}v_{\\parallel}}_{04} & = 
    \\frac{1}{105} \\left[
        y^2 s_1 s_2 - 2y(s_1^2 + s_2^2) + 3s_1 s_2
    \\right]
    \\, , \\\\
    %%%%%%%%%%%%%%%%%%%%%%%%
    \\mathcal{J}^{v_{\\parallel}v_{\\parallel}}_{20} & = \\frac{1}{3} y s^2 
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

## Inputs

- `P1::Point` and `P2::Point`, or `s1` and `s2`: `Point`/comoving distances where the 
  TPCF has to be calculated; they contain all the 
  data of interest needed for this calculus (comoving distance, growth factor and so on).
  
- `y`: the cosine of the angle between the two points `P1` and `P2` wrt the observer

- `cosmo::Cosmology`: cosmology to be used in this computation; it contains all the splines
  used for the conversion `s` -> `Point`, and all the cosmological parameters ``\\Omega_{\\mathrm{M}0}``, ...


See also: [`Point`](@ref), [`Cosmology`](@ref), [`ξ_LD_multipole`](@ref), 
[`map_ξ_LD_multipole`](@ref), [`print_map_ξ_LD_multipole`](@ref)
"""
ξ_LD_Doppler
