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


function ξ_GNC_Newtonian(P1::Point, P2::Point, y, cosmo::Cosmology; 
	b1=nothing, b2=nothing, s_b1=nothing, s_b2=nothing, 𝑓_evo1=nothing, 𝑓_evo2=nothing,
    s_lim=nothing, obs::Union{Bool,Symbol}=:noobsvel)
	
	s1, D1, f1 = P1.comdist, P1.D, P1.f
	s2, D2, f2 = P2.comdist, P2.D, P2.f
	
	b1 = isnothing(b1) ? cosmo.params.b1 : b1
	b2 = isnothing(b2) ? cosmo.params.b2 : b2

	Δs = s(s1, s2, y)

	J00 = 1 / 15 * (15 * b1 * b2 + 5 * (b1 * f2 + b2 * f1) + (2 * y^2 + 1) * f1 * f2)
	J02 = -1 / (21 * Δs^2) * (
		s1^2 * (14 * b2 * f1 + 7 * b1 * f2 * (3 * y^2 - 1) + (11 * y^2 + 1) * f1 * f2)
		+
		s2^2 * (14 * b1 * f2 + 7 * b2 * f1 * (3 * y^2 - 1) + (11 * y^2 + 1) * f1 * f2)
		-
		4 * y * s1 * s2 * (7 * (b2 * f1 + b1 * f2) + (y^2 + 5) * f1 * f2)
	)
	J04 = f1 * f2 / (35 * Δs^4) * (
		4 * (3 * y^2 - 1) * (s1^4 + s2^4)
		+
		3 * (3 + y^2)^2 * s1^2 * s2^2
		-
		8 * y * s1 * s2 * (s1^2 + s2^2) * (3 + y^2)
	)

	I00 = cosmo.tools.I00(Δs)
	I20 = cosmo.tools.I20(Δs)
	I40 = cosmo.tools.I40(Δs)

	res = D1 * D2 * (J00 * I00 + J02 * I20 + J04 * I40)

	return res
end


function ξ_GNC_Newtonian(s1, s2, y, cosmo::Cosmology; kwargs...)
	P1, P2 = Point(s1, cosmo), Point(s2, cosmo)
	return ξ_GNC_Newtonian(P1, P2, y, cosmo; kwargs...)
end

"""
	ξ_GNC_Newtonian(P1::Point, P2::Point, y, cosmo::Cosmology; 
		b1=nothing, b2=nothing, s_b1=nothing, s_b2=nothing, 
		𝑓_evo1=nothing, 𝑓_evo2=nothing, s_lim=nothing,
		obs::Union{Bool, Symbol} = :noobsvel) ::Float64

	ξ_GNC_Newtonian(s1, s2, y, cosmo::Cosmology; 
		kwargs...) ::Float64

Return the Two-Point Correlation Function (TPCF) of the Newtonian auto-correlation effect
arising from the Galaxy Number Counts (GNC).

In the first method, you should pass the two `Point` (`P1` and `P2`) where to 
evaluate the function, while in the second method (that internally recalls the first) 
you must provide the two corresponding comoving distances `s1` and `s2`.
We remember that all the distances are measured in ``h_0^{-1}\\mathrm{Mpc}``.

The analytical expression of this term is the following:

```math
\\begin{split}
    \\xi^{\\delta\\delta} ( s_1 , s_2 , y )  = D_1 D_2 
    \\left[ 
    J_{00}^{\\delta\\delta} I_0^0 (s) + 
    J_{02}^{\\delta\\delta} I_2^0 (s) + 
    J_{04}^{\\delta\\delta} I_4^0 (s) 
    \\right] \\, , 
\\end{split}
```


where

```math
\\begin{split}
    J_{00}^{\\delta\\delta}  &= \\frac{1}{15}  
    \\left[
        15 b_1 b_2 + 5 b_1 f_2 + 5 b_2 f_1  +
        (2 y^2 + 1) f_1 f_2
    \\right]
    \\, , \\\\
    %%%%%%%%%%%%%%%%%%
    J_{02}^{\\delta\\delta}  &= - \\frac{1}{21 s^2}
    \\left\\{
        s_1^2 \\left[
            14 b_2 f_1 + 
            7 b_1 f_2 (3 y^2 - 1) + 
            (11 y^2 + 1) f_1 f_2 
            \\right]  \\right.  \\\\
        &\\left.\\qquad \\qquad  
        + s_2^2 \\left[
            14 b_1 f_2 + 
            7 b_2 f_1 (3 y^2 - 1) + 
            (11 y^2 + 1) f_1 f_2 
            \\right]  \\right. \\\\
        &\\left. \\qquad \\qquad   
        - 4 y s_1 s_2 \\left[
            7 b_2 f_1 + 7 b_1 f_2 +
            (y^2 + 5) f_1 f_2
            \\right]
    \\right\\}
    \\, , \\\\
    %%%%%%%%%%%%%%%%%%
    J_{04}^{\\delta\\delta}  &= \\frac{f_1 f_2}{35 s^4}
    \\left[
        4 (3 y^2 - 1) (s_1^4 + s_2^4) +
        3 (3 + y^2)^2 s_1^2 s_2^2 - 
        8 y s_1 s_2 (s_1^2 + s_2^2) (3 + y^2)
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

## Inputs

- `P1::Point` and `P2::Point`, or `s1` and `s2`: `Point`/comoving distances where the 
  TPCF has to be calculated; they contain all the 
  data of interest needed for this calculus (comoving distance, growth factor and so on).
  
- `y`: the cosine of the angle between the two points `P1` and `P2` wrt the observer

- `cosmo::Cosmology`: cosmology to be used in this computation; it contains all the splines
  used for the conversion `s` -> `Point`, and all the cosmological parameters ``b``, ...

## Keyword Arguments

- `obs::Union{Bool,Symbol} = :noobsvel` : do you want to consider the observer terms in the computation of the 
  chosen GNC TPCF effect?
  - `:yes` or `true` -> all the observer effects will be considered
  - `:no` or `false` -> no observer term will be taken into account
  - `:noobsvel` -> the observer terms related to the observer velocity (that you can find in the CF concerning Doppler)
    will be neglected, the other ones will be taken into account

See also: [`Point`](@ref), [`Cosmology`](@ref), [`ξ_GNC_multipole`](@ref), 
[`map_ξ_GNC_multipole`](@ref), [`print_map_ξ_GNC_multipole`](@ref)
"""
ξ_GNC_Newtonian
