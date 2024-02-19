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


function integrand_ξ_GNCxLD_IntegratedGP_Lensing(
	IP1::Point, IP2::Point, P1::Point, P2::Point, y, cosmo::Cosmology; 
    b1=nothing, b2=nothing, s_b1=nothing, s_b2=nothing,
    𝑓_evo1=nothing, 𝑓_evo2=nothing, s_lim=nothing)

	s1 = P1.comdist
	s2 = P2.comdist
	χ1, D1, a1, f1, ℋ1 = IP1.comdist, IP1.D, IP1.a, IP1.f, IP1.ℋ
	χ2, D2, a2 = IP2.comdist, IP2.D, IP2.a

	Ω_M0 = cosmo.params.Ω_M0
    s_b_s1 = isnothing(s_b1) ? cosmo.params.s_b1 : s_b1
    𝑓_evo_s1 = isnothing(𝑓_evo1) ? cosmo.params.𝑓_evo1 : 𝑓_evo1

    s_lim = isnothing(s_lim) ? cosmo.params.s_lim : s_lim
    ℛ_s1 = func_ℛ_GNC(s1, P1.ℋ, P1.ℋ_p; s_b=s_b_s1, 𝑓_evo=𝑓_evo_s1, s_lim=s_lim)

	Δχ_square = χ1^2 + χ2^2 - 2 * χ1 * χ2 * y
	Δχ = √(Δχ_square) > 1e-8 ? √(Δχ_square) : 1e-8

	prefactor = 9 / 2 * ℋ0^4 * Ω_M0^2
	factor =  - D1 * D2 * χ1 * (χ2 - s2) / (s1 * s2 * a1 * a2)
	parenth = s1 * ℋ1 * ℛ_s1 * (f1 - 1) - 5 * s_b_s1 + 2 

	new_J31 = 2 * y * Δχ^2
	new_J22 = χ1 * χ2 * (y^2 - 1)

	I13 = cosmo.tools.I13(Δχ)
	I22 = cosmo.tools.I22(Δχ)

	res = prefactor * factor * parenth * (new_J22 * I22 + new_J31 * I13)

	return res
end


function integrand_ξ_GNCxLD_IntegratedGP_Lensing(
	χ1::Float64, χ2::Float64, s1::Float64, s2::Float64, y, cosmo::Cosmology;
	kwargs...)

	P1, P2 = Point(s1, cosmo), Point(s2, cosmo)
	IP1, IP2 = Point(χ1, cosmo), Point(χ2, cosmo)
	return integrand_ξ_GNCxLD_IntegratedGP_Lensing(IP1, IP2, P1, P2, y, cosmo; kwargs...)
end




"""
	integrand_ξ_GNCxLD_IntegratedGP_Lensing(
		IP1::Point, IP2::Point, P1::Point, P2::Point, y, cosmo::Cosmology; 
		b1=nothing, b2=nothing, s_b1=nothing, s_b2=nothing,
		𝑓_evo1=nothing, 𝑓_evo2=nothing, s_lim=nothing ) ::Float64

	integrand_ξ_GNCxLD_IntegratedGP_Lensing(
		χ1::Float64, χ2::Float64, s1::Float64, s2::Float64, 
		y, cosmo::Cosmology; kwargs... ) ::Float64

Return the integrand of the Two-Point Correlation Function (TPCF) given by the cross correlation 
between the Integrated Gravitational Potential (GP) effect arising from the 
Galaxy Number Counts (GNC) and the Lensing
one arising from the Luminosity Distance (LD) perturbations.

In the first method, you should pass the two extreme `Point`s (`P1` and `P2`) and the 
two intermediate integrand `Point`s (`IP1` and `IP2`) where to 
evaluate the function. In the second method (that internally recalls the first),
you must provide the three corresponding comoving distances `s1`, `s2`, `χ1`, `χ2`.
We remember that all the distances are measured in ``h_0^{-1}\\mathrm{Mpc}``.

The analytical expression of this integrand is the following:

```math
\\begin{split}
    f^{\\int\\!\\phi \\, \\kappa} (\\chi_1, \\chi_2, s_1 , s_2, y ) = 
    \\mathfrak{J}_{\\alpha}^{\\int\\!\\phi \\, \\kappa} 
    \\left[ 
        \\mathfrak{J}_{31}^{\\int\\!\\phi \\, \\kappa} I_1^3 ( \\Delta \\chi ) +
        \\mathfrak{J}_{22}^{\\int\\!\\phi \\, \\kappa} I_2^2 ( \\Delta \\chi ) 
     \\right] \\,
\\end{split}
```

with

```math
\\begin{split}
    \\mathfrak{J}_{\\alpha}^{\\int\\!\\phi \\, \\kappa} &=
    - \\frac{
        9 \\chi_1 \\mathcal{H}_0^4 \\Omega_{\\mathrm{M}0}^2 D(\\chi_1) D(\\chi_2) 
    }{
        a(\\chi_1)  a(\\chi_2) s_1  s_2
    }
    (\\chi_2 - s_2)
    \\left[
        \\mathcal{H}(\\chi_1) \\mathcal{R}_1 s_1 (f(\\chi_1) - 1) - 5 s_{\\mathrm{b}, 1}  + 2
    \\right]
    \\, , \\\\
    %%%%%%%%%%%%%
    \\mathfrak{J}_{31}^{\\int\\!\\phi \\, \\kappa} &=  y \\Delta\\chi^2
    \\, , \\\\
    %%%%%%%%%%%%%
    \\mathfrak{J}_{22}^{\\int\\!\\phi \\, \\kappa} &= 
    \\frac{1}{2} (y^2 - 1) \\chi_1 \\chi_2 
    \\,  ,
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

- ``\\mathfrak{R}_1 = \\mathfrak{R}(s_1)``, ... is 
  computed by `func_ℛ_LD` in `cosmo::Cosmology` (and evaluated in ``s_1`` );
  the definition of ``\\mathcal{R}(s)`` is the following:
  ```math
  \\mathfrak{R}(s) = 1 - \\frac{1}{\\mathcal{H}(s) s} ;
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


This function is used inside `ξ_GNCxLD_IntegratedGP_Lensing` with the trapz() from the 
[Trapz](https://github.com/francescoalemanno/Trapz.jl) Julia package.

## Inputs

-  `IP1::Point`, `IP2::Point`, `P1::Point`, `P2::Point` or `χ1`,`χ2`,`s1`,`s2`: `Point`/comoving 
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
  If `nothing`, the fault value stored in `cosmo` will be considered.


See also: [`Point`](@ref), [`Cosmology`](@ref), [`ξ_GNCxLD_multipole`](@ref), 
[`map_ξ_GNCxLD_multipole`](@ref), [`print_map_ξ_GNCxLD_multipole`](@ref)
"""
integrand_ξ_GNCxLD_IntegratedGP_Lensing



##########################################################################################92



function ξ_GNCxLD_IntegratedGP_Lensing(P1::Point, P2::Point, y, cosmo::Cosmology;
	en::Float64 = 1e6, N_χs_2::Int = 100, kwargs...)

	χ1s = P1.comdist .* range(1e-6, 1.0, length = N_χs_2)
	χ2s = P2.comdist .* range(1e-6, 1.0, length = N_χs_2 + 7)

	IP1s = [GaPSE.Point(x, cosmo) for x in χ1s]
	IP2s = [GaPSE.Point(x, cosmo) for x in χ2s]

	int_ξs = [
		en * GaPSE.integrand_ξ_GNCxLD_IntegratedGP_Lensing(IP1, IP2, P1, P2, y, cosmo; kwargs...)
		for IP1 in IP1s, IP2 in IP2s
	]

	res = trapz((χ1s, χ2s), int_ξs)
	#println("res = $res")
	return res / en
end


function ξ_GNCxLD_IntegratedGP_Lensing(s1, s2, y, cosmo::Cosmology; kwargs...)
	P1, P2 = Point(s1, cosmo), Point(s2, cosmo)
	return ξ_GNCxLD_IntegratedGP_Lensing(P1, P2, y, cosmo; kwargs...)
end



"""
    ξ_GNCxLD_IntegratedGP_Lensing(
        s1, s2, y, cosmo::Cosmology;
        b1=nothing, b2=nothing, s_b1=nothing, s_b2=nothing,
        𝑓_evo1=nothing, 𝑓_evo2=nothing, s_lim=nothing, 
        en::Float64 = 1e6, N_χs::Int = 100 ) ::Float64

    ξ_GNCxLD_IntegratedGP_Lensing(
        s1, s2, y, cosmo::Cosmology; 
        kwargs... ) ::Float64

Return the Two-Point Correlation Function (TPCF) given by the cross correlation 
between the Integrated Gravitational Potential (GP) effect arising from the 
Galaxy Number Counts (GNC) and the Lensing
one arising from the Luminosity Distance (LD) perturbations.

In the first method, you should pass the two extreme `Point`s (`P1` and `P2`) where to 
evaluate the function. In the second method (that internally recalls the first),
you must provide the two corresponding comoving distances `s1`, `s2`.
We remember that all the distances are measured in ``h_0^{-1}\\mathrm{Mpc}``.

The analytical expression of this TPCF is the following:

```math
\\begin{split}
    \\xi^{\\int\\!\\phi \\, \\kappa} ( s_1 , s_2, y ) = 
    \\int_0^{s_1}\\! \\mathrm{d}\\chi_1 \\int_0^{s_2}\\! \\mathrm{d}\\chi_2 \\;
    \\mathfrak{J}_{\\alpha}^{\\int\\!\\phi \\, \\kappa} 
    \\left[ 
        \\mathfrak{J}_{31}^{\\int\\!\\phi \\, \\kappa} I_1^3 ( \\Delta \\chi ) +
        \\mathfrak{J}_{22}^{\\int\\!\\phi \\, \\kappa} I_2^2 ( \\Delta \\chi ) 
     \\right] \\,
\\end{split}
```

with

```math
\\begin{split}
    \\mathfrak{J}_{\\alpha}^{\\int\\!\\phi \\, \\kappa} &=
    - \\frac{
        9 \\chi_1 \\mathcal{H}_0^4 \\Omega_{\\mathrm{M}0}^2 D(\\chi_1) D(\\chi_2) 
    }{
        a(\\chi_1)  a(\\chi_2) s_1  s_2
    }
    (\\chi_2 - s_2)
    \\left[
        \\mathcal{H}(\\chi_1) \\mathcal{R}_1 s_1 (f(\\chi_1) - 1) - 5 s_{\\mathrm{b}, 1}  + 2
    \\right]
    \\, , \\\\
    %%%%%%%%%%%%%
    \\mathfrak{J}_{31}^{\\int\\!\\phi \\, \\kappa} &=  y \\Delta\\chi^2
    \\, , \\\\
    %%%%%%%%%%%%%
    \\mathfrak{J}_{22}^{\\int\\!\\phi \\, \\kappa} &= 
    \\frac{1}{2} (y^2 - 1) \\chi_1 \\chi_2 
    \\,  ,
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

- ``\\mathfrak{R}_1 = \\mathfrak{R}(s_1)``, ... is 
  computed by `func_ℛ_LD` in `cosmo::Cosmology` (and evaluated in ``s_1`` );
  the definition of ``\\mathcal{R}(s)`` is the following:
  ```math
  \\mathfrak{R}(s) = 1 - \\frac{1}{\\mathcal{H}(s) s} ;
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



The computation is made applying trapz() (see the 
[Trapz](https://github.com/francescoalemanno/Trapz.jl) Julia package) to
the integrand function `integrand_ξ_GNCxLD_IntegratedGP_Lensing`.

## Inputs

- `P1::Point`, `P2::Point` or `s1`,`s2`: `Point`/comoving 
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
  If `nothing`, the fault value stored in `cosmo` will be considered.

- `en::Float64 = 1e6`: just a float number used in order to deal better 
  with small numbers;

- `N_χs_2::Int = 100`: number of points to be used for sampling the integral
  along the ranges `(0, s1)` (for `χ1`) and `(0, s2)` (for `χ2`); it has been checked that
  with `N_χs_2 ≥ 50` the result is stable.

See also: [`Point`](@ref), [`Cosmology`](@ref), [`ξ_GNCxLD_multipole`](@ref), 
[`map_ξ_GNCxLD_multipole`](@ref), [`print_map_ξ_GNCxLD_multipole`](@ref)
"""
ξ_GNCxLD_IntegratedGP_Lensing



##########################################################################################92

##########################################################################################92

##########################################################################################92




"""
    ξ_LDxGNC_Lensing_IntegratedGP(s1, s2, y, cosmo::Cosmology;         
        b1=nothing, b2=nothing, s_b1=nothing, s_b2=nothing, 
        𝑓_evo1=nothing, 𝑓_evo2=nothing, s_lim=nothing,
        obs::Union{Bool, Symbol} = :noobsvel ) ::Float64

Return the Two-Point Correlation Function (TPCF) given by the cross correlation between the Lensing 
effect arising from the Luminosity Distance (LD) perturbations and Integrated
Gravitational Potential (GP) one arising from the Galaxy Number Counts (GNC).

It's computed through the symmetric function `ξ_GNCxLD_IntegratedGP_Lensing`; check its documentation for
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
[`map_ξ_LDxGNC_multipole`](@ref), [`print_map_ξ_LDxGNC_multipole`](@ref),
[`ξ_GNCxLD_IntegratedGP_Lensing`](@ref)
"""
function ξ_LDxGNC_Lensing_IntegratedGP(s1, s2, y, cosmo::Cosmology; 
    b1=nothing, b2=nothing, s_b1=nothing, s_b2=nothing,
    𝑓_evo1=nothing, 𝑓_evo2=nothing, s_lim=nothing, kwargs...)

    b1 = isnothing(b1) ? cosmo.params.b1 : b1
    b2 = isnothing(b2) ? cosmo.params.b2 : b2
    s_b1 = isnothing(s_b1) ? cosmo.params.s_b1 : s_b1
    s_b2 = isnothing(s_b2) ? cosmo.params.s_b2 : s_b2
    𝑓_evo1 = isnothing(𝑓_evo1) ? cosmo.params.𝑓_evo1 : 𝑓_evo1
    𝑓_evo2 = isnothing(𝑓_evo2) ? cosmo.params.𝑓_evo2 : 𝑓_evo2

    ξ_GNCxLD_IntegratedGP_Lensing(s2, s1, y, cosmo; 
        b1=b2, b2=b1, s_b1=s_b2, s_b2=s_b1,
        𝑓_evo1=𝑓_evo2, 𝑓_evo2=𝑓_evo1, s_lim=s_lim, kwargs...)
end
