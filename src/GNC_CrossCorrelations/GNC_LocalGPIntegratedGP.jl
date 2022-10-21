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


"""
     integrand_ξ_GNC_LocalGP_IntegratedGP(
          IP::Point, P1::Point, P2::Point,
          y, cosmo::Cosmology) :: Float64

Return the integrand of the LocalGP-IntegratedGP cross-correlation function 
``\\xi^{\\phi\\int\\phi} (s_1, s_2, \\cos{\\theta})``, i.e. the function 
``f(s_1, s_2, y, \\chi_1, \\chi_2)`` defined as follows:  

```math
f(s_1, s_2, y, \\chi_1, \\chi_2) = 
     \\frac{9 \\mathcal{H}_0^4 \\Omega_{M0}^2 D(s_1) (\\mathcal{R}(s_1) +1)}{2 a(s_1)} 
     \\frac{D(\\chi_2) \\Delta\\chi_2^4}{ a(\\chi_2)}
     \\left(
          \\mathcal{H}(\\chi_2)( f(\\chi_2) - 1) \\mathcal{R}(s_2) - \\frac{1}{s_2}
     \\right) \\tilde{I}^4_0(\\Delta\\chi_2)
```
where ``\\mathcal{H} = a H``, 
``\\Delta\\chi_2 = \\sqrt{s_1^2 + \\chi_2^2 - 2 s_1 \\chi_2 \\cos{\\theta}}``, 
``y = \\cos{\\theta} = \\hat{\\mathbf{s}}_1 \\cdot \\hat{\\mathbf{s}}_2``).

## Inputs

- `IP::Point`: `Point` inside the integration limits, placed 
  at comoving distance `χ1`.

- `P1::Point` and `P2::Point`: extreme `Point` of the integration, placed 
  at comoving distance `s1` and `s2` respectively.

- `y`: the cosine of the angle between the two points `P1` and `P2`

- `cosmo::Cosmology`: cosmology to be used in this computation


See also: [`ξ_GNC_LocalGP_IntegratedGP`](@ref), [`int_on_mu_LocalGP_IntegratedGP`](@ref)
[`integral_on_mu`](@ref), [`ξ_GNC_multipole`](@ref)
"""
function integrand_ξ_GNC_LocalGP_IntegratedGP(
     IP::Point, P1::Point, P2::Point,
     y, cosmo::Cosmology; obs::Bool = true)

     s1, D_s1, f_s1, a_s1, ℋ_s1, ℛ_s1 = P1.comdist, P1.D, P1.f, P1.a, P1.ℋ, P1.ℛ_GNC
     s2, ℛ_s2 = P2.comdist, P2.ℛ_GNC
     χ2, D2, a2, f2, ℋ2 = IP.comdist, IP.D, IP.a, IP.f, IP.ℋ
     s_b_s1, s_b_s2 = cosmo.params.s_b, cosmo.params.s_b
     𝑓_evo_s1 = cosmo.params.𝑓_evo
     Ω_M0 = cosmo.params.Ω_M0

     Δχ2_square = s1^2 + χ2^2 - 2 * s1 * χ2 * y
     Δχ2 = Δχ2_square > 0 ? √(Δχ2_square) : 0

     factor = 3 * ℋ0^2 * Ω_M0 * D2 * (s2 * ℋ2 * ℛ_s2 * (f2 - 1) - 5 * s_b_s2 + 2) / (2 * s2 * a2)
     parenth = 2 * f_s1 * ℋ_s1^2 * a_s1 * (𝑓_evo_s1 - 3) + 3 * ℋ0^2 * Ω_M0 * (f_s1 + ℛ_s1 + 5 * s_b_s1 - 2)
     
     I04_tilde = cosmo.tools.I04_tilde(Δχ2)

     if obs == false
          return Δχ2^4 / a_s1 * D_s1 * factor * parenth * I04_tilde
     else
          #### New observer terms #########

          I04_tilde_χ2 = cosmo.tools.I04_tilde(χ2)

          obs_terms = χ2^4 / s1 * factor * (ℋ0 * s1 * ℛ1 * (2 * f0 - 3 * Ω_M0) + 2 * f0 * (5 * s_b_s1 - 2)) * I04_tilde_χ2
          
          #################################

          return Δχ2^4 / a_s1 * D_s1 * factor * parenth * I04_tilde + obs_terms
     end
end


function integrand_ξ_GNC_LocalGP_IntegratedGP(
     χ2::Float64, s1::Float64, s2::Float64,
     y, cosmo::Cosmology;
     kwargs...)

     P1, P2 = Point(s1, cosmo), Point(s2, cosmo)
     IP = Point(χ2, cosmo)
     return integrand_ξ_GNC_LocalGP_IntegratedGP(IP, P1, P2, y, cosmo; kwargs...)
end


"""
     ξ_GNC_LocalGP_IntegratedGP(s1, s2, y, cosmo::Cosmology;
          en::Float64 = 1e6, N_χs::Int = 100):: Float64

Return the LocalGP-IntegratedGP cross-correlation function 
``\\xi^{v_{\\parallel}\\int \\phi} (s_1, s_2, \\cos{\\theta})`` concerning the perturbed
luminosity distance, defined as follows:
    
```math
\\xi^{v_{\\parallel}\\int \\phi} (s_1, s_2, \\cos{\\theta}) = 
     \\frac{9 \\mathcal{H}_0^4 \\Omega_{M0}^2 D(s_1) (\\mathcal{R}(s_1) +1)}{2 a(s_1)} 
     \\int_0^{s_2} \\mathrm{d}\\chi_2 \\frac{D(\\chi_2) \\Delta\\chi_2^4}{ a(\\chi_2)}
     \\left(
          \\mathcal{H}(\\chi_2)( f(\\chi_2) - 1) \\mathcal{R}(s_2) - \\frac{1}{s_2}
     \\right) \\tilde{I}^4_0(\\Delta\\chi_2)
```
where ``\\mathcal{H} = a H``, 
``\\Delta\\chi_2 = \\sqrt{s_1^2 + \\chi_2^2 - 2 s_1 \\chi_2 \\cos{\\theta}}``, 
``y = \\cos{\\theta} = \\hat{\\mathbf{s}}_1 \\cdot \\hat{\\mathbf{s}}_2``).

The computation is made applying [`trapz`](@ref) (see the 
[Trapz](https://github.com/francescoalemanno/Trapz.jl) Julia package) to
the integrand function `integrand_ξ_GNC_LocalGP_IntegratedGP`.


## Inputs

- `s1` and `s2`: comovign distances where the function must be evaluated

- `y`: the cosine of the angle between the two points `P1` and `P2`

- `cosmo::Cosmology`: cosmology to be used in this computation


## Optional arguments 

- `en::Float64 = 1e6`: just a float number used in order to deal better 
  with small numbers;

- `N_χs::Int = 100`: number of points to be used for sampling the integral
  along the ranges `(0, s1)` (for `χ1`) and `(0, s1)` (for `χ2`); it has been checked that
  with `N_χs ≥ 50` the result is stable.


See also: [`integrand_ξ_GNC_LocalGP_IntegratedGP`](@ref), [`int_on_mu_LocalGP_IntegratedGP`](@ref)
[`integral_on_mu`](@ref), [`ξ_GNC_multipole`](@ref)
"""
function ξ_GNC_LocalGP_IntegratedGP(s1, s2, y, cosmo::Cosmology;
     en::Float64 = 1e6, N_χs::Int = 100)

     χ2s = s2 .* range(1e-6, 1, length = N_χs)

     P1, P2 = GaPSE.Point(s1, cosmo), GaPSE.Point(s2, cosmo)
     IPs = [GaPSE.Point(x, cosmo) for x in χ2s]

     int_ξs = [
          en * GaPSE.integrand_ξ_GNC_LocalGP_IntegratedGP(IP, P1, P2, y, cosmo)
          for IP in IPs
     ]

     res = trapz(χ2s, int_ξs)
     #println("res = $res")
     return res / en
end




##########################################################################################92

##########################################################################################92

##########################################################################################92



function ξ_GNC_IntegratedGP_LocalGP(s1, s2, y, cosmo::Cosmology; kwargs...)
    ξ_GNC_LocalGP_IntegratedGP(s2, s1, y, cosmo; kwargs...)
end

