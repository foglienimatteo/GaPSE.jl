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



function integrand_ξ_GNC_Lensing_Doppler(
     IP::Point, P1::Point, P2::Point,
     y, cosmo::Cosmology; obs::Union{Bool, Symbol} = :noobsvel,
     Δχ_min::Float64 = 1e-1)

     s1 = P1.comdist
     s2, D_s2, f_s2, ℋ_s2, ℛ_s2 = P2.comdist, P2.D, P2.f, P2.ℋ, P2.ℛ_GNC
     χ1, D1, a1 = IP.comdist, IP.D, IP.a
     s_b_s1, s_b_s2 = cosmo.params.s_b, cosmo.params.s_b
     Ω_M0 = cosmo.params.Ω_M0


     Δχ1_square = χ1^2 + s2^2 - 2 * χ1 * s2 * y
     Δχ1 = Δχ1_square > 0 ? √(Δχ1_square) : 0

     common = ℋ0^2 * Ω_M0 * D1 * (χ1 - s1) * (5 * s_b_s1 - 2) / (s1 * a1)
     factor = D_s2 * f_s2 * ℋ_s2 * ℛ_s2

     first_part = if Δχ1 ≥ Δχ_min
          new_J00 = 1 / 15 * (χ1^2 * y + χ1 * s2 * (4 * y^2 - 3) - 2 * y * s2^2)
          new_J02 = 1 / (42 * Δχ1^2) * (
               4 * χ1^4 * y + 4 * χ1^3 * (2 * y^2 - 3) * s2
               + χ1^2 * y * (11 - 23 * y^2) * s2^2
               + χ1 * (23 * y^2 - 3) * s2^3 - 8 * y * s2^4)
          new_J04 = 1 / (70 * Δχ1^2) * (
               2 * χ1^4 * y + 2 * χ1^3 * (2 * y^2 - 3) * s2
               -
               χ1^2 * y * (y^2 + 5) * s2^2
               +
               χ1 * (y^2 + 9) * s2^3 - 4 * y * s2^4)
          new_J20 = y * Δχ1^2

          I00 = cosmo.tools.I00(Δχ1)
          I20 = cosmo.tools.I20(Δχ1)
          I40 = cosmo.tools.I40(Δχ1)
          I02 = cosmo.tools.I02(Δχ1)

          common * factor * (
                 new_J00 * I00 + new_J02 * I20 +
                 new_J04 * I40 + new_J20 * I02
            )
     else
          common * factor * cosmo.tools.σ_2
     end


     if obs == false || obs == :no || obs == :noobsvel
          return first_part
            
     elseif obs == true || obs == :yes
          #### New observer terms #########

          I13_χ1 = cosmo.tools.I13(χ1)

          obs_terms = - 3 * χ1^2 * y * f0 * ℋ0 * (ℛ_s2 - 5 * s_b_s2 + 2) * common * I13_χ1

          #################################     
          
          return first_part + obs_terms
     else
          throw(AssertionError(":$obs is not a valid Symbol for \"obs\"; they are: \n\t"*
               "$(":".*string.(VALID_OBS_VALUES) .* vcat([" , " for i in 1:length(VALID_OBS_VALUES)-1], " .")... )" 
               ))
     end

end


function integrand_ξ_GNC_Lensing_Doppler(
     χ1::Float64, s1::Float64, s2::Float64,
     y, cosmo::Cosmology; kwargs...)

     P1, P2 = Point(s1, cosmo), Point(s2, cosmo)
     IP = Point(χ1, cosmo)
     return integrand_ξ_GNC_Lensing_Doppler(IP, P1, P2, y, cosmo; kwargs...)
end


##########################################################################################92


"""
     ξ_GNC_Lensing_Doppler(s1, s2, y, cosmo::Cosmology;
          en::Float64 = 1e6, N_χs::Int = 100):: Float64

Return the Lensing-Doppler cross-correlation function 
``\\xi^{v_{\\parallel}\\kappa} (s_1, s_2, \\cos{\\theta})`` concerning the perturbed
luminosity distance, defined as follows:
    
```math
\\xi^{v_{\\parallel}\\kappa} (s_1, s_2, \\cos{\\theta}) = 
     \\mathcal{H}_0^2 \\Omega_{M0} D(s_2) f(s_2) \\mathcal{H}(s_2) \\mathcal{R}(s_2) 
     \\int_0^{s_1} \\mathrm{d} \\chi_1 
     \\frac{ D(\\chi_1) (\\chi_1 - s_1) }{a(\\chi_1) s_1} 
     \\left(
          J_{00} I^0_0(\\Delta\\chi_1) + J_{02} I^0_2(\\Delta\\chi_1) 
          + J_{04} I^0_4(\\Delta\\chi_1) + J_{20} I^2_0(\\Delta\\chi_1)
     \\right)
```

where ``\\mathcal{H} = a H``, 
``\\Delta\\chi_1= \\sqrt{\\chi_1^2 + s_2^2 - 2 \\chi_1 s_2 \\cos{\\theta}}``, 
``y = \\cos{\\theta} = \\hat{\\mathbf{s}}_1 \\cdot \\hat{\\mathbf{s}}_2``) 
and the ``J`` coefficients are given by:

```math
\\begin{align*}
     J_{00} & = \\frac{1}{15}(\\chi_1^2 y + \\chi_1(4 y^2 - 3) s_2 - 2 y s_2^2) \\\\
     J_{02} & = \\frac{1}{42 \\Delta\\chi_1^2} 
          (4 \\chi_1^4 y + 4 \\chi_1^3 (2 y^2 - 3) s_2 + \\chi_1^2 y (11 - 23 y^2) s_2^2 + 
          \\chi_1 (23 y^2 - 3) s_2^3 - 8 y s_2^4) \\\\
     J_{04} & = \\frac{1}{70 \\Delta\\chi_1^2}
          (2 \\chi_1^4 y + 2 \\chi_1^3 (2y^2 - 3) s_2 - \\chi_1^2 y (y^2 + 5) s_2^2 + 
          \\chi_1 (y^2 + 9) s_2^3 - 4 y s_2^4) \\\\
     J_{20} & = y \\Delta\\chi_1^2
\\end{align*}
```

The computation is made applying [`trapz`](@ref) (see the 
[Trapz](https://github.com/francescoalemanno/Trapz.jl) Julia package) to
the integrand function `integrand_ξ_GNC_GNC_Lensing_Doppler`.


## Inputs

- `s1` and `s2`: comovign distances where the function must be evaluated

- `y`: the cosine of the angle between the two points `P1` and `P2`

- `cosmo::Cosmology`: cosmology to be used in this computation


## Optional arguments 

- `en::Float64 = 1e6`: just a float number used in order to deal better 
  with small numbers;

- `Δχ_min::Float64 = 1e-6` : when ``\\Delta\\chi = \\sqrt{\\chi_1^2 + \\chi_2^2 - 2 \\, \\chi_1 \\chi_2 y} \\to 0^{+}``,
  some ``I_\\ell^n`` term diverges, but the overall parenthesis has a known limit:

  ```math
     \\lim_{\\chi\\to0^{+}} (J_{00} \\, I^0_0(\\chi) + J_{02} \\, I^0_2(\\chi) + 
          J_{31} \\, I^3_1(\\chi) + J_{22} \\, I^2_2(\\chi)) = 
          \\frac{4}{15} \\, (5 \\, \\sigma_2 + \\frac{2}{3} \\, σ_0 \\,s_1^2 \\, \\chi_2^2)
  ```

  So, when it happens that ``\\chi < \\Delta\\chi_\\mathrm{min}``, the function considers this limit
  as the result of the parenthesis instead of calculating it in the normal way; it prevents
  computational divergences.

- `N_χs::Int = 100`: number of points to be used for sampling the integral
  along the ranges `(0, s1)` (for `χ1`) and `(0, s1)` (for `χ2`); it has been checked that
  with `N_χs ≥ 50` the result is stable.


See also: [`integrand_ξ_GNC_Lensing_Doppler`](@ref), [`int_on_mu_Lensing_Doppler`](@ref)
[`integral_on_mu`](@ref), [`ξ_GNC_multipole`](@ref)
"""
function ξ_GNC_Lensing_Doppler(s1, s2, y, cosmo::Cosmology;
     en::Float64 = 1e6, N_χs::Int = 100, obs::Union{Bool, Symbol} = :noobsvel, 
     Δχ_min::Float64 = 1e-1)

     χ1s = s1 .* range(1e-6, 1, length = N_χs)

     P1, P2 = GaPSE.Point(s1, cosmo), GaPSE.Point(s2, cosmo)
     IPs = [GaPSE.Point(x, cosmo) for x in χ1s]

     int_ξs = [
          en * GaPSE.integrand_ξ_GNC_Lensing_Doppler(IP, P1, P2, y, cosmo; 
               obs = obs, Δχ_min = Δχ_min)
          for IP in IPs
     ]

     res = trapz(χ1s, int_ξs)
     #println("res = $res")
     return res / en
end





##########################################################################################92

##########################################################################################92

##########################################################################################92




"""
     ξ_GNC_Doppler_Lensing(s1, s2, y, cosmo::Cosmology; kwargs...) = 
          ξ_GNC_Lensing_Doppler(s2, s1, y, cosmo; kwargs...)

Return the Two-Point Correlation Function (TPCF) given by the cross correlation between the 
Doppler and the Lensing effects arising from the Galaxy Number Counts (GNC).

It's computed through the symmetric function `ξ_GNC_Lensing_Doppler`; check its documentation for
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
[`ξ_GNC_Lensing_Doppler`](@ref)
"""
function ξ_GNC_Doppler_Lensing(s1, s2, y, cosmo::Cosmology; kwargs...)
     ξ_GNC_Lensing_Doppler(s2, s1, y, cosmo; kwargs...)
end

#=
function integrand_ξ_GNC_Doppler_Lensing(
     χ1::Float64, s1::Float64, s2::Float64,
     y, cosmo::Cosmology; kwargs...)

     P1, P2 = Point(s1, cosmo), Point(s2, cosmo)
     IP = Point(χ1, cosmo)
     return integrand_ξ_GNC_Lensing_Doppler(IP, P2, P1, y, cosmo; kwargs...)
end
=#
