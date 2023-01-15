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



function integrand_ξ_GNC_Doppler_IntegratedGP(
     IP::Point, P1::Point, P2::Point,
     y, cosmo::Cosmology; obs::Union{Bool, Symbol} = :noobsvel)

     s1, D_s1, f_s1, ℋ_s1, ℛ_s1 = P1.comdist, P1.D, P1.f, P1.ℋ, P1.ℛ_GNC
     s2, ℛ_s2 = P2.comdist, P2.ℛ_GNC
     χ2, D2, a2, f2, ℋ2 = IP.comdist, IP.D, IP.a, IP.f, IP.ℋ
     s_b_s1, s_b_s2 = cosmo.params.s_b, cosmo.params.s_b
     Ω_M0 = cosmo.params.Ω_M0

     Δχ2_square = s1^2 + χ2^2 - 2 * s1 * χ2 * y
     Δχ2 = Δχ2_square > 0 ? √(Δχ2_square) : 0

     common = ℋ0^2 * Ω_M0 * D2 / (s2 * a2)
     factor = Δχ2^2 * f_s1 * ℋ_s1 * ℛ_s1 * (χ2 * y - s1) 
     parenth = s2 * ℋ2 * ℛ_s2 * (f2 - 1) - 5 * s_b_s2 + 2

     I00 = cosmo.tools.I00(Δχ2)
     I20 = cosmo.tools.I20(Δχ2)
     I40 = cosmo.tools.I40(Δχ2)
     I02 = cosmo.tools.I02(Δχ2)

     if obs == false || obs == :no || obs == :noobsvel
          return D_s1 * common * factor * parenth * (
                    1 / 15 * I00 + 2 / 21 * I20
                    + 1 / 35 * I40 + 1 * I02
               )
     elseif obs == true || obs == :yes
          #### New observer terms #########

          I13_χ2 = cosmo.tools.I13(χ2)

          obs_terms = - 3 * χ2^3 * y * f0 * ℋ0 * (ℛ_s1 - 5 * s_b_s1 + 2) * common * parenth * I13_χ2

          #################################
          
          return D_s1 * common * factor * parenth * (
                    1 / 15 * I00 + 2 / 21 * I20
                    + 1 / 35 * I40 + 1 * I02
               ) + obs_terms
     else
          throw(AssertionError(":$obs is not a valid Symbol for \"obs\"; they are: \n\t"*
               "$(":".*string.(VALID_OBS_VALUES) .* vcat([" , " for i in 1:length(VALID_OBS_VALUES)-1], " .")... )" 
               ))
     end
end


function integrand_ξ_GNC_Doppler_IntegratedGP(
     χ2::Float64, s1::Float64, s2::Float64,
     y, cosmo::Cosmology; obs::Union{Bool,Symbol}=:noobsvel)

     P1, P2 = Point(s1, cosmo), Point(s2, cosmo)
     IP = Point(χ2, cosmo)
     return integrand_ξ_GNC_Doppler_IntegratedGP(IP, P1, P2, y, cosmo; obs=obs)
end



##########################################################################################92



function ξ_GNC_Doppler_IntegratedGP(s1, s2, y, cosmo::Cosmology;
     en::Float64 = 1e6, N_χs::Int = 100, obs::Union{Bool, Symbol} = :noobsvel)

     χ2s = s2 .* range(1e-6, 1, length = N_χs)

     P1, P2 = GaPSE.Point(s1, cosmo), GaPSE.Point(s2, cosmo)
     IPs = [GaPSE.Point(x, cosmo) for x in χ2s]

     int_ξs = [
          en * GaPSE.integrand_ξ_GNC_Doppler_IntegratedGP(IP, P1, P2, y, cosmo; obs = obs)
          for IP in IPs
     ]

     res = trapz(χ2s, int_ξs)
     #println("res = $res")
     return res / en
end




##########################################################################################92

##########################################################################################92

##########################################################################################92



"""
     ξ_GNC_IntegratedGP_Doppler(s1, s2, y, cosmo::Cosmology; kwargs...) = 
          ξ_GNC_Doppler_IntegratedGP(s2, s1, y, cosmo; kwargs...)

Return the Two-Point Correlation Function (TPCF) given by the cross correlation between the 
Integrated Gravitational Potential (GP) and the Doppler effects arising from the Galaxy Number Counts (GNC).

It's computed through the symmetric function `ξ_GNC_Doppler_IntegratedGP`; check its documentation for
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
[`ξ_GNC_Doppler_IntegratedGP`](@ref)
"""
function ξ_GNC_IntegratedGP_Doppler(s1, s2, y, cosmo::Cosmology; kwargs...)
     ξ_GNC_Doppler_IntegratedGP(s2, s1, y, cosmo; kwargs...)
end

