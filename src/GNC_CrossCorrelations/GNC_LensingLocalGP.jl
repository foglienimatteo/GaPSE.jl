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



function integrand_ξ_GNC_Lensing_LocalGP(
     IP::Point, P1::Point, P2::Point,
     y, cosmo::Cosmology; obs::Union{Bool, Symbol} = :noobsvel)

     s1 = P1.comdist
     s2, D_s2, f_s2, a_s2, ℋ_s2, ℛ_s2 = P2.comdist, P2.D, P2.f, P2.a, P2.ℋ, P2.ℛ_GNC
     χ1, D1, a1 = IP.comdist, IP.D, IP.a
     s_b_s1, s_b_s2 = cosmo.params.s_b, cosmo.params.s_b
     𝑓_evo_s2 = cosmo.params.𝑓_evo
     Ω_M0 = cosmo.params.Ω_M0

     Δχ1_square = χ1^2 + s2^2 - 2 * χ1 * s2 * y
     Δχ1 = Δχ1_square > 0 ? √(Δχ1_square) : 0

     common = D_s2 * ℋ0^2 * Ω_M0 * s2 * D1 * (χ1 - s1) *  (5 * s_b_s1 - 2) * (
                   2 * f_s2 * a_s2 * ℋ_s2^2 * (𝑓_evo_s2 - 3)
                   + 3 * ℋ0^2 * Ω_M0 * (f_s2 + ℛ_s2 + 5 * s_b_s2 - 2)
              ) / (a1 * a_s2 * s1)
     factor = 2 * y * χ1^2 - χ1 * s2 * (y^2 + 3) + 2 * y * s2^2

     J20 = 1/2 * y * Δχ1^2

     I00 = cosmo.tools.I00(Δχ1)
     I20 = cosmo.tools.I20(Δχ1)
     I40 = cosmo.tools.I40(Δχ1)
     I02 = cosmo.tools.I02(Δχ1)

     return common * (
                factor * (1 / 60 * I00 + 1 / 42 * I20 + 1 / 140 * I40)
                + J20 * I02
           )
end


function integrand_ξ_GNC_Lensing_LocalGP(
     χ1::Float64, s1::Float64, s2::Float64,
     y, cosmo::Cosmology; obs::Union{Bool,Symbol}=:noobsvel)

     P1, P2 = Point(s1, cosmo), Point(s2, cosmo)
     IP = Point(χ1, cosmo)
     return integrand_ξ_GNC_Lensing_LocalGP(IP, P1, P2, y, cosmo; obs=obs)
end




##########################################################################################92




function ξ_GNC_Lensing_LocalGP(s1, s2, y, cosmo::Cosmology;
     en::Float64 = 1e6, N_χs::Int = 100, obs::Union{Bool, Symbol} = :noobsvel)

     χ1s = s1 .* range(1e-6, 1, length = N_χs)

     P1, P2 = GaPSE.Point(s1, cosmo), GaPSE.Point(s2, cosmo)
     IPs = [GaPSE.Point(x, cosmo) for x in χ1s]

     int_ξs = [
          en * GaPSE.integrand_ξ_GNC_Lensing_LocalGP(IP, P1, P2, y, cosmo; obs = obs)
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
     ξ_GNC_LocalGP_Lensing(s1, s2, y, cosmo::Cosmology; kwargs...) = 
          ξ_GNC_Lensing_LocalGP(s2, s1, y, cosmo; kwargs...)

Return the Two-Point Correlation Function (TPCF) given by the cross correlation between the 
Local Gravitational Potential (GP) and the Lensing effects arising from the Galaxy Number Counts (GNC).

It's computed through the symmetric function `ξ_GNC_Lensing_LocalGP`; check its documentation for
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
[`ξ_GNC_Lensing_LocalGP`](@ref)
"""
function ξ_GNC_LocalGP_Lensing(s1, s2, y, cosmo::Cosmology; kwargs...)
     ξ_GNC_Lensing_LocalGP(s2, s1, y, cosmo; kwargs...)
end

