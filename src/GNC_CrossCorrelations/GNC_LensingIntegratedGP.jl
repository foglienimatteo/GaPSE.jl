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
     y, cosmo::Cosmology; obs::Union{Bool, Symbol} = :noobsvel)

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

#=
function func_Δχ_min(s1, s2, y; frac = 1e-4)
     Δs = s(s1, s2, y)
     return frac * Δs
end
=#


##########################################################################################92



function ξ_GNC_Lensing_IntegratedGP(P1::Point, P2::Point, y, cosmo::Cosmology;
     en::Float64 = 1e6, N_χs_2::Int = 100, obs::Union{Bool, Symbol} = :noobsvel)

     χ1s = P1.comdist .* range(1e-6, 1, length = N_χs_2)
     χ2s = P2.comdist .* range(1e-6, 1, length = N_χs_2)

     IP1s = [GaPSE.Point(x, cosmo) for x in χ1s]
     IP2s = [GaPSE.Point(x, cosmo) for x in χ2s]

     int_ξs = [
          en * GaPSE.integrand_ξ_GNC_Lensing_IntegratedGP(IP1, IP2, P1, P2, y, cosmo; obs = obs)
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

