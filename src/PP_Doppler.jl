
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

function ξ_PP_Doppler_L0(P::Point, cosmo::Cosmology)
     P1 = Point(cosmo.s_eff, cosmo)
     D, f, ℋ, ℛ = P1.D, P1.f, P1.ℋ, P1.ℛ
     s = P.comdist

     1.0 / 3.0 * f^2 * ℋ^2 * ℛ^2 * D^2 * s^2 * cosmo.tools.I02(s)
end

function ξ_PP_Doppler_L0(s1, cosmo::Cosmology)
     P1 = Point(s1, cosmo)
     return ξ_PP_Doppler_L0(P1, cosmo)
end

function ξ_PP_Doppler_L2(P::Point, cosmo::Cosmology)
     P1 = Point(cosmo.s_eff, cosmo)
     D, f, ℋ, ℛ = P1.D, P1.f, P1.ℋ, P1.ℛ
     s = P.comdist
     -2.0 / 3.0 * f^2 * ℋ^2 * ℛ^2 * D^2 * s^2 * cosmo.tools.I22(s)
end

function ξ_PP_Doppler_L2(s1, cosmo::Cosmology)
     P1 = Point(s1, cosmo)
     return ξ_PP_Doppler_L2(P1, cosmo)
end

function ξ_PP_Doppler(s1, y, cosmo::Cosmology)
     return ξ_PP_Doppler_L0(s1, cosmo) + ξ_PP_Doppler_L2(s1, cosmo) * Pl(y, 2)
end



##########################################################################################92



function PP_integrand_ξ_multipole(s, μ, cosmo::Cosmology;
     L::Integer=0, use_windows::Bool=true)

     res = if use_windows == true
          ϕ_s = ϕ(s, cosmo.s_min, cosmo.s_max)
          (ϕ_s > 0.0) || (return 0.0)
          #println("s1 = $s1 \\t s2 = $(s2(s1, s, μ)) \\t  y=$(y(s1, s, μ))")
          int = ξ_PP_Doppler(s, μ, cosmo)
          #println("int = $int")
          int .* (ϕ_s * spline_F(s / cosmo.s_eff, μ, cosmo.windowF) * Pl(μ, L))
     else
          #println("s1 = $s1 \\t s2 = $(s2(s1, s, μ)) \\t  y=$(y(s1, s, μ))")
          int = ξ_PP_Doppler(s, μ, cosmo)
          #println("int = $int")
          #println( "Pl(μ, L) = $(Pl(μ, L))")
          int .* Pl(μ, L)
     end

     #println("res = $res")
     return (2.0 * L + 1.0) / 2.0 * res
end

function PP_ξ_multipole(
     s, cosmo::Cosmology;
     L::Integer=0,
     use_windows::Bool=true,
     enhancer::Float64=1e6,
     μ_atol::Float64=0.0,
     μ_rtol::Float64=1e-2)


     orig_f(μ) = enhancer * PP_integrand_ξ_multipole(s, μ, cosmo;
          L=L, use_windows=use_windows)

     int = quadgk(μ -> orig_f(μ), -1.0, 1.0; atol=μ_atol, rtol=μ_rtol)[1]


     return int / enhancer
end
