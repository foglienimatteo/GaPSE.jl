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

function ξ_integratedGP(P1::Point, P2::Point, y, cosmo::Cosmology; enhancer = 1.0)
     s1, D1, f1, a1, ℛ1 = P1.comdist, P1.D, P1.f, P1.a, P1.ℛ
     s2, D2, f2, a2, ℛ2 = P2.comdist, P2.D, P2.f, P2.a, P2.ℛ

     Δs = s(s1, s2, y)
     prefac = 2.25 * ℋ0^4 * cosmo.params.Ω_M0^2 * D1 * D2 * Δs^4 / (a1 * a2)
     parenth = 1.0 + ℛ1 + ℛ2 + ℛ1 * ℛ2

     I04 = cosmo.tools.I04(Δs)

     res = prefac * I04 * parenth

     return enhancer * res
end



function integrand_on_mu_integratedGP(s1, s, μ,
     cosmo::Cosmology; L::Integer = 0, enhancer = 1.0,
     use_windows::Bool = true)

     s2_value = s2(s1, s, μ)
     y_value = y(s1, s, μ)

     if use_windows == true
          ϕ_s2 = ϕ(s2_value)
          (ϕ_s2 > 0.0) || (return 0.0)
          P1, P2 = Point(s1, cosmo), Point(s2_value, cosmo)
          val = ξ_integratedGP(P1, P2, y_value, cosmo; enhancer = enhancer)
          return val * ϕ_s2 * spline_F(s / s1, μ, cosmo.windowF) * Pl(μ, L)
     else
          P1, P2 = Point(s1, cosmo), Point(s2_value, cosmo)
          val = ξ_integratedGP(P1, P2, y_value, cosmo; enhancer = enhancer)
          return val * Pl(μ, L)
     end
end


##########################################################################################92
