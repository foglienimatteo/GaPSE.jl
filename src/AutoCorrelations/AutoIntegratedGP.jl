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


function integrand_ξ_integratedGP(IP1::Point, IP2::Point,
     P1::Point, P2::Point,
     y, cosmo::Cosmology;
     enhancer::Float64 = 1.0, Δχ_min = 0.0)

     s1, ℛ_s1 = P1.comdist, P1.ℛ
     s2, ℛ_s2 = P2.comdist, P2.ℛ
     χ1, D1, a_χ1, ℋ1, f1 = IP1.comdist, IP1.D, IP1.a, IP1.ℋ, IP1.f
     χ2, D2, a_χ2, ℋ2, f2 = IP2.comdist, IP2.D, IP2.a, IP2.ℋ, IP2.f
     Ω_M0 = cosmo.params.Ω_M0

     Δχ = √(χ1^2 + χ2^2 - 2 * χ1 * χ2 * y)

     denomin = s1 * s2 * a_χ1 * a_χ2
     factor = ℋ0^4 * Ω_M0^2 * D1 * D2 * Δχ^4
     par_1 = s1 * ℋ1 * ℛ_s1 * (f1 - 1) - 1
     par_2 = s2 * ℋ2 * ℛ_s2 * (f2 - 1) - 1
     #println("factor = $factor")
     #println("denomin = $denomin")
     

     return enhancer * factor / denomin * par_1 * par_2

     #=
     Δχ_min = 0.0
     first_res = if Δχ > Δχ_min
          χ1χ2 = χ1 * χ2

          new_J00 = -0.75 * χ1χ2^2 / Δχ^4 * (y^2 - 1) * (8 * y * (χ1^2 + χ2^2) - χ1χ2 * (9 * y^2 + 7))
          new_J02 = -1.5 * χ1χ2^2 / Δχ^4 * (y^2 - 1) * (4 * y * (χ1^2 + χ2^2) - χ1χ2 * (3 * y^2 + 5))
          new_J31 = 9 * y * Δχ^2
          new_J22 = 2.25 * χ1χ2 / Δχ^4 * (
               2 * (χ1^4 + χ2^4) * (7 * y^2 - 3)
               -
               16 * y * χ1χ2 * (y^2 + 1) * (χ1^2 + χ2^2)
               +
               χ1χ2^2 * (11y^4 + 14y^2 + 23)
          )

          I00 = cosmo.tools.I00(Δχ)
          I20 = cosmo.tools.I20(Δχ)
          I13 = cosmo.tools.I13(Δχ)
          I22 = cosmo.tools.I22(Δχ)

          #println("J00 = $new_J00, \t I00(Δχ) = $(I00)")
          #println("J02 = $new_J02, \t I20(Δχ) = $(I20)")
          #println("J31 = $new_J31, \t I13(Δχ) = $(I13)")
          #println("J22 = $new_J22, \t I22(Δχ) = $(I22)")

          (
               new_J00 * I00 + new_J02 * I20 +
               new_J31 * I13 + new_J22 * I22
          )
     else

          lim = 4.0 / 15.0 * (5.0 * cosmo.tools.σ_2 + 6.0 * cosmo.tools.σ_0 * χ2^2)
          #println("lim = $lim")
          9.0 / 4.0 * lim
     end

     res = enhancer * factor / denomin * first_res
     #println("res = ", res, "\n")
     return res
     =#
end



@doc raw"""
```math
\xi^{\int\phi\int\phi} (s_1, s_2, \cos{\theta}) = 
     \int_0^{s_1} \mathrm{d} \chi_1 \int_0^{s_2}\mathrm{d} \chi_2 
     J_{40} \tilde{I}^4_0(\chi)
```
where
```math
J_{40} = 
    \frac{
          9 \mathcal{H}_0^4 \Omega_{M0}^2 D(\chi_1) D(\chi_2) \chi^4
    }{    a(\chi_1) a(\chi_2) s_1 s_2} 
    (s_2 \mathcal{H}(\chi_2) \mathcal{R}(s_2) (f(\chi_2)-1) - 1) 
    (s_1 \mathcal{H}(\chi_1) \mathcal{R}(s_1) (f(\chi_1)-1) - 1)
```

"""
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
