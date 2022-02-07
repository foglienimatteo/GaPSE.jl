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


@doc raw"""
     integrand_ξ_integratedGP(IP1::Point, IP2::Point,
          P1::Point, P2::Point,
          y, cosmo::Cosmology;
          enhancer::Float64 = 1.0) :: Float64

Return the integrand of the integrated gravitational potential 
auto-correlation function ``\xi^{\int\phi\int\phi} (s_1, s_2, \cos{\theta})``, 
i.e. the function ``f(s_1, s_2, y, \chi_1, \chi_2)`` defined as follows:  

```math
f(s_1, s_2, y, \chi_1, \chi_2) = J_{40}(s_1, s_2, y, \chi_1, \chi_2) \tilde{I}^4_0(\chi)
```
where ``\chi = \sqrt{\chi_1^2 + \chi_2^2 - 2 \, \chi_1 \, \chi_2 \, y} ``,
``y = \cos{\theta} = \hat{\mathbf{s}}_1 \dot \hat{\mathbf{s}}_2`` and:
```math
\begin{split}
     &J_{40}(s_1, s_2, y, \chi_1, \chi_2)  = 
          \frac{
                    9 \mathcal{H}_0^4 \Omega_{M0}^2 D(\chi_1) D(\chi_2) \chi^4
          }{    a(\chi_1) a(\chi_2) s_1 s_2} 
          (s_2 \mathcal{H}(\chi_2) \mathcal{R}(s_2) (f(\chi_2)-1) - 1) 
          (s_1 \mathcal{H}(\chi_1) \mathcal{R}(s_1) (f(\chi_1)-1) - 1)\\[5pt]
     &\tilde{I}^4_0 (s) = \int_0^\infty \frac{\mathrm{d}q}{2\pi^2} 
          q^2 \, P(q) \, \frac{j_0(q s) - 1}{(q s)^4}
\end{split}
```


## Inputs

- `IP1::Point` and `IP2::Point`: `Point` inside the integration limits, placed 
  at comoving distance `χ1` and `χ2` respectively.

- `P1::Point` and `P2::Point`: extreme `Point` of the integration, placed 
  at comoving distance `s1` and `s2` respectively.

- `y`: the cosine of the angle between the two points `P1` and `P2`

- `cosmo::Cosmology`: cosmology to be used in this computation


## Optional arguments 

- `enhancer::Float64 = 1.0` : multiply the resulting ``f(s_1, s_2, y, \chi_1, \chi_2)`` value; it
  is very useful for interal computations in other functions (for instance 
  `map_integral_on_mu_lensing`), in order to deal better with small float numbers.


See also: [`ξ_integratedGP`](@ref), [`integrand_on_mu_integratedGP`](@ref)
[`integral_on_mu`](@ref), [`ξ_multipole`](@ref)
"""
function integrand_ξ_integratedGP(IP1::Point, IP2::Point,
     P1::Point, P2::Point,
     y, cosmo::Cosmology;
     enhancer::Float64 = 1.0)

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
end



@doc raw"""
     ξ_integratedGP(P1::Point, P2::Point, y, cosmo::Cosmology; 
          enhancer::Float64 = 1.0, 
          N_χs::Integer = 100) :: Float64

Return the integrated gravitational potential auto-correlation function 
``\xi^{\int\phi\int\phi}(s_1, s_2, \cos{\theta})``, defined as follows:
    
```math
\xi^{\int\phi\int\phi} (s_1, s_2, \cos{\theta}) = 
     \int_0^{s_1} \mathrm{d} \chi_1 \int_0^{s_2}\mathrm{d} \chi_2 \;
     J_{40}(s_1, s_2, y, \chi_1, \chi_2) \, \tilde{I}^4_0(\chi)
```
where ``\chi = \sqrt{\chi_1^2 + \chi_2^2 - 2 \, \chi_1 \, \chi_2 \, y} ``,
``y = \cos{\theta} = \hat{\mathbf{s}}_1 \dot \hat{\mathbf{s}}_2`` and:
```math
\begin{align*}
     J_{40}(s_1, s_2, y, \chi_1, \chi_2) &= 
          \frac{
                    9 \mathcal{H}_0^4 \Omega_{M0}^2 D(\chi_1) D(\chi_2) \chi^4
          }{    a(\chi_1) a(\chi_2) s_1 s_2} 
          (s_2 \mathcal{H}(\chi_2) \mathcal{R}(s_2) (f(\chi_2)-1) - 1) 
          (s_1 \mathcal{H}(\chi_1) \mathcal{R}(s_1) (f(\chi_1)-1) - 1) \\[5pt]
     \tilde{I}^4_0 (s) &= \int_0^\infty \frac{\mathrm{d}q}{2\pi^2} 
          q^2 \, P(q) \, \frac{j_0(q s) - 1}{(q s)^4}
\end{aling*}
```
and ``P(q)`` is the input power spectrum.


The computation is made applying [`trapz`](@ref) (see the 
[Trapz](https://github.com/francescoalemanno/Trapz.jl) Julia package) to
the integrand function `integrand_ξ_lensing`.


## Inputs

- `s1` and `s2`: comovign distances where the function must be evaluated

- `y`: the cosine of the angle between the two points `P1` and `P2`

- `cosmo::Cosmology`: cosmology to be used in this computation


## Optional arguments 

- `enhancer::Float64 = 1.0` : multiply the resulting value; it
  is very useful for interal computations in other functions (for instance 
  `map_integral_on_mu`), in order to deal better with small float numbers.

- `N_χs::Integer = 100`: number of points to be used for sampling the integral
  along the ranges `(0, s1)` (for `χ1`) and `(0, s1)` (for `χ2`); it has been checked that
  with `N_χs ≥ 50` the result is stable.


See also: [`integrand_ξ_integratedGP`](@ref), [`integrand_on_mu_integratedGP`](@ref)
[`integral_on_mu`](@ref), [`ξ_multipole`](@ref)
"""
function ξ_integratedGP(P1::Point, P2::Point, y, cosmo::Cosmology;
     enhancer::Float64 = 1.0, N_χs::Integer = 100)
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
