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
     integrand_ξ_lensing(
          IP1::Point, IP2::Point,
          P1::Point, P2::Point,
          y, cosmo::Cosmology;
          Δχ_min::Float64 = 1e-4) :: Float64

Return the integrand of the lensing auto-correlation function 
``\xi^{v_{\parallel}\kappa} (s_1, s_2, \cos{\theta})``, i.e. the function 
``f(s_1, s_2, y, \chi_1, \chi_2)`` defined as follows:  

```math
f(s_1, s_2, y, \chi_1, \chi_2) = 
     \mathcal{H}_0^2 \Omega_{M0} D(s_2) f(s_2) \mathcal{H}(s_2) \mathcal{R}(s_2) 
     \frac{ D(\chi_1) (\chi_1 - s_1) }{a(\chi_1) s_1} 
     \left(
          J_{00} I^0_0(\chi) + J_{02} I^0_2(\chi) + J_{04} I^0_4(\chi) + J_{20} I^2_0(\chi)
     \right)
```

where ``D_1 = D(\chi_1)``, ``D_2 = D(\chi_2)`` and so on, ``\mathcal{H} = a H``, 
``\chi = \sqrt{\chi_1^2 + \chi_2^2 - 2\chi_1\chi_2\cos{\theta}}``, 
``y = \cos{\theta} = \hat{\mathbf{s}}_1 \dot \hat{\mathbf{s}}_2``) 
and the ``J`` coefficients are given by 

```math
\begin{align*}
     J_{00} & = \frac{1}{15}(\chi_1^2 y + \chi_1(4 y^2 - 3) s_2 - 2 y s_2^2) \\
     J_{02} & = \frac{1}{42 \chi^2} 
          (4 \chi_1^4 y + 4 \chi_1^3 (2 y^2 - 3) s_2 + \chi_1^2 y (11 - 23 y^2) s_2^2 + 
          \chi_1 (23 y^2 - 3) s_2^3 - 8 y s_2^4) \\
     J_{04} & = \frac{1}{70 \chi^2}
          (2 \chi_1^4 y + 2 \chi_1^3 (2y^2 - 3) s_2 - \chi_1^2 y (y^2 + 5) s_2^2 + 
          \chi_1 (y^2 + 9) s_2^3 - 4 y s_2^4) \\
     J_{20} & = y \chi^2
\end{align*}
```

## Inputs

- `IP1::Point` and `IP2::Point`: `Point` inside the integration limits, placed 
  at comoving distance `χ1` and `χ2` respectively.

- `P1::Point` and `P2::Point`: extreme `Point` of the integration, placed 
  at comoving distance `s1` and `s2` respectively.

- `y`: the cosine of the angle between the two points `P1` and `P2`

- `cosmo::Cosmology`: cosmology to be used in this computation


## Optional arguments

- `Δχ_min::Float64 = 1e-6` : when ``\Delta\chi = \sqrt{\chi_1^2 + \chi_2^2 - 2 \, \chi_1 \chi_2 y} \to 0^{+}``,
  some ``I_\ell^n`` term diverges, but the overall parenthesis has a known limit:

  ```math
     \lim_{\chi\to0^{+}} (J_{00} \, I^0_0(\chi) + J_{02} \, I^0_2(\chi) + 
          J_{31} \, I^3_1(\chi) + J_{22} \, I^2_2(\chi)) = 
          \frac{4}{15} \, (5 \, \sigma_2 + \frac{2}{3} \, σ_0 \,s_1^2 \, \chi_2^2)
  ```

  So, when it happens that ``\chi < \Delta\chi_\mathrm{min}``, the function considers this limit
  as the result of the parenthesis instead of calculating it in the normal way; it prevents
  computational divergences.


See also: [`ξ_lensing`](@ref), [`integrand_on_mu_lensing`](@ref)
[`integral_on_mu`](@ref), [`ξ_multipole`](@ref)
"""
function integrand_ξ_lensingdoppler(
     IP::Point, P1::Point, P2::Point,
     y, cosmo::Cosmology)

     s1 = P1.comdist, P1.D
     s2, D_s2, f_s2, ℋ_s2, ℛ_s2 = P2.comdist, P2.D, P2.f, P2.ℋ, P2.ℛ
     χ1, D1, a_χ1 = IP.comdist, IP.D, IP.a
     Ω_M0 = cosmo.params.Ω_M0

     Δχ = √(χ1^2 + χ2^2 - 2 * χ1 * χ2 * y)

     factor = D1 * (χ1 - s1) / (s1 * a_χ1)
     prefactor = ℋ0^2 * Ω_M0^2 * D_s2 * f_s2 * ℋ_s2 * ℛ_s2
     #println("factor = $factor")
     #println("denomin = $denomin")
     #factor = ℋ0^4 * Ω_M0^2 * D1 * (χ1 - s1) * D2 * (χ2 - s2) * psb1 * psb2

     new_J00 = 1.0 / 15.0 * (χ1^2 * y + χ1 * (4 * y^2 - 3) * s2 - 2 * y * s2^2)
     new_J02 = 1.0 / (42 * Δχ^2) * (
          4 * χ1^4 * y + 4 * χ1^3 * (2 * y^2 - 3) * s2
          +
          χ1^2 * y * (11 - 23 * y^2) * s2^2 + χ1 * (23 * y^2 - 3) * s2^2 - 8 * y * s2^4)
     new_J04 = 1.0 / (70 * Δχ^2) *(
          2 * χ1^4 * y + 2 * χ1^3 * (2 * y^2 - 3) * s2
          -
          χ1^2 * y * (y^2 + 5) * s2^2 + χ1 * (y^2 + 9) * s2^2 - 4 * y * s2^4)

     new_J20 = y * Δχ^2

     I00 = cosmo.tools.I00(Δχ)
     I20 = cosmo.tools.I20(Δχ)
     I40 = cosmo.tools.I40(Δχ)
     I02 = cosmo.tools.I02(Δχ)

     #println("J00 = $new_J00, \t I00(Δχ) = $(I00)")
     #println("J02 = $new_J02, \t I20(Δχ) = $(I20)")
     #println("J31 = $new_J31, \t I13(Δχ) = $(I13)")
     #println("J22 = $new_J22, \t I22(Δχ) = $(I22)")

     parenth = (
          new_J00 * I00 + new_J02 * I20 +
          new_J04 * I40 + new_J20 * I02
     )


     res =  prefactor * factor * parenth
     #println("res = ", res, "\n")
     return res
end

#=
function func_Δχ_min(s1, s2, y; frac = 1e-4)
     Δs = s(s1, s2, y)
     return frac * Δs
end
=#

#=
@doc raw"""
     ξ _lensing(s1, s2, y, cosmo::Cosmology;
          enhancer::Float64 = 1.0,
          Δχ_min::Float64 = 1e-3,
          N_χs::Integer = 100) :: Float64

Return the lensing auto-correlation function 
``\xi^{\kappa\kappa} (s_1, s_2, \cos{\theta})``, defined as follows:
    
```math
\xi^{\kappa\kappa} (s_1, s_2, \cos{\theta}) = 
\int_0^{s_1} \mathrm{d} \chi_1 \int_0^{s_2} \mathrm{d} \chi_2 
\frac{1}{2}
\frac{
     \mathcal{H}_0^4 \Omega_{ \mathrm{M0}}^2 D_1 D_2 (\chi_1 - s_1)(\chi_2 - s_2)
}{
     s_1 s_2 a(\chi_1) a(\chi_2) }
(J_{00} \, I^0_0(\chi) + J_{02} \, I^0_2(\chi) + 
     J_{31} \, I^3_1(\chi) + J_{22} \, I^2_2(\chi))
```

where ``D_1 = D(\chi_1)``, ``D_2 = D(\chi_2)`` and so on, ``\mathcal{H} = a H``, 
``\chi = \sqrt{\chi_1^2 + \chi_2^2 - 2\chi_1\chi_2\cos{\theta}}``, 
``y = \cos{\theta} = \hat{\mathbf{s}}_1 \dot \hat{\mathbf{s}}_2``) 
and the ``J`` coefficients are given by 

```math
\begin{align*}
    J_{00} & = - \frac{3 \chi_1^2 \chi_2^2}{4 \chi^4} (y^2 - 1) 
               (8 y (\chi_1^2 + \chi_2^2) - 9 \chi_1 \chi_2 y^2 - 7 \chi_1 \chi_2) \\
    J_{02} & = - \frac{3 \chi_1^2 \chi_2^2}{2 \chi^4} (y^2 - 1)
               (4 y (\chi_1^2 + \chi_2^2) - 3 \chi_1 \chi_2 y^2 - 5 \chi_1 \chi_2) \\
    J_{31} & = 9 y \chi^2 \\
    J_{22} & = \frac{9 \chi_1 \chi_2}{4 \chi^4}
               [ 2 (\chi_1^4 + \chi_2^4) (7 y^2 - 3) 
                 - 16 y \chi_1 \chi_2 (\chi_1^2 + \chi_2^2) (y^2+1) 
               + \chi_1^2 \chi_2^2 (11 y^4 + 14 y^2 + 23)]
\end{align*}
```

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

- `Δχ_min::Float64 = 1e-6` : when ``\Delta\chi = \sqrt{\chi_1^2 + \chi_2^2 - 2 \, \chi_1 \chi_2 y} \to 0^{+}``,
  some ``I_\ell^n`` term diverges, but the overall parenthesis has a known limit:

  ```math
     \lim_{\chi\to0^{+}} (J_{00} \, I^0_0(\chi) + J_{02} \, I^0_2(\chi) + 
          J_{31} \, I^3_1(\chi) + J_{22} \, I^2_2(\chi)) = 
          \frac{4}{15} \, (5 \, \sigma_2 + \frac{2}{3} \, σ_0 \,s_1^2 \, \chi_2^2)
  ```

  So, when it happens that ``\chi < \Delta\chi_\mathrm{min}``, the function considers this limit
  as the result of the parenthesis instead of calculating it in the normal way; it prevents
  computational divergences.

- `N_χs::Integer = 100`: number of points to be used for sampling the integral
  along the ranges `(0, s1)` (for `χ1`) and `(0, s1)` (for `χ2`); it has been checked that
  with `N_χs ≥ 50` the result is stable.


See also: [`integrand_ξ_lensing`](@ref), [`integrand_on_mu_lensing`](@ref)
[`integral_on_mu`](@ref), [`ξ_multipole`](@ref)
"""
function ξ_lensing(s1, s2, y, cosmo::Cosmology;
     enhancer::Float64 = 1.0, N_χs::Integer = 100, Δχ_min::Float64 = 1e-3)

     adim_χs = range(1e-6, 1.0, N_χs)
     #Δχ_min = func_Δχ_min(s1, s2, y; frac = frac_Δχ_min)

     P1, P2 = GaPSE.Point(s1, cosmo), GaPSE.Point(s2, cosmo)
     χ1s = adim_χs .* s1
     χ2s = adim_χs .* s2

     IP1s = [GaPSE.Point(x, cosmo) for x in χ1s]
     IP2s = [GaPSE.Point(x, cosmo) for x in χ2s]

     int_ξ_lensings = [
          GaPSE.integrand_ξ_lensing(IP1, IP2, P1, P2, y, cosmo;
               enhancer = enhancer, Δχ_min = Δχ_min)
          for IP1 in IP1s, IP2 in IP2s
     ]

     res = trapz((χ1s, χ2s), int_ξ_lensings)
     #println("res = $res")
     return res
end
=#
