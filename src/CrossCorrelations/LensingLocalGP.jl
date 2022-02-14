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
     integrand_ξ_lensinglocalgp(
          IP::Point, P1::Point, P2::Point,
          y, cosmo::Cosmology) :: Float64

Return the integrand of the Lensing-LocalGP cross-correlation function 
``\xi^{\kappa \phi} (s_1, s_2, \cos{\theta})``, i.e. the function 
``f(s_1, s_2, y, \chi_1, \chi_2)`` defined as follows:  

```math
f(s_1, s_2, y, \chi_1, \chi_2) = 
     \frac{
          9 \mathcal{H}_0^4 \Omega_{M0}^2 D(s_2) (1 + \mathcal{R}(s_2)) s_2
     }{4 a(s_2) s_1} 
     \frac{D(\chi_1)(s_1 - \chi_1) }{a(\chi_1)}
     \left( J_{31} I^3_1(\Delta\chi_1) +  J_{22} I^2_2(\Delta\chi_1) \right)
```

where ``\mathcal{H} = a H``, 
``\Delta\chi_1 = \sqrt{\chi_1^2 + s_2^2 - 2 \chi_1 s_2\cos{\theta}}``, 
``y = \cos{\theta} = \hat{\mathbf{s}}_1 \dot \hat{\mathbf{s}}_2``) 
and the ``J`` coefficients are given by 

```math
\begin{align*}
     J_{31} & = -2 y \Delta\chi_1^2 \\
     J_{22} & = \chi_1 s_2 (1 - y^2)
\end{align*}
```

## Inputs

- `IP::Point`: `Point` inside the integration limits, placed 
  at comoving distance `χ1`.

- `P1::Point` and `P2::Point`: extreme `Point` of the integration, placed 
  at comoving distance `s1` and `s2` respectively.

- `y`: the cosine of the angle between the two points `P1` and `P2`

- `cosmo::Cosmology`: cosmology to be used in this computation


See also: [`ξ_lensinglocalgp`](@ref), [`int_on_mu_lensinglocalgp`](@ref)
[`integral_on_mu`](@ref), [`ξ_multipole`](@ref)
"""
function integrand_ξ_lensinglocalgp(
     IP::Point, P1::Point, P2::Point,
     y, cosmo::Cosmology)

     s1 = P1.comdist
     s2, D_s2, a_s2, ℛ_s2 = P2.comdist, P2.D, P2.a, P2.ℛ
     χ1, D1, a1 = IP.comdist, IP.D, IP.a
     Ω_M0 = cosmo.params.Ω_M0

     Δχ1_square = χ1^2 + s2^2 - 2 * χ1 * s2 * y
     Δχ1 = Δχ1_square > 0.0 ? √(Δχ1_square) : 0.0

     common = 9 * ℋ0^4 * Ω_M0^2 * D_s2 * (1 + ℛ_s2) * s2 / (4 * a_s2 * s1)
     factor = D1 * (s1 - χ1) / a1

     new_J31 = -2 * y * Δχ1^2
     new_J22 = χ1 * s2 * (1 - y^2)

     I13 = cosmo.tools.I13(Δχ1)
     I22 = cosmo.tools.I22(Δχ1)

     #println("J00 = $new_J00, \t I00(Δχ1) = $(I00)")
     #println("J02 = $new_J02, \t I20(Δχ1) = $(I20)")
     #println("J31 = $new_J31, \t I13(Δχ1) = $(I13)")
     #println("J22 = $new_J22, \t I22(Δχ1) = $(I22)")

     parenth = (new_J31 * I13 + new_J22 * I22)

     first = common * factor * parenth

     return first
end


function integrand_ξ_lensinglocalgp(
     χ1::Float64, s1::Float64, s2::Float64,
     y, cosmo::Cosmology)

     P1, P2 = Point(s1, cosmo), Point(s2, cosmo)
     IP = Point(χ1, cosmo)
     return integrand_ξ_lensinglocalgp(IP, P1, P2, y, cosmo)
end


@doc raw"""
     ξ_lensinglocalgp(s1, s2, y, cosmo::Cosmology;
          en::Float64 = 1e6, N_χs::Integer = 100):: Float64

Return the Lensing-LocalGP cross-correlation function 
``\xi^{\kappa \phi} (s_1, s_2, \cos{\theta})``, defined as follows:
    
```math
\xi^{\kappa \phi} (s_1, s_2, \cos{\theta}) = 
     \frac{
          9 \mathcal{H}_0^4 \Omega_{M0}^2 D(s_2) (1 + \mathcal{R}(s_2)) s_2
     }{4 a(s_2) s_1} 
     \int_0^{s_1} \mathrm{d}\chi_1 \frac{D(\chi_1)(s_1 - \chi_1) }{a(\chi_1)}
     \left( J_{31} I^3_1(\Delta\chi_1) +  J_{22} I^2_2(\Delta\chi_1) \right)
```

where ``\mathcal{H} = a H``, 
``\Delta\chi_1 = \sqrt{\chi_1^2 + s_2^2 - 2 \chi_1 s_2\cos{\theta}}``, 
``y = \cos{\theta} = \hat{\mathbf{s}}_1 \dot \hat{\mathbf{s}}_2``) 
and the ``J`` coefficients are given by 

```math
\begin{align*}
     J_{31} & = -2 y \Delta\chi_1^2 \\
     J_{22} & = \chi_1 s_2 (1 - y^2)
\end{align*}
```

The computation is made applying [`trapz`](@ref) (see the 
[Trapz](https://github.com/francescoalemanno/Trapz.jl) Julia package) to
the integrand function `integrand_ξ_lensinglocalgp`.


## Inputs

- `s1` and `s2`: comovign distances where the function must be evaluated

- `y`: the cosine of the angle between the two points `P1` and `P2`

- `cosmo::Cosmology`: cosmology to be used in this computation


## Optional arguments 

- `en::Float64 = 1e6`: just a float number used in order to deal better 
  with small numbers;

- `N_χs::Integer = 100`: number of points to be used for sampling the integral
  along the ranges `(0, s1)` (for `χ1`) and `(0, s1)` (for `χ2`); it has been checked that
  with `N_χs ≥ 50` the result is stable.


See also: [`integrand_ξ_lensinglocalgp`](@ref), [`int_on_mu_lensinglocalgp`](@ref)
[`integral_on_mu`](@ref), [`ξ_multipole`](@ref)
"""
function ξ_lensinglocalgp(s1, s2, y, cosmo::Cosmology;
     en::Float64 = 1e6, N_χs::Integer = 100)

     adim_χs = range(1e-6, 1.0, N_χs)
     χ1s = adim_χs .* s1

     P1, P2 = GaPSE.Point(s1, cosmo), GaPSE.Point(s2, cosmo)
     IPs = [GaPSE.Point(x, cosmo) for x in χ1s]

     int_ξs = [
          en * GaPSE.integrand_ξ_lensinglocalgp(IP, P1, P2, y, cosmo)
          for IP in IPs
     ]

     res = trapz(χ1s, int_ξs)
     #println("res = $res")
     return res / en
end



##########################################################################################92



@doc raw"""
     int_on_mu_lensinglocalgp(s1, s, μ, cosmo::Cosmology;
          L::Integer = 0, 
          use_windows::Bool = true, 
          en::Float64 = 1e6,
          N_χs::Integer = 100) :: Float64

Return the integrand on ``\mu = \hat{\mathbf{s}}_1 \dot \hat{\mathbf{s}}`` 
of the Lensing-LocalGP cross-correlation function, i.e.
the following function ``f(s_1, s, \mu)``:

```math
     f(s_1, s, \mu) = \xi^{\kappa \phi} (s_1, s_2, \cos{\theta}) 
          \, \mathcal{L}_L(\mu) \,  \phi(s_2) \, F\left(\frac{s}{s_1}, \mu \right)
```
where ``y =  \cos{\theta} = \hat{\mathbf{s}}_1 \dot \hat{\mathbf{s}}_2`` and
``s = \sqrt{s_1^2 + s_2^2 - 2 \, s_1 \, s_2 \, y}``.

In case `use_windows` is set to `false`, the window functions ``\phi`` and ``F``
are removed, i.e is returned the following function ``f^{'}(s_1, s, \mu)``:

```math
     f^{'}(s_1, s, \mu) = \xi^{\kappa \phi} (s_1, s_2, \cos{\theta}) 
          \, \mathcal{L}_L(\mu) 
```

The function ``\xi^{\kappa \phi}(s_1, s_2, \cos{\theta})`` is calculated
from `ξ_lensing`; note that these is an internal conversion of coordiate sistems
from `(s1, s, μ)` to `(s1, s2, y)` thorugh the functions `y` and `s2`

## Inputs

- `s1`: the comoving distance where must be evaluated the integral

- `s`: the comoving distance from `s1` where must be evaluated the integral

- `μ`: the cosine between `s1` and `s` where must be evaluated the integral

- `cosmo::Cosmology`: cosmology to be used in this computation


## Optional arguments 

- `L::Integer = 0`: order of the Legendre polynomial to be used

- `en::Float64 = 1e6`: just a float number used in order to deal better 
  with small numbers;

- `use_windows::Bool = false`: tells if the integrand must consider the two
   window function ``\phi`` and ``F``

- `N_χs::Integer = 100`: number of points to be used for sampling the integral
  along the ranges `(0, s1)` (for `χ1`) and `(0, s1)` (for `χ2`); it has been checked that
  with `N_χs ≥ 50` the result is stable.

See also: [`integrand_ξ_lensinglocalgp`](@ref), [`ξ_lensinglocalgp`](@ref),
[`integral_on_mu`](@ref), [`map_integral_on_mu`](@ref),
[`spline_F`](@ref), [`ϕ`](@ref), [`Cosmology`](@ref), 
[`y`](@ref), [`s2`](@ref)
"""
function int_on_mu_lensinglocalgp(s1, s, μ, cosmo::Cosmology;
     L::Integer = 0,
     use_windows::Bool = true,
     en::Float64 = 1e6,
     N_χs::Integer = 100)

     s2_value = s2(s1, s, μ)
     y_value = y(s1, s, μ)
     res = if use_windows == true
          ϕ_s2 = ϕ(s2_value; s_min = cosmo.s_min, s_max = cosmo.s_max)
          (ϕ_s2 > 0.0) || (return 0.0)
          #println("s1 = $s1 \t s2 = $(s2(s1, s, μ)) \t  y=$(y(s1, s, μ))")
          int = ξ_lensinglocalgp(s1, s2_value, y_value, cosmo;
               en = en, N_χs = N_χs)
          #println("int = $int")
          int .* (ϕ_s2 * spline_F(s / s1, μ, cosmo.windowF) * Pl(μ, L))
     else
          #println("s1 = $s1 \t s2 = $(s2(s1, s, μ)) \t  y=$(y(s1, s, μ))")
          int = ξ_lensinglocalgp(s1, s2_value, y_value, cosmo;
               en = en, N_χs = N_χs)
          #println("int = $int")
          #println( "Pl(μ, L) = $(Pl(μ, L))")
          int .* Pl(μ, L)
     end

     #println("res = $res")
     return res
end



##########################################################################################92

##########################################################################################92

##########################################################################################92



function ξ_localgplensing(s1, s2, y, cosmo::Cosmology; kwargs...)
     ξ_lensinglocalgp(s2, s1, y, cosmo; kwargs...)
end


function int_on_mu_localgplensing(s1, s, μ, cosmo::Cosmology;
     L::Integer = 0,
     use_windows::Bool = true,
     en::Float64 = 1e6,
     N_χs::Integer = 100)

     s2_value = s2(s1, s, μ)
     y_value = y(s1, s, μ)
     res = if use_windows == true
          ϕ_s2 = ϕ(s2_value; s_min = cosmo.s_min, s_max = cosmo.s_max)
          (ϕ_s2 > 0.0) || (return 0.0)
          int = ξ_localgplensing(s1, s2_value, y_value, cosmo;
               en = en, N_χs = N_χs)
          int .* (ϕ_s2 * spline_F(s / s1, μ, cosmo.windowF) * Pl(μ, L))
     else
          int = ξ_localgplensing(s1, s2_value, y_value, cosmo;
               en = en, N_χs = N_χs)
          int .* Pl(μ, L)
     end

     return res
end





