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
     integrand_ξ_localgpintegratedgp(
          IP::Point, P1::Point, P2::Point,
          y, cosmo::Cosmology) :: Float64

Return the integrand of the LocalGP-IntegratedGP cross-correlation function 
``\xi^{\phi\int\phi} (s_1, s_2, \cos{\theta})``, i.e. the function 
``f(s_1, s_2, y, \chi_1, \chi_2)`` defined as follows:  

```math
f(s_1, s_2, y, \chi_1, \chi_2) = 
     \frac{9 \mathcal{H}_0^4 \Omega_{M0}^2 D(s_1) (\mathcal{R}(s_1) +1)}{2 a(s_1)} 
     \frac{D(\chi_2) \Delta\chi_2^4}{ a(\chi_2)}
     \left(
          \mathcal{H}(\chi_2)( f(\chi_2) - 1) \mathcal{R}(s_2) - \frac{1}{s_2}
     \right) \tilde{I}^4_0(\Delta\chi_2)
```
where ``\mathcal{H} = a H``, 
``\Delta\chi_2 = \sqrt{s_1^2 + \chi_2^2 - 2 s_1 \chi_2 \cos{\theta}}``, 
``y = \cos{\theta} = \hat{\mathbf{s}}_1 \dot \hat{\mathbf{s}}_2``).

## Inputs

- `IP::Point`: `Point` inside the integration limits, placed 
  at comoving distance `χ1`.

- `P1::Point` and `P2::Point`: extreme `Point` of the integration, placed 
  at comoving distance `s1` and `s2` respectively.

- `y`: the cosine of the angle between the two points `P1` and `P2`

- `cosmo::Cosmology`: cosmology to be used in this computation


See also: [`ξ_localgpintegratedgp`](@ref), [`int_on_mu_localgpintegratedgp`](@ref)
[`integral_on_mu`](@ref), [`ξ_multipole`](@ref)
"""
function integrand_ξ_localgpintegratedgp(
     IP::Point, P1::Point, P2::Point,
     y, cosmo::Cosmology)

     s1, a_s1, D_s1, ℛ_s1 = P1.comdist, P1.a, P1.D, P1.ℛ
     s2, ℛ_s2 = P2.comdist, P2.ℛ
     χ2, D2, a2, f2, ℋ2 = IP.comdist, IP.D, IP.a, IP.f, IP.ℋ
     Ω_M0 = cosmo.params.Ω_M0

     Δχ2_square = s1^2 + χ2^2 - 2 * s1 * χ2 * y
     Δχ2 = Δχ2_square > 0 ? √(Δχ2_square) : 0.0

     prefactor = 9 * ℋ0^4 * Ω_M0^2 * D_s1 * (ℛ_s1 + 1) / (2 * a_s1)
     factor = D2 * Δχ2^4 / a2 * (ℋ2 * ℛ_s2 * (f2 - 1) - 1 / s2)
     I04_tilde = cosmo.tools.I04_tilde(Δχ2)

     res = prefactor * factor * I04_tilde

     return res
end


function integrand_ξ_localgpintegratedgp(
     χ2::Float64, s1::Float64, s2::Float64,
     y, cosmo::Cosmology;
     kwargs...)

     P1, P2 = Point(s1, cosmo), Point(s2, cosmo)
     IP = Point(χ2, cosmo)
     return integrand_ξ_localgpintegratedgp(IP, P1, P2, y, cosmo; kwargs...)
end


@doc raw"""
     ξ_localgpintegratedgp(s1, s2, y, cosmo::Cosmology;
          en::Float64 = 1e6, N_χs::Integer = 100):: Float64

Return the LocalGP-IntegratedGP cross-correlation function 
``\xi^{v_{\parallel}\int \phi} (s_1, s_2, \cos{\theta})``, defined as follows:
    
```math
\xi^{v_{\parallel}\int \phi} (s_1, s_2, \cos{\theta}) = 
     \frac{9 \mathcal{H}_0^4 \Omega_{M0}^2 D(s_1) (\mathcal{R}(s_1) +1)}{2 a(s_1)} 
     \int_0^{s_2} \mathrm{d}\chi_2 \frac{D(\chi_2) \Delta\chi_2^4}{ a(\chi_2)}
     \left(
          \mathcal{H}(\chi_2)( f(\chi_2) - 1) \mathcal{R}(s_2) - \frac{1}{s_2}
     \right) \tilde{I}^4_0(\Delta\chi_2)
```
where ``\mathcal{H} = a H``, 
``\Delta\chi_2 = \sqrt{s_1^2 + \chi_2^2 - 2 s_1 \chi_2 \cos{\theta}}``, 
``y = \cos{\theta} = \hat{\mathbf{s}}_1 \dot \hat{\mathbf{s}}_2``).

The computation is made applying [`trapz`](@ref) (see the 
[Trapz](https://github.com/francescoalemanno/Trapz.jl) Julia package) to
the integrand function `integrand_ξ_localgpintegratedgp`.


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


See also: [`integrand_ξ_localgpintegratedgp`](@ref), [`int_on_mu_localgpintegratedgp`](@ref)
[`integral_on_mu`](@ref), [`ξ_multipole`](@ref)
"""
function ξ_localgpintegratedgp(s1, s2, y, cosmo::Cosmology;
     en::Float64 = 1e6, N_χs::Integer = 100)

     #=
     f(χ2) = en * integrand_ξ_localgpintegratedgp(χ2, s1, s2, y, cosmo)

     return quadgk(f, 1e-6, s2; rtol=1e-3)[1] / en
     =#

     adim_χs = range(1e-6, 1.0, N_χs)
     χ2s = adim_χs .* s2

     P1, P2 = GaPSE.Point(s1, cosmo), GaPSE.Point(s2, cosmo)
     IPs = [GaPSE.Point(x, cosmo) for x in χ2s]

     int_ξs = [
          en * GaPSE.integrand_ξ_localgpintegratedgp(IP, P1, P2, y, cosmo)
          for IP in IPs
     ]

     res = trapz(χ2s, int_ξs)
     #println("res = $res")
     return res / en
end



##########################################################################################92



@doc raw"""
     int_on_mu_localgpintegratedgp(s1, s, μ, cosmo::Cosmology;
          L::Integer = 0, 
          use_windows::Bool = true, 
          en::Float64 = 1e6,
          N_χs::Integer = 100) :: Float64

Return the integrand on ``\mu = \hat{\mathbf{s}}_1 \dot \hat{\mathbf{s}}`` 
of the LocalGP-IntegratedGP cross-correlation function, i.e.
the following function ``f(s_1, s, \mu)``:

```math
     f(s_1, s, \mu) = \xi^{\phi\int\phi} (s_1, s_2, \cos{\theta}) 
          \, \mathcal{L}_L(\mu) \,  \phi(s_2) \, F\left(\frac{s}{s_1}, \mu \right)
```
where ``y =  \cos{\theta} = \hat{\mathbf{s}}_1 \dot \hat{\mathbf{s}}_2`` and
``s = \sqrt{s_1^2 + s_2^2 - 2 \, s_1 \, s_2 \, y}``.

In case `use_windows` is set to `false`, the window functions ``\phi`` and ``F``
are removed, i.e is returned the following function ``f^{'}(s_1, s, \mu)``:

```math
     f^{'}(s_1, s, \mu) = \xi^{\phi\int\phi} (s_1, s_2, \cos{\theta}) 
          \, \mathcal{L}_L(\mu) 
```

The function ``\xi^{\phi\int\phi}(s_1, s_2, \cos{\theta})`` is calculated
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

See also: [`integrand_ξ_localgpintegratedgp`](@ref), [`ξ_localgpintegratedgp`](@ref),
[`integral_on_mu`](@ref), [`map_integral_on_mu`](@ref),
[`spline_F`](@ref), [`ϕ`](@ref), [`Cosmology`](@ref), 
[`y`](@ref), [`s2`](@ref)
"""
function int_on_mu_localgpintegratedgp(s1, s, μ, cosmo::Cosmology;
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
          int = ξ_localgpintegratedgp(s1, s2_value, y_value, cosmo;
               en = en, N_χs = N_χs)
          #println("int = $int")
          int .* (ϕ_s2 * spline_F(s / s1, μ, cosmo.windowF) * Pl(μ, L))
     else
          #println("s1 = $s1 \t s2 = $(s2(s1, s, μ)) \t  y=$(y(s1, s, μ))")
          int = ξ_localgpintegratedgp(s1, s2_value, y_value, cosmo;
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



function ξ_integratedgplocalgp(s1, s2, y, cosmo::Cosmology; kwargs...)
    ξ_localgpintegratedgp(s2, s1, y, cosmo; kwargs...)
end


function int_on_mu_integratedgplocalgp(s1, s, μ, cosmo::Cosmology;
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
          int = ξ_integratedgplocalgp(s1, s2_value, y_value, cosmo;
               en = en, N_χs = N_χs)
          #println("int = $int")
          int .* (ϕ_s2 * spline_F(s / s1, μ, cosmo.windowF) * Pl(μ, L))
     else
          #println("s1 = $s1 \t s2 = $(s2(s1, s, μ)) \t  y=$(y(s1, s, μ))")
          int = ξ_integratedgplocalgp(s1, s2_value, y_value, cosmo;
               en = en, N_χs = N_χs)
          #println("int = $int")
          #println( "Pl(μ, L) = $(Pl(μ, L))")
          int .* Pl(μ, L)
     end

     #println("res = $res")
     return res
end
