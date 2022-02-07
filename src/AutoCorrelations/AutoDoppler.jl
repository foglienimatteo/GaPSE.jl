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
     ξ_doppler(P1::Point, P2::Point, y, cosmo::Cosmology; enhancer = 1.0) :: Float64

Return the Doppler auto-correlation function, defined as follows:
```math
\xi^{v_{\parallel}v_{\parallel}} (s_1, s_2, \cos{\theta}) 
= D_1 D_2 f_1 f_2 \mathcal{H}_1 \mathcal{H}_2 \mathcal{R}_1 \mathcal{R}_2 
(J_{00}\, I^0_0(s) + J_{02}\,I^0_2(s) + J_{04}\,I^0_4(s) + J_{20}\,I^2_0(s))
```
where ``D_1 = D(s_1)``, ``D_2 = D(s_2)`` and so on, ``\mathcal{H} = a H`` and 
the J coefficients are given by (with ``y = \cos{\theta}``):
```math
\begin{align*}
    J_{00} (s_1, s_2, y) & = \frac{1}{45} (y^2 s_1 s_2 - 2y(s_1^2 + s_2^2) + 3s_1 s_2) \\
    J_{02} (s_1, s_2, y) & = \frac{2}{63} (y^2 s_1 s_2 - 2y(s_1^2 + s_2^2) + 3s_1 s_2) \\
    J_{04} (s_1, s_2, y) & = \frac{1}{105} (y^2 s_1 s_2 - 2y(s_1^2 + s_2^2) + 3s_1 s_2) \\
    J_{20} (s_1, s_2, y) & = \frac{1}{3} y s^2
\end{align*}
```

## Inputs

- `P1::Point` and `P2::Point`: `Point` where the CF has to be calculated; they contain all the 
  data of interest needed for this calculus (comoving distance, growth factor and so on)
     
- `y`: the cosine of the angle between the two points `P1` and `P2`

- `cosmo::Cosmology`: cosmology to be used in this computation


## Optional arguments

- `enhancer::Float64 = 1.0`: just a float number used in order to deal better with small numbers; the returned
  value is actually `enhancer * xi_doppler`, where `xi_doppler` is the true value calculated
  as shown before.

See also: [`Point`](@ref), [`Cosmology`](@ref)
"""
function ξ_doppler(P1::Point, P2::Point, y, cosmo::Cosmology; enhancer::Float64 = 1.0)
     s1, D1, f1, ℋ1, ℛ1 = P1.comdist, P1.D, P1.f, P1.ℋ, P1.ℛ
     s2, D2, f2, ℋ2, ℛ2 = P2.comdist, P2.D, P2.f, P2.ℋ, P2.ℛ

     Δs = s(P1.comdist, P2.comdist, y)
     prefac = D1 * D2 * f1 * f2 * ℛ1 * ℛ2 * ℋ1 * ℋ2
     c1 = 3 * s1 * s2 - 2 * y * (s1^2 + s2^2) + s1 * s2 * y^2

     I00 = cosmo.tools.I00(Δs)
     I20 = cosmo.tools.I20(Δs)
     I40 = cosmo.tools.I40(Δs)
     I02 = cosmo.tools.I02(Δs)

     parenth = I00 / 45.0 + I20 / 31.5 + I40 / 105.0

     first = prefac * (c1 * parenth + I02 * y * Δs^2 / 3.0)

     return enhancer * first
end



@doc raw"""
     integrand_on_mu_doppler(s1, s, μ, cosmo::Cosmology; 
          L::Integer = 0, 
          enhancer = 1.0,
          use_windows::Bool = true) ::Float64

Return the integrand on ``\mu = \hat{\mathbf{s}}_1 \dot \hat{\mathbf{s}}`` 
of the Doppler auto-correlation function, i.e.
the following function ``f(s_1, s, \mu)``:

```math
     f(s_1, s, \mu) = \xi^{v_{\parallel}v_{\parallel}} (s_1, s_2, \cos{\theta}) 
          \, \mathcal{L}_L(\mu) \,  \phi(s_2) \, F\left(\frac{s}{s_1}, \mu \right)
```
where ``y =  \cos{\theta} = \hat{\mathbf{s}}_1 \dot \hat{\mathbf{s}}_2``

In case `use_windows` is set to `false`, the window function ``\phi`` and ``F``
are removed, i.e is returned the following function ``f^{'}(s_1, s, \mu)``:

```math
     f^{'}(s_1, s, \mu) = \xi^{v_{\parallel}v_{\parallel}} (s_1, s_2, \cos{\theta}) 
          \, \mathcal{L}_L(\mu) 
```

The function ``\xi^{v_{\parallel}v_{\parallel}} (s_1, s_2, \cos{\theta})`` is calculated
from `ξ_doppler`; note that these is an internal conversion of coordiate sistems
from `(s1, s, μ)` to `(s1, s2, y)` thorugh the functions `y` and `s2`


## Inputs

- `s1`: the comoving distance where must be evaluated the integral

- `s`: the comoving distance from `s1` where must be evaluated the integral

- `μ`: the cosine between `s1` and `s` where must be evaluated the integral

- `cosmo::Cosmology`: cosmology to be used in this computation


## Optional arguments

- `L::Integer = 0`: order of the Legendre polynomial to be used

- `enhancer::Float64 = 1.0`: just a float number used in order to deal better with small numbers; the returned
  value is actually `enhancer * f`, where `f` is the true value calculated
  as shown before.

- `use_windows::Bool = false`: tells if the integrand must consider the two
   window function ``\phi`` and ``F``

See also: [`ξ_doppler`](@ref), [`Cosmology`](@ref), [`y`](@ref)
[`s2`](@ref)
"""
function integrand_on_mu_doppler(s1, s, μ,
     cosmo::Cosmology; L::Integer = 0,
     enhancer::Float64 = 1.0,
     use_windows::Bool = true)

     s2_value = s2(s1, s, μ)
     y_value = y(s1, s, μ)

     if use_windows == true
          ϕ_s2 = ϕ(s2_value)
          (ϕ_s2 > 0.0) || (return 0.0)
          P1, P2 = Point(s1, cosmo), Point(s2_value, cosmo)
          val = ξ_doppler(P1, P2, y_value, cosmo; enhancer = enhancer)
          return val * ϕ_s2 * spline_F(s / s1, μ, cosmo.windowF) * Pl(μ, L)
     else
          P1, P2 = Point(s1, cosmo), Point(s2_value, cosmo)
          val = ξ_doppler(P1, P2, y_value, cosmo; enhancer = enhancer)
          return val * Pl(μ, L)
     end
end


##########################################################################################92


#=
function ξ_doppler(s1, s2, y; enhancer = 1.0)
     delta_s = s(s1, s2, y)
     D1, D2 = D(s1), D(s2)

     f1, ℋ1 = f(s1), ℋ(s1)
     f2, ℋ2 = f(s2), ℋ(s2)
     #f1, ℋ1, ℋ1_p, s_b1 = f(s1), ℋ(s1), ℋ_p(s1), s_b(s1)
     #f2, ℋ2, ℋ2_p, s_b2 = f(s2), ℋ(s2), ℋ_p(s2), s_b(s2)
     #ℛ1, ℛ2 = ℛ(s1, ℋ1, ℋ1_p, s_b1), ℛ(s2, ℋ2, ℋ2_p, s_b2)
     ℛ1, ℛ2 = ℛ(s1, ℋ1), ℛ(s2, ℋ2)

     prefac = D1 * D2 * f1 * f2 * ℛ1 * ℛ2 * ℋ1 * ℋ2
     c1 = 3 * s1 * s2 - 2 * y * (s1^2 + s2^2) + s1 * s2 * y^2

     parenth = I00(delta_s) / 45.0 + I20(delta_s) / 31.5 + I40(delta_s) / 105.0

     first = prefac * (c1 * parenth + I02(delta_s) * y * delta_s^2 / 3.0)

     return enhancer * first
end
=#

#=
function integrand_on_mu_doppler(s1, s, μ; L::Integer = 0, enhancer = 1.0,
     use_windows::Bool = true)
     if use_windows == true
          ϕ_s2 = ϕ(s2(s1, s, μ))
          (ϕ_s2 > 0.0) || (return 0.0)
          val = ξ_doppler(s1, s2(s1, s, μ), y(s1, s, μ), enhancer = enhancer)
          return val * ϕ_s2 * spline_F(s / s1, μ) * Pl(μ, L)
     else
          val = ξ_doppler(s1, s2(s1, s, μ), y(s1, s, μ), enhancer = enhancer)
          return val * Pl(μ, L)
     end
end
=#
