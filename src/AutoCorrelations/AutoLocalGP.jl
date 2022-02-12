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
     ξ_localGP(P1::Point, P2::Point, y, cosmo::Cosmology) :: Float64

Return the local gravitational potential auto-correlation function, 
defined as follows:
```math
\xi^{\phi\phi} (s_1, s_2, \cos{\theta}) = 
     \frac{9 \mathcal{H}_0^4 \Omega_{M0}^2 D(s_1) D(s_2)s^4}{4 a(s_1) a(s_2)}
     (1 + \mathcal{R}_1 + \mathcal{R}_2 + \mathcal{R}_1\mathcal{R}_2)
     \tilde{I}^4_0(s)
```
where ``D_1 = D(s_1)``, ``D_2 = D(s_2)`` and so on, ``\mathcal{H} = a H``, 
``y = \cos{\theta} = \hat{\mathbf{s}}_1 \dot \hat{\mathbf{s}}_2`` and:
```math
\tilde{I}^4_0 (s) &= \int_0^\infty \frac{\mathrm{d}q}{2\pi^2} 
          q^2 \, P(q) \, \frac{j_0(q s) - 1}{(q s)^4}
```

## Inputs

- `P1::Point` and `P2::Point`: `Point` where the CF has to be calculated; they contain all the 
  data of interest needed for this calculus (comoving distance, growth factor and so on)
     
- `y`: the cosine of the angle between the two points `P1` and `P2`

- `cosmo::Cosmology`: cosmology to be used in this computation

See also: [`Point`](@ref), [`Cosmology`](@ref)
"""
function ξ_localGP(P1::Point, P2::Point, y, cosmo::Cosmology)
     s1, D1, f1, a1, ℛ1 = P1.comdist, P1.D, P1.f, P1.a, P1.ℛ
     s2, D2, f2, a2, ℛ2 = P2.comdist, P2.D, P2.f, P2.a, P2.ℛ

     Δs = s(s1, s2, y)
     prefac = 2.25 * ℋ0^4 * cosmo.params.Ω_M0^2 * D1 * D2 * Δs^4 / (a1 * a2)
     parenth = 1.0 + ℛ1 + ℛ2 + ℛ1 * ℛ2

     I04_tilde = cosmo.tools.I04_tilde(Δs)

     res = prefac * I04_tilde * parenth

     return res
end


function ξ_localGP(s1, s2, y, cosmo::Cosmology)
     P1, P2 = Point(s1, cosmo), Point(s2, cosmo)
     return ξ_localGP(P1, P2, y, cosmo)
end


##########################################################################################92



@doc raw"""
     integrand_on_mu_localGP(s1, s, μ, cosmo::Cosmology;
          L::Integer = 0, 
          use_windows::Bool = true) :: Float64

Return the integrand on ``\mu = \hat{\mathbf{s}}_1 \dot \hat{\mathbf{s}}`` 
of the local gravitational potential auto-correlation function, i.e.
the following function ``f(s_1, s, \mu)``:

```math
     f(s_1, s, \mu) = \xi^{\phi\phi} (s_1, s_2, \cos{\theta}) 
          \, \mathcal{L}_L(\mu) \,  \phi(s_2) \, F\left(\frac{s}{s_1}, \mu \right)
```
where ``y =  \cos{\theta} = \hat{\mathbf{s}}_1 \dot \hat{\mathbf{s}}_2`` and
``s = \sqrt{s_1^2 + s_2^2 - 2 \, s_1 \, s_2 \, y}``.

In case `use_windows` is set to `false`, the window functions ``\phi`` and ``F``
are removed, i.e is returned the following function ``f^{'}(s_1, s, \mu)``:

```math
     f^{'}(s_1, s, \mu) = \xi^{\phi\phi} (s_1, s_2, \cos{\theta}) 
          \, \mathcal{L}_L(\mu) 
```

The function ``\xi^{\phi\phi}(s_1, s_2, \cos{\theta})`` is calculated
from `ξ_localGP`; note that these is an internal conversion of coordiate sistems
from `(s1, s, μ)` to `(s1, s2, y)` thorugh the functions `y` and `s2`

## Inputs

- `s1`: the comoving distance where must be evaluated the integral

- `s`: the comoving distance from `s1` where must be evaluated the integral

- `μ`: the cosine between `s1` and `s` where must be evaluated the integral

- `cosmo::Cosmology`: cosmology to be used in this computation


## Optional arguments 

- `L::Integer = 0`: order of the Legendre polynomial to be used

- `use_windows::Bool = false`: tells if the integrand must consider the two
   window function ``\phi`` and ``F``.


See also: [`ξ_localGP`](@ref),
[`integral_on_mu`](@ref), [`map_integral_on_mu`](@ref),
[`spline_F`](@ref), [`ϕ`](@ref), [`Cosmology`](@ref), 
[`y`](@ref), [`s2`](@ref)
"""
function integrand_on_mu_localGP(s1, s, μ,
     cosmo::Cosmology; L::Integer = 0,
     use_windows::Bool = true)

     s2_value = s2(s1, s, μ)
     y_value = y(s1, s, μ)

     if use_windows == true
          ϕ_s2 = ϕ(s2_value)
          (ϕ_s2 > 0.0) || (return 0.0)
          val = ξ_localGP(s1, s2_value, y_value, cosmo)
          return val * ϕ_s2 * spline_F(s / s1, μ, cosmo.windowF) * Pl(μ, L)
     else
          val = ξ_localGP(s1, s2_value, y_value, cosmo)
          return val * Pl(μ, L)
     end
end

