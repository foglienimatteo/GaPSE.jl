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
          y, cosmo::Cosmology) :: Float64

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


See also: [`ξ_integratedGP`](@ref), [`integrand_on_mu_integratedGP`](@ref)
[`integral_on_mu`](@ref), [`ξ_multipole`](@ref)
"""
function integrand_ξ_integratedGP(IP1::Point, IP2::Point,
     P1::Point, P2::Point,
     y, cosmo::Cosmology)

     s1, ℛ_s1 = P1.comdist, P1.ℛ
     s2, ℛ_s2 = P2.comdist, P2.ℛ
     χ1, D1, a_χ1, ℋ1, f1 = IP1.comdist, IP1.D, IP1.a, IP1.ℋ, IP1.f
     χ2, D2, a_χ2, ℋ2, f2 = IP2.comdist, IP2.D, IP2.a, IP2.ℋ, IP2.f
     Ω_M0 = cosmo.params.Ω_M0

     Δχ = √(χ1^2 + χ2^2 - 2 * χ1 * χ2 * y)

     factor = 9 * ℋ0^4 * Ω_M0^2 * D1 * D2 * Δχ^4 / (s1 * s2 * a_χ1 * a_χ2)
     par_1 = s1 * ℋ1 * ℛ_s1 * (f1 - 1) - 1
     par_2 = s2 * ℋ2 * ℛ_s2 * (f2 - 1) - 1
     #println("factor = $factor")
     #println("denomin = $denomin")

     #println(Δχ)
     I04_t = cosmo.tools.I04_tilde(Δχ)

     return factor * par_1 * par_2 * I04_t
end



@doc raw"""
     ξ_integratedGP(s1, s2, y, cosmo::Cosmology; 
          en::Float64 = 1e10,
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

- `en::Float64 = 1e10`: just a float number used in order to deal better 
  with small numbers.

- `N_χs::Integer = 100`: number of points to be used for sampling the integral
  along the ranges `(0, s1)` (for `χ1`) and `(0, s1)` (for `χ2`); it has been checked that
  with `N_χs ≥ 50` the result is stable.

See also: [`integrand_ξ_integratedGP`](@ref), [`integrand_on_mu_integratedGP`](@ref)
[`integral_on_mu`](@ref), [`ξ_multipole`](@ref)
"""
function ξ_integratedGP(P1::Point, P2::Point, y, cosmo::Cosmology;
     en::Float64 = 1e10, N_χs::Integer = 100, focus::Float64 = 10.0)

     #adim_χs = range(1e-12, 1.0, N_χs)
     adim_χs = range(0.0, 1.0, length = N_χs)[begin+1:end]
     #Δχ_min = func_Δχ_min(s1, s2, y; frac = frac_Δχ_min)

     χ1s = adim_χs .* P1.comdist
     χ2s = adim_χs .* P2.comdist

     IP1s = [GaPSE.Point(x, cosmo) for x in χ1s]
     IP2s = [GaPSE.Point(x, cosmo) for x in χ2s]

     int_ξ_igp = [
          en * GaPSE.integrand_ξ_integratedGP(IP1, IP2, P1, P2, y, cosmo)
          for IP1 in IP1s, IP2 in IP2s
     ]

     res = trapz((χ1s, χ2s), int_ξ_igp)
     #println("res = $res")

     #=
     χ1s = [x for x in range(0.0, P1.comdist, length = N_χs)[begin+1:end]]
     l = Int(floor(N_χs/2))
     matrix_χ2s = [begin
          a = [x for x in range(0.0, P2.comdist, length=l)[begin+1:end]];
          b = [x for x in range(x1-focus, x1+focus, length=l)];
          vcat(a[a.<x1-focus], b, a[a.>x1+focus])
          end for x1 in χ1s]

     IP1s = [GaPSE.Point(x, cosmo) for x in χ1s]
     matrix_IP2s = [[GaPSE.Point(x, cosmo) for x in y] for y in matrix_χ2s]
     matrix_int_ξs = [
          [en * GaPSE.integrand_ξ_integratedGP(IP1, IP2, P1, P2, y, cosmo) 
          for IP2 in matrix_IP2s[i]]
          for (i,IP1) in enumerate(IP1s)]
     
     vec_trapz = [trapz(χ2s,int_ξs) for (χ2s,int_ξs) in zip(matrix_χ2s, matrix_int_ξs)]
     res = trapz(χ1s, vec_trapz)
     =#
     return res / en
end


function ξ_integratedGP(s1, s2, y, cosmo::Cosmology; kwargs...)
     P1, P2 = Point(s1, cosmo), Point(s2, cosmo)
     return ξ_integratedGP(P1, P2, y, cosmo; kwargs...)
end


##########################################################################################92



@doc raw"""
     integrand_on_mu_integratedGP(s1, s, μ, cosmo::Cosmology;
          L::Integer = 0, 
          use_windows::Bool = true,
          en::Float64 = 1e10,
          N_χs::Integer = 100) :: Float64

Return the integrand on ``\mu = \hat{\mathbf{s}}_1 \dot \hat{\mathbf{s}}`` 
of the integrated gravitational potential auto-correlation function, i.e.
the following function ``f(s_1, s, \mu)``:

```math
     f(s_1, s, \mu) = \xi^{\int\phi\int\phi} (s_1, s_2, \cos{\theta}) 
          \, \mathcal{L}_L(\mu) \,  \phi(s_2) \, F\left(\frac{s}{s_1}, \mu \right)
```
where ``y =  \cos{\theta} = \hat{\mathbf{s}}_1 \dot \hat{\mathbf{s}}_2`` and
``s = \sqrt{s_1^2 + s_2^2 - 2 \, s_1 \, s_2 \, y}``.

In case `use_windows` is set to `false`, the window functions ``\phi`` and ``F``
are removed, i.e is returned the following function ``f^{'}(s_1, s, \mu)``:

```math
     f^{'}(s_1, s, \mu) = \xi^{\int\phi\int\phi} (s_1, s_2, \cos{\theta}) 
          \, \mathcal{L}_L(\mu) 
```

The function ``\xi^{\int\phi\int\phi}(s_1, s_2, \cos{\theta})`` is calculated
from `ξ_integratedGP`; note that these is an internal conversion of coordiate sistems
from `(s1, s, μ)` to `(s1, s2, y)` thorugh the functions `y` and `s2`

## Inputs

- `s1`: the comoving distance where must be evaluated the integral

- `s`: the comoving distance from `s1` where must be evaluated the integral

- `μ`: the cosine between `s1` and `s` where must be evaluated the integral

- `cosmo::Cosmology`: cosmology to be used in this computation


## Optional arguments 

- `L::Integer = 0`: order of the Legendre polynomial to be used

- `en::Float64 = 1e10`: just a float number used in order to deal better 
  with small numbers;

- `use_windows::Bool = false`: tells if the integrand must consider the two
   window function ``\phi`` and ``F``.

- `N_χs::Integer = 100`: number of points to be used for sampling the integral
  along the ranges `(0, s1)` (for `χ1`) and `(0, s1)` (for `χ2`); it has been checked that
  with `N_χs ≥ 50` the result is stable.

See also: [`integrand_ξ_integratedGP`](@ref), [`ξ_integratedGP`](@ref),
[`integral_on_mu`](@ref), [`map_integral_on_mu`](@ref),
[`spline_F`](@ref), [`ϕ`](@ref), [`Cosmology`](@ref), 
[`y`](@ref), [`s2`](@ref)
"""
function integrand_on_mu_integratedGP(s1, s, μ, cosmo::Cosmology;
     L::Integer = 0, en::Float64 = 1e10,
     use_windows::Bool = true,
     N_χs::Integer = 50)

     s2_value = s2(s1, s, μ)
     y_value = y(s1, s, μ)
     res = if use_windows == true
          ϕ_s2 = ϕ(s2_value; s_min = cosmo.s_min, s_max = cosmo.s_max)
          (ϕ_s2 > 0.0) || (return 0.0)
          #println("s1 = $s1 \t s2 = $(s2(s1, s, μ)) \t  y=$(y(s1, s, μ))")
          int = ξ_integratedGP(s1, s2_value, y_value, cosmo; en = en, N_χs = N_χs)
          #println("int = $int")
          int .* (ϕ_s2 * spline_F(s / s1, μ, cosmo.windowF) * Pl(μ, L))
     else
          #println("s1 = $s1 \t s2 = $(s2(s1, s, μ)) \t  y=$(y(s1, s, μ))")
          int = ξ_integratedGP(s1, s2_value, y_value, cosmo; en = en, N_χs = N_χs)
          #println("int = $int")
          int .* Pl(μ, L)
     end

     #println("res = $res")
     return res
end

