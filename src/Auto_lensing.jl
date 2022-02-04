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
     integrand_ξ_lensing(χ1, χ2, s1, s2, y; tol = 0.5, 
          enhancer = 1, Δχ_min = 1e-6) :: Float64

Return the integarnd of the lensing auto-correlation function 
``\xi^{\kappa\kappa} (s_1, s_2, \cos{\theta})``, i.e. the function 
``f(s_1, s_2, y, \chi_1, \chi_2)`` defined as follows:  

    
```math
f(s_1, s_2, y, \chi_1, \chi_2) = 
\frac{1}{2}
\frac{
     \mathcal{H}_0^4 \Omega_{ \mathrm{M0}}^2 D_1 D_2 (\chi_1 - s_1)(\chi_2 - s_2)
}{
     s_1 s_2 a(\chi_1) a(\chi_2) }
(J_{00} \, I^0_0(\chi) + J_{02} \, I^0_2(\chi) + 
     J_{31} \, I^3_1(\chi) + J_{22} \, I^2_2(\chi))
```

where ``D_1 = D(\chi_1)``, ``D_2 = D(\chi_2)`` and so on, ``\mathcal{H} = a H``, 
``\chi = \sqrt{\chi_1^2 + \chi_2^2 - 2\chi_1\chi_2\cos{\theta}}`` 
and the ``J`` coefficients are given by (with ``y = \cos{\theta}``)

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

## Optional arguments 

- `tol = 0.5` : if ``s = \sqrt{s_1^2 + s_2^2 - 2 \, s_1 s_2 y} \leq \mathrm{tol}``,
  then the returned value is `0.0`; it prevents computational problems conserning too
  close points which are insignificants (cosnidering that `tol` is a distance, so it
  is measured in ``h_0^{-1}\,\mathrm{Mpc}``).

- `enhancer = 1` : multiply the resulting ``f(s_1, s_2, y, \chi_1, \chi_2)`` value; it
  is very useful for interal computations in other functions (for instance 
  `map_integral_on_mu_lensing`), in order to deal better with small float numbers.

- ` Δχ_min = 1e-6` : when ``\Delta\chi = \sqrt{\chi_1^2 + \chi_2^2 - 2 \, \chi_1 \chi_2 y} \to 0^{+}``,
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
[`integral_on_mu_lensing`](@ref), [`map_integral_on_mu_lensing`](@ref)
"""
function integrand_ξ_lensing(
     IP1::Point, IP2::Point,
     P1::Point, P2::Point,
     y, cosmo::Cosmology;
     enhancer = 1, Δχ_min = 1)

     s1 = P1.comdist
     s2 = P2.comdist
     χ1, D1, a_χ1 = IP1.comdist, IP1.D, IP1.a
     χ2, D2, a_χ2 = IP2.comdist, IP2.D, IP2.a
     Ω_M0 = cosmo.params.Ω_M0

     Δχ = √(χ1^2 + χ2^2 - 2 * χ1 * χ2 * y)

     denomin = s1 * s2 * a_χ1 * a_χ2
     factor = ℋ0^4 * Ω_M0^2 * D1 * (χ1 - s1) * D2 * (χ2 - s2)
     #println("factor = $factor")
     #println("denomin = $denomin")
     #factor = ℋ0^4 * Ω_M0^2 * D1 * (χ1 - s1) * D2 * (χ2 - s2) * psb1 * psb2

     res = if Δχ > Δχ_min
          χ1χ2 = χ1 * χ2

          new_J00 = -0.75 * χ1χ2^2 * (y^2 - 1) * (8 * y * (χ1^2 + χ2^2) - χ1χ2 * (9 * y^2 + 7))
          new_J02 = -1.5 * χ1χ2^2 * (y^2 - 1) * (4 * y * (χ1^2 + χ2^2) - χ1χ2 * (3 * y^2 + 5))
          new_J31 = 9 * y * Δχ^6
          new_J22 = 2.25 * χ1χ2 * (
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

          enhancer * factor / denomin / Δχ^4 * (
               new_J00 * I00 + new_J02 * I20 +
               new_J31 * I13 + new_J22 * I22
          )
     else

          lim = 4.0 / 15.0 * (5.0 * cosmo.tools.σ_2 + 2.0 / 3.0 * cosmo.tools.σ_0 * χ2^2)
          #println("lim = $lim")
          9.0 / 4.0 * enhancer * factor / denomin * lim
     end

     #println("res = ", res, "\n")
     return res
end


function func_Δχ_min(s1, s2, y; frac = 1e-4)
     Δs = s(s1, s2, y)
     return frac * Δs
end


@doc raw"""
     ξ_lensing(s1, s2, y; tol = 0.5, 
          enhancer = 1, Δχ_min = 1e-6, kwargs...) :: Tuple{Float64, Float64}

Return the Lensing Auto-correlation function 
``\xi^{\kappa\kappa} (s_1, s_2, \cos{\theta})`` , defined as follows:
    
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
``\chi = \sqrt{\chi_1^2 + \chi_2^2 - 2\chi_1\chi_2\cos{\theta}}`` 
and the ``J`` coefficients are given by (with ``y = \cos{\theta}``)

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

The computation is made applying [`hcubature`](@ref) (see the 
[Hcubature](https://github.com/JuliaMath/HCubature.jl) Julia package) to
the integrand function `integrand_ξ_lensing`.


## Optional arguments 

- `tol = 0.5` : if ``s = \sqrt{s_1^2 + s_2^2 - 2 \, s_1 s_2 y} \leq \mathrm{tol}``,
  then the returned value is `0.0`; it prevents computational problems conserning too
  close points which are insignificants (cosnidering that `tol` is a distance, so it
  is measured in ``h_0^{-1}\,\mathrm{Mpc}``).

- `enhancer = 1` : multiply the resulting value; it
  is very useful for interal computations in other functions (for instance 
  `map_integral_on_mu_lensing`), in order to deal better with small float numbers.

- ` Δχ_min = 1e-6` : when ``\Delta\chi = \sqrt{\chi_1^2 + \chi_2^2 - 2 \, \chi_1 \chi_2 y} \to 0^{+}``,
  some ``I_\ell^n`` term diverges, but the overall parenthesis has a known limit:

  ```math
     \lim_{\chi\to0^{+}} (J_{00} \, I^0_0(\chi) + J_{02} \, I^0_2(\chi) + 
          J_{31} \, I^3_1(\chi) + J_{22} \, I^2_2(\chi)) = 
          \frac{4}{15} \, (5 \, \sigma_2 + \frac{2}{3} \, σ_0 \,s_1^2 \, \chi_2^2)
  ```

  So, when it happens that ``\chi < \Delta\chi_\mathrm{min}``, the function considers this limit
  as the result of the parenthesis instead of calculating it in the normal way; it prevents
  computational divergences.

- `kwargs...` : keyword arguments which should be passed to `hcubature` when performing
  the 2-dims integral; we shall recomend to set `atol≥1e-4` and `rtol≥-3`, as a consequence
  of the fast increase for the evaluation time needed for this integral.


## Return

A `Tuple{Float64, Float64}` : the former is the integral result, the latter its error.

See also: [`integrand_ξ_lensing`](@ref), [`integrand_on_mu_lensing`](@ref)
[`integral_on_mu_lensing`](@ref), [`map_integral_on_mu_lensing`](@ref)
"""
function ξ_lensing(s1, s2, y, cosmo::GaPSE.Cosmology;
     N_χs = 100, enhancer = 1, frac_Δχ_min = 1e-3)

     adim_χs = range(1e-6, 1.0, N_χs)
     Δχ_min = func_Δχ_min(s1, s2, y; frac = frac_Δχ_min)

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



##########################################################################################92



@doc raw"""
     integrand_on_mu_lensing(s1, s, μ; L::Integer = 0, enhancer = 1,
          Δχ_min = 1e-6, tol = 0.5, 
          χ_atol = 1e-3, χ_rtol = 1e-3) :: Tuple{Float64, Float64}

Evaluate the following function ``f^{\kappa\kappa}(s_1, s, \mu)``:

```math
f^{\kappa\kappa}(s_1, s, \mu) = \xi^{\kappa\kappa}(s_1, s, \mu) \, \phi(s_2) \,
        \mathcal{L}_L(\mu) \, F\left(\frac{s}{s_1}, \mu \right)
```

The evaluation of ``\xi^{\kappa\kappa}(s_1, s, \mu)`` is made through `ξ_lensing`.

## Optional arguments 

- `tol = 0.5` : if ``s \leq \mathrm{tol}``, then the integrand value is `0.0`; 
  it prevents computational problems conserning too
  close points which are insignificants (cosnidering that `tol` is a distance, so it
  is measured in ``h_0^{-1}\,\mathrm{Mpc}``).

- `enhancer = 1` : multiply the resulting ``f(s_1, s_2, y, \chi_1, \chi_2)`` value; it
  is very useful for interal computations in other functions (for instance 
  `map_integral_on_mu_lensing`), in order to deal better with small float numbers.

- ` Δχ_min = 1e-6` : a Float64 parameter used inside `integrand_ξ_lensing` in order to
  avoid computatinal divergences; it should be `0<Δχ_min<<1`, see the `integrand_ξ_lensing`
  docstring for more informations.

- `χ_atol = 1e-3` : absolute tolerance to be used in the computation of the 2-dims integral
  on ``\chi_1`` and ``\chi_2``, made inside `ξ_lensing`; for computational time reasons,
  it's better to use `χ_atol ≥ 1e-5`

- `χ_rtol = 1e-3` : relative tolerance to be used in the computation of the 2-dims integral
  on ``\chi_1`` and ``\chi_2``, made inside `ξ_lensing`; for computational time reasons,
  it's better to use `χ_rtol ≥ 1e-4`.


## Returns

A `Tuple{Float64, Float64}` : the former is the integral ressult, the latter its error.


See also: [`integrand_ξ_lensing`](@ref), [`ξ_lensing`](@ref)
[`integral_on_mu_lensing`](@ref), [`map_integral_on_mu_lensing`](@ref),
[`spline_F`](@ref), [`ϕ`](@ref)
"""
function integrand_on_mu_lensing(s1, s, μ, cosmo::Cosmology;
     L::Integer = 0, use_windows::Bool = true, kwargs...)

     s2_value = s2(s1, s, μ)
     y_value = y(s1, s, μ)
     res = if use_windows == true
          ϕ_s2 = ϕ(s2_value; s_min = cosmo.s_min, s_max = cosmo.s_max)
          (ϕ_s2 > 0.0) || (return 0.0)
          #println("s1 = $s1 \t s2 = $(s2(s1, s, μ)) \t  y=$(y(s1, s, μ))")
          int = ξ_lensing(s1, s2_value, y_value, cosmo; kwargs...)
          #println("int = $int")
          int .* (ϕ_s2 * spline_F(s / s1, μ, cosmo.windowF) * Pl(μ, L))
     else
          #println("s1 = $s1 \t s2 = $(s2(s1, s, μ)) \t  y=$(y(s1, s, μ))")
          int = ξ_lensing(s1, s2_value, y_value, cosmo; kwargs...)
          #println("int = $int")
          int .* Pl(μ, L)
     end

     #println("res = $res")
     return res
end





