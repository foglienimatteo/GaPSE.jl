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
function integrand_ξ_lensing(χ1, χ2, s1, s2, y; tol = 0.5,
     enhancer = 1, Δχ_min = 1)
     (s(s1, s2, y) >= tol) || (return 0.0)

     Δχ = √(χ1^2 + χ2^2 - 2 * χ1 * χ2 * y)
     #println("Δχ = ", Δχ)
     #println("χ1 = ", χ1)
     #println("χ2 = ", χ2)

     D1, D2 = D(χ1), D(χ2)
     a_χ1, a_χ2 = 1 / (1 + z_of_s(χ1)), 1 / (1 + z_of_s(χ2))
     #psb1, psb2 = -1.0, -1.0
     #s_b1, s_b2 = s_b(s1), s_b(s2)
     #psb1, psb2 = 5 * s_b1 - 2, 5 * s_b2 - 2

     denomin = s1 * s2 * a_χ1 * a_χ2
     factor = ℋ0^4 * Ω_M0^2 * D1 * (χ1 - s1) * D2 * (χ2 - s2)
     #println("factor = $factor")
     #println("denomin = $denomin")
     #factor = ℋ0^4 * Ω_M0^2 * D1 * (χ1 - s1) * D2 * (χ2 - s2) * psb1 * psb2

     res = if Δχ > Δχ_min #if abs(y^2 - 1) > y_min #if Δχ > Δχ_min
          χ1χ2 = χ1 * χ2
     
          #abs(y^2 - 1) < y_min ? 0.0 :
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
     
          #println("J00 = $new_J00, \t I00(Δχ) = $(I00(Δχ))")
          #println("J02 = $new_J02, \t I20(Δχ) = $(I20(Δχ))")
          #println("J31 = $new_J31, \t I13(Δχ) = $(I13(Δχ))")
          #println("J22 = $new_J22, \t I22(Δχ) = $(I22(Δχ))")
     
          0.5 * enhancer * factor / denomin / Δχ^4 * (
               new_J00 * I00(Δχ) + new_J02 * I20(Δχ) +
               new_J31 * I13(Δχ) + new_J22 * I22(Δχ)
          )
     else
          lim = 4.0 / 15.0 * (5.0 * σ_2 + 2.0 / 3.0 * σ_0 * χ2^2)
          #println("lim = $lim")
          0.5 * 9.0 / 4.0 * enhancer * factor / denomin * lim
     end

     #println("res = ", res, "\n")
     return res
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
function ξ_lensing(s1, s2, y; tol = 0.5, enhancer = 1, Δχ_min = 1, kwargs...)
      
     my_int(var) = integrand_ξ_lensing(var[1], var[2], s1, s2, y;
          tol = tol, Δχ_min = Δχ_min, enhancer = enhancer)
     a = [0.0, 0.0]
     b = [s1, s2]
     int = hcubature(my_int, a, b; kwargs...)
     #println(int)
     return int
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
function integrand_on_mu_lensing(s1, s, μ; 
               L::Integer = 0, 
               enhancer::Float64 = 1.0,
               tol::Float64 = 0.5,
               Δχ_min::Float64 = 1e-6, 
               χ_atol::Float64 = 1e-3, 
               χ_rtol::Float64 = 1e-3, 
               use_F::Bool=true
          )

     if ϕ(s2(s1, s, μ)) > 0.0
          #println("s1 = $s1 \t s2 = $(s2(s1, s, μ)) \t  y=$(y(s1, s, μ))")
          int = ξ_lensing(s1, s2(s1, s, μ), y(s1, s, μ); Δχ_min = Δχ_min,
               tol = tol, rtol = χ_rtol, atol = χ_atol, enhancer = enhancer)
          #println("int = $int")
          if use_F == true
               return int .* (spline_F(s / s1, μ) * Pl(μ, L))
          else
               return int .* Pl(μ, L)
          end
     else
          return (0.0, 0.0)
     end
end



@doc raw"""
     integral_on_mu_lensing(s1, s;  L::Integer = 0, enhancer = 1e6, 
          Δχ_min = 1e-6, tol = 0.5, χ_atol = 1e-3, χ_rtol = 1e-3, 
          kwargs...) :: Tuple{Float64, Float64}

Evaluate the integral on ``\mu`` of the lensing auto-correlation function 
``\xi^{\kappa\kappa} (s_1, s_2, \cos{\theta})``, weighted with the window functions
``F(x, \mu)`` and `` \phi(s_2)`` and the Legendre polynomial ``\mathcal{L}(\mu)``, 
i.e. the following function ``f_\mathrm{in}(s_1, s)``:

```math
f^{\kappa\kappa}(s_1, s) =  \int_{-1}^{+1} \mathrm{d} \mu \;
        \xi^{\kappa\kappa}(s_1, s, \mu) \, \phi(s_2) \,
        \mathcal{L}_L(\mu) \, F\left(\frac{s}{s_1}, \mu \right)
```

The computation is made applying [`quadgk`](@ref) (see the 
[QuadGK](https://github.com/JuliaMath/QuadGK.jl) Julia package) to
the integrand function `integrand_on_mu_lensing`.


## Optional arguments 

- `L::Integer = 0` : Lagrange polynomial degree to be used in the computation.

- `tol = 0.5` : if ``s \leq \mathrm{tol}``, then the integrand value is `0.0`; 
  it prevents computational problems conserning too
  close points which are insignificants (cosnidering that `tol` is a distance, so it
  is measured in ``h_0^{-1}\,\mathrm{Mpc}``).

- `enhancer = 1e6` : inside the integration, multiply the results in order to deal with
  values not too close to zero; at the end, the result is also divided for this number, 
  restoring the correct value.

- ` Δχ_min = 1e-6` : a Float64 parameter used inside `integrand_ξ_lensing` in order to
  avoid computatinal divergences; it should be `0<Δχ_min<<1`, see the `integrand_ξ_lensing`
  docstring for more informations.

- `χ_atol = 1e-3` : absolute tolerance to be used in the computation of the 2-dims integral
  on ``\chi_1`` and ``\chi_2``, made inside `ξ_lensing`; for computational time reasons,
  it's better to use `χ_atol ≥ 1e-5`

- `χ_rtol = 1e-3` : relative tolerance to be used in the computation of the 2-dims integral
  on ``\chi_1`` and ``\chi_2``, made inside `ξ_lensing`; for computational time reasons,
  it's better to use `χ_rtol ≥ 1e-4`.

- `kwargs...` : keyword arguments which should be passed to `quadgk` when performing
  the 1-dim integral; we shall recomend to set `atol≥1e-4` and `rtol≥-3`, as a consequence
  of the fast increase for the evaluation time needed for this integral.


## Return

A `Tuple{Float64, Float64}` : the former is the integral result, the latter its error.

See also: [`integrand_ξ_lensing`](@ref), [`ξ_lensing`](@ref)
[`integrand_on_mu_lensing`](@ref), [`map_integral_on_mu_lensing`](@ref),
[`spline_F`](@ref), [`ϕ`](@ref)
"""
function integral_on_mu_lensing(s1, s; 
          L::Integer = 0, 
          enhancer::Float64 = 1e6,
          tol::Float64 = 0.5,
          Δχ_min::Float64 = 1e-6, 
          χ_atol::Float64 = 1e-3, 
          χ_rtol::Float64 = 1e-3, 
          use_F::Bool = true,     
          kwargs...
     )

     f(μ) = integrand_on_mu_lensing(s1, s, μ; enhancer = enhancer, tol = tol,
          L = L, Δχ_min = Δχ_min, χ_atol = χ_atol, χ_rtol = χ_rtol, use_F = use_F)[1]
     #println("s1 = $s1 \t s = $s")
     int = quadgk(μ -> f(μ), -1.0, 1.0; kwargs...)
     #println("s1 = $s1 \t s2 = $s \t int = $int")
     return int ./ enhancer
end


function map_integral_on_mu_lensing(
     s1::Float64 = s_eff, v_ss::Union{Vector{Float64}, Nothing}=nothing;
     L::Integer = 0, pr::Bool = true, Δχ_min = 1e-6,
     χ_atol = 1e-3, χ_rtol = 1e-3,
     enhancer = 1e6, tol = 0.5,
     kwargs...)

     t1 = time()
     ss = isnothing(v_ss) ? 10 .^ range(-1, 3, length = 100) : v_ss
     f(s) = integral_on_mu_lensing(s1, s; L = L, enhancer = enhancer,
          Δχ_min = Δχ_min, χ_atol = χ_atol, χ_rtol = χ_rtol, tol = tol, kwargs...)
     vec = @showprogress [f(s) for s in ss]
     xis, xis_err = [x[1] for x in vec], [x[2] for x in vec]
     t2 = time()
     pr && println("\ntime needed for map_integral_on_mu_lensing [in s] = $(t2-t1)")
     return (ss[ss.>tol], xis[ss.>tol], xis_err[ss.>tol])
end


@doc raw"""
     map_integral_on_mu_lensing(s1 = s_eff; L::Integer = 0, pr::Bool = true, 
          Δχ_min = 1e-6, χ_atol = 1e-3, χ_rtol = 1e-3, enhancer = 1e6, 
          tol = 0.5, kwargs...
          ) :: Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}}

Evaluate `integral_on_mu_lensing` in a range of distance values, returning value
and error for each integral.

## Arguments

- `s1 = s_eff` : as already mentioned, in this program we use the effective redshift
  approximation, so the integral on ``s_1`` is replaced with an evaluation of the integrand
  in ``s1 = s_eff``.

## Optional arguments 

- `L::Integer = 0` : Lagrange polynomial degree to be used in the computation.

- `pr::Bool = true` : tells if the println messages should be printed or not.

- `tol = 0.5` : if during the evaluation of the integral inside 
  `integral_on_mu_lensing` happens that 
  ``s =  \sqrt{s_1^2 + s_2^2 - 2 \, s_1 s_2 y} \leq \mathrm{tol}``, 
  then the integrand value is `0.0`; 
  it prevents computational problems conserning too
  close points which are insignificants (cosnidering that `tol` is a distance, so it
  is measured in ``h_0^{-1}\,\mathrm{Mpc}``).

- `enhancer = 1e6` : inside the integration, multiply the results in order to deal with
  values not too close to zero; at the end, the result is also divided for this number, 
  restoring the correct value.

- ` Δχ_min = 1e-6` : a Float64 parameter used inside `integrand_ξ_lensing` in order to
  avoid computatinal divergences; it should be `0<Δχ_min<<1`, see the `integrand_ξ_lensing`
  docstring for more informations.

- `χ_atol = 1e-3` : absolute tolerance to be used in the computation of the 2-dims integral
  on ``\chi_1`` and ``\chi_2``, made inside `ξ_lensing`; for computational time reasons,
  it's better to use `χ_atol ≥ 1e-5`

- `χ_rtol = 1e-3` : relative tolerance to be used in the computation of the 2-dims integral
  on ``\chi_1`` and ``\chi_2``, made inside `ξ_lensing`; for computational time reasons,
  it's better to use `χ_rtol ≥ 1e-4`.

- `kwargs...` : keyword arguments which should be passed to `quadgk` when performing
  the 1-dim integral in `integral_on_mu_lensing`; we shall recomend to set `atol≥1e-4` 
  and `rtol≥-3`, as a consequence of the fast increase for the evaluation time needed for this integral.



## Return

A `Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}}`: 
- the first `Vector{Float64}` contains the distance `s` where the 
  integarl is evaluated;
- the second one the integral results at the corresponding distance; 
- the third and last one contains the error in the integral estimation,
  as returned by `integral_on_mu_lensing`

See also: [`z_eff`](@ref), [`z_eff`](@ref), [`integrand_ξ_lensing`](@ref), 
[`ξ_lensing`](@ref), [`integrand_on_mu_lensing`](@ref), 
[`integral_on_mu_lensing`](@ref), [`print_map_integral_on_mu_lensing`](@ref)
"""
map_integral_on_mu_lensing



#=
function PS_lensing(; int_s_min = 1e-2, int_s_max = 2 * s_max, N = 128,
     L::Integer = 0, pr::Bool = true, kwargs...)

     if ϕ(s_eff) > 0
          t1 = time()
          ks, pks = xicalc(s -> 2 * π^2 * integral_on_mu_lensing(s_eff, s; pr = pr, L = L, kwargs...)[1], L, 0;
               N = N, kmin = int_s_min, kmax = int_s_max, r0 = 1 / int_s_max)
          t2 = time()
          pr && println("\ntime needed for PS_lensing [in s] = $(t2-t1)\n")

          if iseven(L)
               return ks, ((2 * L + 1) / A(s_min, s_max, θ_MAX) * ϕ(s_eff) * (-1)^(L / 2)) .* pks
          else
               return ks, ((2 * L + 1) / A(s_min, s_max, θ_MAX) * ϕ(s_eff) * (-im)^L) .* pks
          end
     else
          throw(ErrorException(
               "ϕ(s_eff) should be >0, but in s_eff=$s_eff " *
               "is ϕ(s_eff) = $(ϕ(s_eff))! \n"
          ))
     end
end
=#



##########################################################################################92



@doc raw"""
     print_map_int_on_mu_lensing(out::String; L::Integer = 0,
          s1 = s_eff, pr::Bool = true, kwargs...)

Print the output of `map_integral_on_mu_lensing` to the given input
filename `out`, with some meta-data at the beginning of it.
If that file already exists, it will be destroyed an re-created.

## Arguments

- `s1 = s_eff` : as already mentioned, in this program we use the effective redshift
  approximation, so the integral on ``s_1`` is replaced with an evaluation of the integrand
  in ``s1 = s_eff``.

## Optional arguments

- `kwargs...` : all the keyword arguments supported by `map_integral_on_mu_lensing`;
  see its docstring for more details.


See also: [`z_eff`](@ref), [`z_eff`](@ref), [`integrand_ξ_lensing`](@ref), 
[`ξ_lensing`](@ref), [`integrand_on_mu_lensing`](@ref), 
[`integral_on_mu_lensing`](@ref), [`map_integral_on_mu_lensing`](@ref)
"""
function print_map_int_on_mu_lensing(out::String,
     v_ss::Union{Vector{Float64},Nothing} = nothing, s1::Float64 = s_eff;
     kwargs...)

     t1 = time()
     vec = map_integral_on_mu_lensing(s1, v_ss; kwargs...)
     t2 = time()

     isfile(out) && run(`rm $out`)
     open(out, "w") do io
          println(io, "# This is an integration map on mu of xi_lensing.")
          parameters_used(io)
          println(io, "# computational time needed (in s) : $(@sprintf("%.4f", t2-t1))")
          print(io, "# kwards passed: ")

          if isempty(kwargs)
               println(io, "none")
          else
               print(io, "\n")
               for (i, key) in enumerate(keys(kwargs))
                    println(io, "# \t\t$(key) = $(kwargs[key])")
               end
          end

          println(io, "\ns \t xi \t xi_error")
          for (s, xi, xi_err) in zip(vec[1], vec[2], vec[3])
               println(io, "$s \t $xi \t $(xi_err)")
          end
     end
end
