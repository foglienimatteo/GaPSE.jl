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
     integrand_on_mu_lensing(s1, s, μ, cosmo::Cosmology;
          L::Integer = 0, 
          use_windows::Bool = true, 
          en::Float64 = 1e6,
          Δχ_min::Float64 = 1e-4,
          N_χs::Integer = 100) :: Float64

Return the integrand on ``\mu = \hat{\mathbf{s}}_1 \dot \hat{\mathbf{s}}`` 
of the lensing auto-correlation function, i.e.
the following function ``f(s_1, s, \mu)``:

```math
     f(s_1, s, \mu) = \xi^{\kappa\kappa} (s_1, s_2, \cos{\theta}) 
          \, \mathcal{L}_L(\mu) \,  \phi(s_2) \, F\left(\frac{s}{s_1}, \mu \right)
```
where ``y =  \cos{\theta} = \hat{\mathbf{s}}_1 \dot \hat{\mathbf{s}}_2`` and
``s = \sqrt{s_1^2 + s_2^2 - 2 \, s_1 \, s_2 \, y}``.

In case `use_windows` is set to `false`, the window functions ``\phi`` and ``F``
are removed, i.e is returned the following function ``f^{'}(s_1, s, \mu)``:

```math
     f^{'}(s_1, s, \mu) = \xi^{\kappa\kappa} (s_1, s_2, \cos{\theta}) 
          \, \mathcal{L}_L(\mu) 
```

The function ``\xi^{\kappa\kappa}(s_1, s_2, \cos{\theta})`` is calculated
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

- ` Δχ_min::Float64 = 1e-4` : parameter used inside `integrand_ξ_lensing` in order to
  avoid computatinal divergences; it should be `0<Δχ_min<<1`, see the `integrand_ξ_lensing`
  docstring for more informations.

- `N_χs::Integer = 100`: number of points to be used for sampling the integral
  along the ranges `(0, s1)` (for `χ1`) and `(0, s1)` (for `χ2`); it has been checked that
  with `N_χs ≥ 50` the result is stable.

See also: [`integrand_ξ_lensing`](@ref), [`ξ_lensing`](@ref),
[`integral_on_mu`](@ref), [`map_integral_on_mu`](@ref),
[`spline_F`](@ref), [`ϕ`](@ref), [`Cosmology`](@ref), 
[`y`](@ref), [`s2`](@ref)
"""
function integrand_on_mu(s1, s, μ, integrand::Function, cosmo::Cosmology;
     L::Integer = 0, use_windows::Bool = true, kwargs...)

     s2_value = s2(s1, s, μ)
     y_value = y(s1, s, μ)
     res = if use_windows == true
          ϕ_s2 = ϕ(s2_value, cosmo.s_min, cosmo.s_max)
          (ϕ_s2 > 0.0) || (return 0.0)
          #println("s1 = $s1 \t s2 = $(s2(s1, s, μ)) \t  y=$(y(s1, s, μ))")
          int = integrand(s1, s2_value, y_value, cosmo; kwargs...)
          #println("int = $int")
          int .* (ϕ_s2 * spline_F(s / s1, μ, cosmo.windowF) * Pl(μ, L))
     else
          #println("s1 = $s1 \t s2 = $(s2(s1, s, μ)) \t  y=$(y(s1, s, μ))")
          int = integrand(s1, s2_value, y_value, cosmo; kwargs...)
          #println("int = $int")
          #println( "Pl(μ, L) = $(Pl(μ, L))")
          int .* Pl(μ, L)
     end

     #println("res = $res")
     return res
end


function integrand_on_mu(s1, s, μ, effect::String, cosmo::Cosmology; kwargs...)
     error = "$effect is not a valid GR effect name.\n" *
             "Valid GR effect names are the following:\n" *
             string(IMPLEMENTED_GR_EFFECTS .* " , "...)
     @assert (effect ∈ IMPLEMENTED_GR_EFFECTS) error

     return integrand_on_mu(s1, s, μ, DICT_GR_ξs[effect], cosmo; kwargs...)
end



##########################################################################################92


function integral_on_mu(
     s1, s, integrand::Function, cosmo::Cosmology;
     L::Integer = 0,
     use_windows::Bool = true,
     enhancer::Float64 = 1e6,
     N_μs::Integer = 50,
     μ_atol::Float64 = 0.0,
     μ_rtol::Float64 = 1e-2,
     SPLINE::Bool = false,
     kwargs...)

     orig_f(μ) = enhancer * integrand_on_mu(s1, s, μ, integrand, cosmo;
          L = L, use_windows = use_windows, kwargs...)

     μs = union(range(-1.0, -0.95, length = N_μs),
               range(-0.95, +0.95, length = N_μs),
               range(+0.95, +1.0, length = N_μs))
     int =
          if s > 1.0 && SPLINE == false
               quadgk(μ -> orig_f(μ), -1.0, 1.0; atol = μ_atol, rtol = μ_rtol)[1]

          elseif s > 1.0 && SPLINE == true
               orig_fs = orig_f.(μs)
               spline_orig_f = Spline1D(μs, orig_fs)
               quadgk(μ -> spline_orig_f(μ), -1.0, 1.0; atol = μ_atol, rtol = μ_rtol)[1]

          else
               orig_fs = orig_f.(μs)
               trapz(μs, orig_fs)
          end

     return int / enhancer
end



function integral_on_mu(s1, s, effect::String, cosmo::Cosmology; kwargs...)
     error = "$effect is not a valid GR effect name.\n" *
             "Valid GR effect names are the following:\n" *
             string(IMPLEMENTED_GR_EFFECTS .* " , "...)

     @assert (effect ∈ IMPLEMENTED_GR_EFFECTS) error
     integral_on_mu(s1, s, DICT_GR_ξs[effect], cosmo; kwargs...)
end



@doc raw"""
     integral_on_mu(s1, s, integrand::Function, cosmo::Cosmology;
          L::Integer = 0,
          enhancer::Float64 = 1e6,
          use_windows::Bool = true,
          μ_atol::Float64 = 1e-4,
          μ_rtol::Float64 = 1e-1,
          kwargs...
          )

     integral_on_mu(s1, s, effect::String, cosmo::Cosmology; kwargs...)

Evaluate the integral on ``\mu`` of the chosen correlation function term, 
through the `quadgk` function (see the [QuadGK](https://github.com/JuliaMath/QuadGK.jl) 
Julia package).

In the former method you have to pass as an input the `integrand` function you want 
to integrate, while in the (recommended) latter one it's necessary to specify the
name of the CF term among the following 
`$(string(IMPLEMENTED_GR_EFFECTS .* " , "...))`
to which correspond the following functions:
`$(string(string.(IMPLEMENTED_INTEGRANDS) .* " , "...))`

The integral evaluated is then the following:

```math
     f(s_1, s, \mu) = \int_{-1}^{+1} \mathrm{d}\mu \; \xi (s_1, s_2, \cos{\theta}) 
          \, \mathcal{L}_L(\mu) \,  \phi(s_2) \, F\left(\frac{s}{s_1}, \mu \right)
```
for `use_windows==true` and 

```math
     f^{'}(s_1, s, \mu) = \int_{-1}^{+1} \mathrm{d}\mu \;  
          \xi^{v_{\parallel}v_{\parallel}} (s_1, s_2, \cos{\theta}) 
          \, \mathcal{L}_L(\mu) 
```
for `use_windows==false`, where ``y =  \cos{\theta} = \hat{\mathbf{s}}_1 \dot \hat{\mathbf{s}}_2``
and ``\xi`` is the chosen CF effect. 

## Inputs

- `s1`: the comoving distance where must be evaluated the integral

- `s`: the comoving distance from `s1` where must be evaluated the integral

- `cosmo::Cosmology`: cosmology to be used in this computation


## Optional arguments

- `L::Integer = 0`: order of the Legendre polynomial to be used

- `enhancer::Float64 = 1e6`: just a float number used in order to deal better with small numbers; 
  the returned value is NOT modified by this value, because after a multiplication
  the internal result is divided by `enhancer`.

- `use_windows::Bool = false`: tells if the integrand must consider the two
   window function ``\phi`` and ``F``

- `μ_atol::Float64 = 1e-3` and `μ_rtol::Float64 = 1e-3`: absolute and relative tolerance
  to be passed to `quadgk`; it's recommended not to set `μ_rtol < 1e-2` because
  of long time for evaluations

"""
integral_on_mu



##########################################################################################92



function ξ_multipole(s1, s, effect::Function, cosmo::Cosmology; L::Integer = 0, kwargs...)
     error = "$(string(effect)) is not a valid GR effect function.\n" *
             "Valid GR effect functions are the following:\n" *
             string(string.(IMPLEMENTED_ξs) .* " , "...)
     @assert (effect ∈ IMPLEMENTED_ξs) error

     return (2.0 * L + 1.0) / 2.0 * integral_on_mu(s1, s, effect, cosmo; L = L, kwargs...)
end


function ξ_multipole(s1, s, effect::String, cosmo::Cosmology; L::Integer = 0, kwargs...)
     error = "$effect is not a valid GR effect name.\n" *
             "Valid GR effect names are the following:\n" *
             string(IMPLEMENTED_GR_EFFECTS .* " , "...)
     @assert (effect ∈ IMPLEMENTED_GR_EFFECTS) error

     return (2L + 1) / 2 * integral_on_mu(s1, s, DICT_GR_ξs[effect], cosmo; L = L, kwargs...)
end


##########################################################################################92


function map_integral_on_mu(
     cosmo::Cosmology,
     effect::Union{String,Function},
     v_ss::Union{Vector{Float64},Nothing} = nothing;
     s_1::Union{Float64,Nothing} = nothing,
     pr::Bool = true, 
     N_log::Integer = 1000,
     L::Integer = 0,
     kwargs...)

     s1 = isnothing(s_1) ? cosmo.s_eff : s_1

     t1 = time()
     ss = isnothing(v_ss) ? 10 .^ range(-1, 3, length =  N_log) : v_ss
     xis = pr ? begin
          @showprogress "$effect, L=$L: " [
          integral_on_mu(s1, s, effect, cosmo; L = L, kwargs...) for s in ss
          ] end : [
          integral_on_mu(s1, s, effect, cosmo; L = L, kwargs...) for s in ss
          ]

     t2 = time()
     pr && println("\ntime needed for map_integral_on_mu for $effect " *
                   "[in s] = $(@sprintf("%.5f", t2-t1)) ")
     return (ss, xis)
end


function map_ξ_multipole(
     cosmo::Cosmology,
     effect::Union{String,Function},
     v_ss::Union{Vector{Float64},Nothing} = nothing;
     L::Integer = 0, kwargs...)

     ss, xis = map_integral_on_mu(cosmo, effect, v_ss; L = L, kwargs...)

     return (ss, (2.0 * L + 1.0) / 2.0 * xis)
end


##########################################################################################92


function print_map_int_on_mu(
     cosmo::Cosmology,
     out::String,
     effect::Union{String,Function},
     v_ss::Union{Vector{Float64},Nothing} = nothing;
     s_1::Union{Float64,Nothing} = nothing,
     kwargs...)

     s1 = isnothing(s_1) ? cosmo.s_eff : s_1
     t1 = time()
     vec = map_integral_on_mu(cosmo, effect, v_ss; s_1 = s1, kwargs...)
     t2 = time()

     isfile(out) && run(`rm $out`)
     open(out, "w") do io
          println(io, "# This is an integration map on mu of the $effect GR effect.")
          parameters_used(io, cosmo)
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
          println(io, "# ")
          println(io, "# s [Mpc/h_0] \t \t xi")
          for (s, xi) in zip(vec[1], vec[2])
               println(io, "$s \t $xi")
          end
     end
end

function print_map_ξ_multipole(
     cosmo::Cosmology,
     out::String,
     effect::Union{String,Function},
     v_ss::Union{Vector{Float64},Nothing} = nothing;
     s_1::Union{Float64,Nothing} = nothing,
     kwargs...)

     s1 = isnothing(s_1) ? cosmo.s_eff : s_1
     t1 = time()
     vec = map_ξ_multipole(cosmo, effect, v_ss; s_1 = s1, kwargs...)
     t2 = time()

     isfile(out) && run(`rm $out`)
     open(out, "w") do io
          println(io, "# This is an integration map on mu of the ξ multipole $effect GR effect.")
          parameters_used(io, cosmo)
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
          println(io, "# ")
          println(io, "# s [Mpc/h_0] \t \t xi")
          for (s, xi) in zip(vec[1], vec[2])
               println(io, "$s \t $xi")
          end
     end
end

