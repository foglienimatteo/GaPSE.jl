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


function integrand_ξ_multipole(s1, s, μ, effect::Function, cosmo::Cosmology;
     L::Integer = 0, use_windows::Bool = true, kwargs...)

     s2_value = s2(s1, s, μ)
     y_value = y(s1, s, μ)
     res = if use_windows == true
          ϕ_s2 = ϕ(s2_value, cosmo.s_min, cosmo.s_max)
          (ϕ_s2 > 0.0) || (return 0.0)
          #println("s1 = $s1 \\t s2 = $(s2(s1, s, μ)) \\t  y=$(y(s1, s, μ))")
          int = effect(s1, s2_value, y_value, cosmo; kwargs...)
          #println("int = $int")
          int .* (ϕ_s2 * spline_F(s / s1, μ, cosmo.windowF) * Pl(μ, L))
     else
          #println("s1 = $s1 \\t s2 = $(s2(s1, s, μ)) \\t  y=$(y(s1, s, μ))")
          int = effect(s1, s2_value, y_value, cosmo; kwargs...)
          #println("int = $int")
          #println( "Pl(μ, L) = $(Pl(μ, L))")
          int .* Pl(μ, L)
     end

     #println("res = $res")
     return (2.0 * L + 1.0) / 2.0 * res
end


function integrand_ξ_multipole(s1, s, μ, effect::String, cosmo::Cosmology; kwargs...)
     error = "$effect is not a valid GR effect name.\n" *
             "Valid GR effect names are the following:\n" *
             string(GaPSE.IMPLEMENTED_GR_EFFECTS .* " , "...)
     @assert (effect ∈ GaPSE.IMPLEMENTED_GR_EFFECTS) error

     return integrand_ξ_multipole(s1, s, μ, DICT_GR_ξs[effect], cosmo; kwargs...)
end



"""
     integrand_ξ_multipole(s1, s, μ, integrand::Function, cosmo::Cosmology;
          L::Integer = 0, 
          use_windows::Bool = true, 
          kwargs...) ::Float64

     integrand_ξ_multipole(s1, s, μ, effect::String, cosmo::Cosmology; 
          L::Integer = 0, 
          use_windows::Bool = true, 
          kwargs...) ::Float64

Return the integrand on ``\\mu = \\hat{\\mathbf{s}}_1 \\cdot \\hat{\\mathbf{s}}`` 
of the of the chosen correlation function term, i.e.
the following function ``f(s_1, s, \\mu)``:

```math
     f(s_1, s, \\mu) = \\xi (s_1, s_2, \\cos{\\theta}) 
          \\, \\mathcal{L}_L(\\mu) \\,  \\phi(s_2) \\, F\\left(\\frac{s}{s_1}, \\mu \\right)
```
where ``y =  \\cos{\\theta} = \\hat{\\mathbf{s}}_1 \\cdot \\hat{\\mathbf{s}}_2``,
``s = \\sqrt{s_1^2 + s_2^2 - 2 \\, s_1 \\, s_2 \\, y}`` and ``\\xi`` is the corresponding
CF effect.

In the former method you have to pass as an input the `integrand` function you want 
to integrate, while in the (recommended) latter one it's necessary to specify the
name of the CF term among the following: 

`$(string(GaPSE.IMPLEMENTED_GR_EFFECTS .* " , "...))`

to which correspond the following functions:

`$(string(string.(GaPSE.IMPLEMENTED_ξs) .* " , "...))`

In case `use_windows` is set to `false`, the window functions ``\\phi`` and ``F``
are removed, i.e is returned the following function ``f^{'}(s_1, s, \\mu)``:

```math
     f^{'}(s_1, s, \\mu) = \\xi (s_1, s_2, \\cos{\\theta}) 
          \\, \\mathcal{L}_L(\\mu) 
```

The function ``\\xi(s_1, s_2, \\cos{\\theta})`` is calculated
from, depending on the value of `effect`:
- `effect == auto_doppler` => [`ξ_Doppler`](@ref)
- `effect == auto_lensing` => [`ξ_Lensing`](@ref)
- `effect == auto_localgp` => [`ξ_LocalGP`](@ref)
- `effect == auto_integratedgp` => [`ξ_IntegratedGP`](@ref)
- `effect == doppler_lensing` => [`ξ_Doppler_lensing`](@ref)
- `effect == lensing_doppler` => [`ξ_Lensing_Doppler`](@ref)
- `effect == doppler_localgp` => [`ξ_Doppler_LocalGP`](@ref)
- `effect == localgp_doppler` => [`ξ_LocalGP_Doppler`](@ref)
- `effect == doppler_integratedgp` => [`ξ_Doppler_IntegratedGP`](@ref)
- `effect == integratedgp_doppler` => [`ξ_IntegratedGP_Doppler`](@ref)
- `effect == lensing_localgp` => [`ξ_Lensing_LocalGP`](@ref)
- `effect == localgp_lensing` => [`ξ_LocalGP_Lensing`](@ref)
- `effect == lensing_integratedgp` => [`ξ_Lensing_IntegratedGP`](@ref)
- `effect == integratedgp_lensing` => [`ξ_IntegratedGP_Lensing`](@ref)
- `effect == localgp_integratedgp` => [`ξ_LocalGP_IntegratedGP`](@ref)
- `effect == integratedgp_localgp` => [`ξ_IntegratedGP_LocalGP`](@ref)

Note that these is an internal conversion of coordiate sistems
from `(s1, s, μ)` to `(s1, s2, y)` thorugh the functions `y` and `s2`

## Inputs

- `s1`: the comoving distance where must be evaluated the integral

- `s`: the comoving distance from `s1` where must be evaluated the integral

- `μ`: the cosine between `s1` and `s` where must be evaluated the integral

- `cosmo::Cosmology`: cosmology to be used in this computation


## Optional arguments 

- `L::Integer = 0`: order of the Legendre polynomial to be used

- `use_windows::Bool = false`: tells if the integrand must consider the two
   window function ``\\phi`` and ``F``

- `kwargs...` : other keyword arguments that will be passed to the selected 
  GR TPCF effect (`ξ_Doppler`, `ξ_Lensing`, ...)


See also: [`integral_on_mu`](@ref), [`map_integral_on_mu`](@ref),
[`ξ_multipole`](@ref), [`map_ξ_multipole`](@ref),
[`spline_F`](@ref), [`ϕ`](@ref), [`Cosmology`](@ref), 
[`y`](@ref), [`s2`](@ref)
"""
integrand_ξ_multipole



##########################################################################################92



function ξ_multipole(
     s1, s, effect::Function, cosmo::Cosmology;
     L::Integer = 0,
     use_windows::Bool = true,
     enhancer::Float64 = 1e6,
     N_μs::Integer = 50,
     μ_atol::Float64 = 0.0,
     μ_rtol::Float64 = 1e-2,
     SPLINE::Bool = false,
     kwargs...)

     error = "$(string(effect)) is not a valid GR effect function.\n" *
             "Valid GR effect functions are the following:\n" *
             string(string.(IMPLEMENTED_ξs) .* " , "...)
     @assert (effect ∈ IMPLEMENTED_ξs) error

     orig_f(μ) = enhancer * integrand_ξ_multipole(s1, s, μ, effect, cosmo;
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



function ξ_multipole(s1, s, effect::String, cosmo::Cosmology; kwargs...)
     error = "$effect is not a valid GR effect name.\n" *
             "Valid GR effect names are the following:\n" *
             string(GaPSE.IMPLEMENTED_GR_EFFECTS .* " , "...)

     @assert (effect ∈ GaPSE.IMPLEMENTED_GR_EFFECTS) error
     ξ_multipole(s1, s, DICT_GR_ξs[effect], cosmo; kwargs...)
end



"""
     ξ_multipole(s1, s, effect::Function, cosmo::Cosmology; 
          L::Integer = 0, 
          enhancer::Float64 = 1e6,
          use_windows::Bool = true,
          μ_atol::Float64 = 1e-4,
          μ_rtol::Float64 = 1e-1, 
          kwargs...) ::Float64

     ξ_multipole(s1, s, effect::String, cosmo::Cosmology; 
          L::Integer = 0, 
          enhancer::Float64 = 1e6,
          use_windows::Bool = true,
          μ_atol::Float64 = 1e-4,
          μ_rtol::Float64 = 1e-1, 
          kwargs...) ::Float64

Evaluate the multipole of order `L` of the chosen correlation function term, 
through the `integral_on_mu` function.

In the former method you have to pass as an input the `integrand` function you want 
to integrate, while in the (recommended) latter one it's necessary to specify the
name of the CF term among the following:

`$(string(GaPSE.IMPLEMENTED_GR_EFFECTS .* " , "...))`

to which correspond the following functions:

`$(string(string.(GaPSE.IMPLEMENTED_ξs) .* " , "...))`

The function evaluated is then the following:

```math
\\xi_L(s_1, s, \\mu) = \\frac{2 L + 1}{2} \\int_{-1}^{+1} \\mathrm{d}\\mu \\; 
    \\xi (s_1, s_2, \\cos{\\theta}) \\, \\mathcal{L}_L(\\mu) \\,  \\times
\\begin{cases}  
    \\phi(s_2) \\, F\\left(\\frac{s}{s_1}, \\mu \\right) \\;,
        \\quad \\mathrm{use_windows == true}\\\\
    1\\;, \\quad \\mathrm{use_windows == false}
\\end{cases}
```
where ``y =  \\cos{\\theta} = \\hat{\\mathbf{s}}_1 
\\cdot \\hat{\\mathbf{s}}_2`` and ``\\xi`` is the chosen CF effect. 

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
   window function ``\\phi`` and ``F``

- `μ_atol::Float64 = 1e-3` and `μ_rtol::Float64 = 1e-3`: absolute and relative tolerance
  to be passed to `quadgk`; it's recommended not to set `μ_rtol < 1e-2` because
  of long time for evaluations

- `kwargs...` : other keyword arguments that will be passed to the selected 
  GR TPCF effect (`ξ_Doppler`, `ξ_Lensing`, ...)

See also: [`integrand_on_mu`](@ref),  [`integral_on_mu`](@ref),
[`map_ξ_multipole`](@ref), [`print_map_ξ_multipole`](@ref)
"""
ξ_multipole


##########################################################################################92


function map_ξ_multipole(cosmo::Cosmology,
     effect::Union{String,Function},
     v_ss = nothing;
     s_1 = nothing,
     pr::Bool = true,
     N_log::Integer = 1000,
     L::Integer = 0,
     kwargs...)

     s1 = isnothing(s_1) ? cosmo.s_eff : s_1

     t1 = time()
     ss = isnothing(v_ss) ? 10 .^ range(-1, 3, length = N_log) : v_ss
     xis = pr ? begin
          @showprogress "$effect, L=$L: " [
               ξ_multipole(s1, s, effect, cosmo; L = L, kwargs...) for s in ss
          ]
     end : [
          ξ_multipole(s1, s, effect, cosmo; L = L, kwargs...) for s in ss
     ]

     t2 = time()
     pr && println("\ntime needed for map_ξ_multipole for $effect " *
                   "[in s] = $(@sprintf("%.5f", t2-t1)) ")
     return (ss, xis)
end


##########################################################################################92



function print_map_ξ_multipole(
     cosmo::Cosmology,
     out::String,
     effect::Union{String,Function},
     v_ss = nothing;
     s_1 = nothing,
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
               for key in keys(kwargs)
                    println(io, "# \t\t$(key) = $(kwargs[key])")
               end
          end
          isnothing(s_1) || println(io, "#\n# NOTE: the computation is done not in " *
                                        "s1 = s_eff, because you specified in input s1 = $s_1 !")
          println(io, "# ")
          println(io, "# s [Mpc/h_0] \t \t xi")
          for (s, xi) in zip(vec[1], vec[2])
               println(io, "$s \t $xi")
          end
     end
end

