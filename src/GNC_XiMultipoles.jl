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
# MERCHANTABILITY or FITNEss FOR A PARTICULAR PURPOSE. See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with GaPSE. If not, see <http://www.gnu.org/licenses/>.
#


function integrand_ξ_GNC_multipole(s1, s, μ, effect::Function, cosmo::Cosmology;
     L::Int = 0, use_windows::Bool = true, kwargs...)

     s2_value = s2(s1, s, μ)
     y_value = y(s1, s, μ)

     res = if use_windows == true
          #println("s1 = $s1 \\t s2 = $(s2(s1, s, μ)) \\t  y=$(y(s1, s, μ))")
          int = effect(s1, s2_value, y_value, cosmo; kwargs...)
          #print(" F_val = ", spline_integrF(s, μ, cosmo.windowFint))
          #isapprox(cost, 3.1749350588720085e10; rtol=1e-5) || print(" \t cost = " , cosmo.WFI_norm)
          #isapprox(Pl(μ, L), 1.0; rotl=1e-4) || println(" \t Lp = " , Pl(μ, L) )
          #println("IF = $val")

          int .* (spline_integrF(s, μ, cosmo.windowFint)/cosmo.WFI_norm * Pl(μ, L))
          #=
          ϕ_s2 = ϕ(s2_value, cosmo.s_min, cosmo.s_max)
          (ϕ_s2 > 0.0) || (return 0.0)
          int = effect(s1, s2_value, y_value, cosmo; kwargs...)
          int .* (ϕ_s2 * spline_F(s / s1, μ, cosmo.windowF) * Pl(μ, L))
          =#
     else
          int = effect(s1, s2_value, y_value, cosmo; kwargs...)
          int .* Pl(μ, L)
     end

     return (2.0 * L + 1.0) / 2.0 * res
end


function integrand_ξ_GNC_multipole(s1, s, μ, effect::String, cosmo::Cosmology; kwargs...)
     error = "$effect is not a valid GR effect name.\n" *
             "Valid GR effect names are the following:\n" *
             string(GaPSE.GR_EFFECTS_GNC .* " , "...)
     @assert (effect ∈ GaPSE.GR_EFFECTS_GNC) error

     return integrand_ξ_GNC_multipole(s1, s, μ, DICT_GR_ξs_GNC[effect], cosmo; kwargs...)
end



"""
     integrand_ξ_GNC_multipole(s1, s, μ, effect::Function, cosmo::Cosmology;
          L::Int = 0, 
          use_windows::Bool = true, 
          kwargs...) ::Float64

     integrand_ξ_GNC_multipole(s1, s, μ, effect::String, cosmo::Cosmology; 
          L::Int = 0, 
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

In the former method you have to pass as an input the `effect` function you want 
to integrate, while in the (recommended) latter one it's necessary to specify the
name of the CF term among the following: 

`$(string(GaPSE.GR_EFFECTS_GNC .* " , "...))`

to which correspond the following functions:

`$(string(string.(GaPSE.VEC_ξs_GNC) .* " , "...))`

The window functions ``F(x, \\mu)`` and ``\\phi(s)`` are calculated for the given
Cosmology `cosmo` through the functions `spline_F` and `ϕ` respectivelly.

In case `use_windows` is set to `false`, the window functions ``\\phi`` and ``F``
are removed, i.e is returned the following function ``f^{'}(s_1, s, \\mu)``:

```math
     f^{'}(s_1, s, \\mu) = \\xi (s_1, s_2, \\cos{\\theta}) 
          \\, \\mathcal{L}_L(\\mu) 
```

The function ``\\xi(s_1, s_2, \\cos{\\theta})`` is calculated
from, depending on the value of `effect`:
- `effect == auto_doppler` => [`ξ_GNC_Doppler`](@ref)
- `effect == auto_lensing` => [`ξ_GNC_Lensing`](@ref)
- `effect == auto_localgp` => [`ξ_GNC_LocalGP`](@ref)
- `effect == auto_integratedgp` => [`ξ_GNC_IntegratedGP`](@ref)
- `effect == doppler_lensing` => [`ξ_GNC_Doppler_lensing`](@ref)
- `effect == lensing_doppler` => [`ξ_GNC_Lensing_Doppler`](@ref)
- `effect == doppler_localgp` => [`ξ_GNC_Doppler_LocalGP`](@ref)
- `effect == localgp_doppler` => [`ξ_GNC_LocalGP_Doppler`](@ref)
- `effect == doppler_integratedgp` => [`ξ_GNC_Doppler_IntegratedGP`](@ref)
- `effect == integratedgp_doppler` => [`ξ_GNC_IntegratedGP_Doppler`](@ref)
- `effect == lensing_localgp` => [`ξ_GNC_Lensing_LocalGP`](@ref)
- `effect == localgp_lensing` => [`ξ_GNC_LocalGP_Lensing`](@ref)
- `effect == lensing_integratedgp` => [`ξ_GNC_Lensing_IntegratedGP`](@ref)
- `effect == integratedgp_lensing` => [`ξ_GNC_IntegratedGP_Lensing`](@ref)
- `effect == localgp_integratedgp` => [`ξ_GNC_LocalGP_IntegratedGP`](@ref)
- `effect == integratedgp_localgp` => [`ξ_GNC_IntegratedGP_LocalGP`](@ref)

Note that these is an internal conversion of coordiate sistems
from `(s1, s, μ)` to `(s1, s2, y)` thorugh the functions `y` and `s2`

## Inputs

- `s1`: the comoving distance where must be evaluated the integral

- `s`: the comoving distance from `s1` where must be evaluated the integral

- `μ`: the cosine between `s1` and `s` where must be evaluated the integral

- `cosmo::Cosmology`: cosmology to be used in this computation


## Optional arguments 

- `L::Int = 0`: order of the Legendre polynomial to be used

- `use_windows::Bool = false`: tells if the integrand must consider the two
   window function ``\\phi`` and ``F``

- `kwargs...` : other keyword arguments that will be passed to the selected 
  GR TPCF effect (`ξ_GNC_Doppler`, `ξ_GNC_Lensing`, ...)


See also: [`ξ_GNC_multipole`](@ref), [`map_ξ_GNC_multipole`](@ref),
[`print_map_ξ_GNC_multipole`](@ref),
[`spline_F`](@ref), [`ϕ`](@ref), [`Cosmology`](@ref), 
[`y`](@ref), [`s2`](@ref), [`GR_EFFECTS_GNC`](@ref)
"""
integrand_ξ_GNC_multipole


##########################################################################################92




function ξ_GNC_multipole(
     s1, s, effect::Function, cosmo::Cosmology;
     alg::Symbol = :lobatto, obs::Union{Bool, Symbol} = :noobsvel,
     enhancer::Float64 = 1e6, N_lob::Int = 100, L::Int = 0,
     N_trap::Int = 50, atol_quad::Float64 = 0.0, rtol_quad::Float64 = 1e-2,
     kwargs...)

     error = "$(string(effect)) is not a valid GR effect function for galaxy number counts.\n" *
             "Valid GR effect functions for galaxy number counts are the following:\n" *
             string(string.(VEC_ξs_GNC) .* " , "...)
     @assert (effect ∈ VEC_ξs_GNC) error
     @assert typeof(obs) == Bool || obs ∈ VALID_OBS_VALUES ":$obs is not a valid Symbol for \"obs\"; they are: \n\t"*
          "$(":".*string.(VALID_OBS_VALUES) .* vcat([" , " for i in 1:length(VALID_OBS_VALUES)-1], " .")... )" 
     @assert alg ∈ VALID_INTEGRATION_ALGORITHM ":$alg is not a valid Symbol for \"alg\"; they are: \n\t"*
          "$(":".*string.(VALID_INTEGRATION_ALGORITHM) .* vcat([" , " for i in 1:length(VALID_INTEGRATION_ALGORITHM)-1], " .")... )" 

     @assert N_trap > 2 "N_trap must be >2,  N_trap = $N_trap is not!"  
     @assert N_lob > 2 "N_lob must be >2,  N_lob = $N_lob is not!"  
     @assert atol_quad ≥ 0.0 "atol_quad must be ≥ 0.0,  atol_quad = $atol_quad is not!"  
     @assert rtol_quad ≥ 0.0  "rtol_trap must be ≥ 0.0,  rtol_quad = $rtol_quad is not!" 
     @assert L ≥ 0 "L must be ≥ 0, L = $L is not!"       


     orig_f(μ) = enhancer * integrand_ξ_GNC_multipole(s1, s, μ, effect, cosmo; L = L, obs = obs, kwargs...)
     
     int = if alg == :lobatto 
               xs, ws = gausslobatto(N_lob)
               dot(ws, orig_f.(xs))

          elseif alg == :quad 
               quadgk(μ -> orig_f(μ), -1.0, 1.0; atol = atol_quad, rtol = rtol_quad)[1]

          elseif  alg == :trap
               μs = union(
                    range(-1.0, -0.98, length = Int(ceil(N_trap/3) + 1)),
                    range(-0.98, 0.98, length = Int(ceil(N_trap/3) + 1)),
                    range(0.98, 1.0, length = Int(ceil(N_trap/3) + 1))
               )
               #μs = range(-1.0 + 1e-6, 1.0 - 1e-6, length=N_trap)
               orig_fs = orig_f.(μs)
               trapz(μs, orig_fs)

          else
               throw(AssertionError("how the hell did you arrive here?"))
          end

     return int / enhancer
end



function ξ_GNC_multipole(s1, s, effect::String, cosmo::Cosmology; kwargs...)
     error = "$effect is not a valid GR effect name for galaxy number counts.\n" *
             "Valid GR effect names for galaxy number counts are the following:\n" *
             string(GaPSE.GR_EFFECTS_GNC .* " , "...)

     @assert (effect ∈ GaPSE.GR_EFFECTS_GNC) error
     ξ_GNC_multipole(s1, s, DICT_GR_ξs_GNC[effect], cosmo; kwargs...)
end



"""
     ξ_GNC_multipole(s1, s, effect::Function, cosmo::Cosmology; 
          L::Int = 0, 
          enhancer::Float64 = 1e6,
          use_windows::Bool = true,
          atol_quad::Float64 = 1e-4,
          rtol_quad::Float64 = 1e-1, 
          kwargs...) ::Float64

     ξ_GNC_multipole(s1, s, effect::String, cosmo::Cosmology; 
          L::Int = 0, 
          enhancer::Float64 = 1e6,
          use_windows::Bool = true,
          atol_quad::Float64 = 1e-4,
          rtol_quad::Float64 = 1e-1, 
          kwargs...) ::Float64

Evaluate the multipole of order `L` of the chosen correlation function term, 
through the `integrand_ξ_GNC_multipole` function.

In the former method you have to pass as an input the `effect` function you want 
to integrate, while in the (recommended) latter one it's necessary to specify the
name of the CF term among the following:

`$(string(GaPSE.GR_EFFECTS_GNC .* " , "...))`

to which correspond the following functions:

`$(string(string.(GaPSE.VEC_ξs_GNC) .* " , "...))`

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

The window functions ``F(x, \\mu)`` and ``\\phi(s)`` are calculated for the given
Cosmology `cosmo` through the functions `spline_F` and `ϕ` respectivelly.

The function ``\\xi(s_1, s_2, \\cos{\\theta})`` is calculated
from, depending on the value of `effect`:
- `effect == auto_doppler` => [`ξ_GNC_Doppler`](@ref)
- `effect == auto_lensing` => [`ξ_GNC_Lensing`](@ref)
- `effect == auto_localgp` => [`ξ_GNC_LocalGP`](@ref)
- `effect == auto_integratedgp` => [`ξ_GNC_IntegratedGP`](@ref)
- `effect == doppler_lensing` => [`ξ_GNC_Doppler_lensing`](@ref)
- `effect == lensing_doppler` => [`ξ_GNC_Lensing_Doppler`](@ref)
- `effect == doppler_localgp` => [`ξ_GNC_Doppler_LocalGP`](@ref)
- `effect == localgp_doppler` => [`ξ_GNC_LocalGP_Doppler`](@ref)
- `effect == doppler_integratedgp` => [`ξ_GNC_Doppler_IntegratedGP`](@ref)
- `effect == integratedgp_doppler` => [`ξ_GNC_IntegratedGP_Doppler`](@ref)
- `effect == lensing_localgp` => [`ξ_GNC_Lensing_LocalGP`](@ref)
- `effect == localgp_lensing` => [`ξ_GNC_LocalGP_Lensing`](@ref)
- `effect == lensing_integratedgp` => [`ξ_GNC_Lensing_IntegratedGP`](@ref)
- `effect == integratedgp_lensing` => [`ξ_GNC_IntegratedGP_Lensing`](@ref)
- `effect == localgp_integratedgp` => [`ξ_GNC_LocalGP_IntegratedGP`](@ref)
- `effect == integratedgp_localgp` => [`ξ_GNC_IntegratedGP_LocalGP`](@ref)

Note that these is an internal conversion of coordiate sistems
from `(s1, s, μ)` to `(s1, s2, y)` thorugh the functions `y` and `s2`

## Inputs

- `s1`: the comoving distance where must be evaluated the integral

- `s`: the comoving distance from `s1` where must be evaluated the integral

- `cosmo::Cosmology`: cosmology to be used in this computation


## Optional arguments

- `L::Int = 0`: order of the Legendre polynomial to be used

- `enhancer::Float64 = 1e6`: just a float number used in order to deal better with small numbers; 
  the returned value is NOT modified by this value, because after a multiplication
  the internal result is divided by `enhancer`.

- `use_windows::Bool = false`: tells if the integrand must consider the two
   window function ``\\phi`` and ``F``

- `atol_quad::Float64 = 1e-3` and `rtol_quad::Float64 = 1e-3`: absolute and relative tolerance
  to be passed to `quadgk`; it's recommended not to set `rtol_quad < 1e-2` because
  of long time for evaluations

- `kwargs...` : other keyword arguments that will be passed to the selected 
  GR TPCF effect (`ξ_GNC_Doppler`, `ξ_GNC_Lensing`, ...)

See also: [`integrand_ξ_GNC_multipole`](@ref), 
[`map_ξ_GNC_multipole`](@ref), [`print_map_ξ_GNC_multipole`](@ref)
[`spline_F`](@ref), [`ϕ`](@ref), [`Cosmology`](@ref), 
[`y`](@ref), [`s2`](@ref), [`GR_EFFECTS_GNC`](@ref)
"""
ξ_GNC_multipole


##########################################################################################92


"""
     map_ξ_GNC_multipole(cosmo::Cosmology,
          effect::Union{String,Function},
          ss = nothing;
          s1 = nothing,
          pr::Bool = true,
          N_log::Int = 1000,
          L::Int = 0,
          kwargs... ) ::Tuple{Vector{Float64}, Vector{Float64}}

Evaluate the multipole of order `L` of the chosen correlation function term, 
through the `ξ_GNC_multipole` function, for all the `s` values stored inside `ss`.
If `ss = nothing`, it is set `ss = 10 .^ range(0, 3, length = N_log)`.
If `s1 = nothing`, it is set `s1 = cosmo.s_eff`.

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
\\cdot \\hat{\\mathbf{s}}_2`` and ``\\xi`` is the chosen CF effect, for all the 
comoving distances `s` inside `ss`. 

The window functions ``F(x, \\mu)`` and ``\\phi(s)`` are calculated for the given
Cosmology `cosmo` through the functions `spline_F` and `ϕ` respectivelly.

## Inputs

- `cosmo::Cosmology` : cosmology to be used in this computation

- `effect::Union{String,Function}` : the GR effect TPCF you want to consider; you may
  specify the name of the effect as one of the following strings (recommended):

  `$(string(GaPSE.GR_EFFECTS_GNC .* " , "...))`
  
  or directly the name of the function among the following: 
  
  `$(string(string.(GaPSE.VEC_ξs_GNC) .* " , "...))`

- ``ss` : vector/range of `s` values where the function must be evaluated; if `ss = nothing`, 
  it is set `ss = 10 .^ range(0, 3, length = N_log)`. This is why it is returned 
  also the vector of the "input" values.


## Optional arguments

- `s1 = nothing` : comoving distance from the observer where the TPCF should be evaluated;
  if `s1 = nothing`, it is automatically set `s1 = cosmo.s_eff` from the given input `Cosmology`.

- `L::Int = 0`: order of the Legendre polynomial to be used

- `pr::Bool = true` : do you want the progress bar showed on screen, in order to 
  check the time needed for the computation? (`true` recommended)

- `N_log::Int = 1000` : number of points to be used in the default logaritmically-spaced 
  range for `ss`, i.e. `range(0, 3, N_log)`; it is ignored if `ss ≠ nothing` 

- `kwargs...` : other keyword arguments that will be passed to `ξ_GNC_multipole`

# Returns

A `Tuple{Vector{Float64}, Vector{Flaot64}}`, which has as first element the `ss` vector
and as second one the corresponding ξ value evaluated.

See also: [`integrand_ξ_GNC_multipole`](@ref), [`ξ_GNC_multipole`](@ref),
[`print_map_ξ_GNC_multipole`](@ref),
[`spline_F`](@ref), [`ϕ`](@ref), [`Cosmology`](@ref), 
[`y`](@ref), [`s2`](@ref), [`GR_EFFECTS_GNC`](@ref)
"""
function map_ξ_GNC_multipole(cosmo::Cosmology,
     effect::Union{String,Function}, ss = nothing;
     s1 = nothing, L::Int = 0, alg::Symbol = :lobatto,
     obs::Union{Bool, Symbol} = :noobsvel,
     N_lob::Int = 100, N_trap::Int = 50, 
     atol_quad::Float64 = 0.0, rtol_quad::Float64 = 1e-2, 
     pr::Bool = true, N_log::Int = 1000, sum_xi::Bool = false, 
     enhancer::Float64 = 1e6,
     kwargs...)

     @assert typeof(obs) == Bool || obs ∈ VALID_OBS_VALUES ":$obs is not a valid Symbol for \"obs\"; they are: \n\t"*
          "$(":".*string.(VALID_OBS_VALUES) .* vcat([" , " for i in 1:length(VALID_OBS_VALUES)-1], " .")... )" 
     @assert alg ∈ VALID_INTEGRATION_ALGORITHM ":$alg is not a valid Symbol for \"alg\"; they are: \n\t"*
          "$(":".*string.(VALID_INTEGRATION_ALGORITHM) .* vcat([" , " for i in 1:length(VALID_INTEGRATION_ALGORITHM)-1], " .")... )" 

     @assert N_trap > 2 "N_trap must be >2,  N_trap = $N_trap is not!"  
     @assert N_lob > 2 "N_lob must be >2,  N_lob = $N_lob is not!"  
     @assert atol_quad ≥ 0.0 "atol_quad must be ≥ 0.0,  atol_quad = $atol_quad is not!"  
     @assert rtol_quad ≥ 0.0  "rtol_trap must be ≥ 0.0,  rtol_quad = $rtol_quad is not!"
     @assert L ≥ 0 "L must be ≥ 0, L = $L is not!"      

     t1 = time()


     s_1 = isnothing(s1) ? cosmo.s_eff : s1
     v_ss = isnothing(ss) ? 10 .^ range(0, log10(2*cosmo.s_max), length = N_log) : ss
     xis = [0.0 for i in 1:length(v_ss)]

     orig_f(μ, s) = enhancer * integrand_ξ_GNC_multipole(s_1, s, μ, effect, cosmo; 
          L = L, obs = obs, kwargs...)

     if alg == :lobatto 
          μs, ws = gausslobatto(N_lob)

          global xis = pr ? begin
               @showprogress "$effect, L=$L: " [
                    dot(ws, [orig_f(μ, s) for μ in μs])/enhancer for s in v_ss
                    ]
               end : [
                    dot(ws, [orig_f(μ, s) for μ in μs])/enhancer for s in v_ss
               ]

     elseif alg == :quad 

          global xis = pr ? begin
               @showprogress "$effect, L=$L: " [
                    quadgk(μ -> orig_f(μ, s), -1.0, 1.0; 
                         atol = atol_quad, rtol = rtol_quad)[1]/enhancer for s in v_ss
                    ]
               end : [
                    quadgk(μ -> orig_f(μ, s), -1.0, 1.0; 
                         atol = atol_quad, rtol = rtol_quad)[1]/enhancer for s in v_ss
               ]

     elseif  alg == :trap

          μs = union(
               range(-1.0, -0.98, length = Int(ceil(N_trap/3) + 1)),
               range(-0.98, 0.98, length = Int(ceil(N_trap/3) + 1)),
               range(0.98, 1.0, length = Int(ceil(N_trap/3) + 1))
          )
          #μs = range(-1.0 + 1e-6, 1.0 - 1e-6, length=N_trap)

          global xis = pr ? begin
               @showprogress "$effect, L=$L: " [
                    trapz(μs, [orig_f(μ, s) for μ in μs])/enhancer for s in v_ss
                    ]
               end : [
                    trapz(μs, [orig_f(μ, s) for μ in μs])/enhancer for s in v_ss
               ]

     else
          throw(AssertionError("how the hell did you arrive here?"))
     end

     #=
     xis = pr ? begin
          @showprogress "$effect, L=$L: " [
               ξ_GNC_multipole(s_1, s, effect, cosmo; L=L, kwargs...) for s in v_ss
          ]
     end : [
          ξ_GNC_multipole(s_1, s, effect, cosmo; L=L, kwargs...) for s in v_ss
     ]
     =#

     t2 = time()
     (!sum_xi) && (pr) && println("\ntime needed for map_ξ_GNC_multipole for $effect, L=$L " *
                   "[in s] = $(@sprintf("%.5f", t2-t1)) \n")
     return (v_ss, xis)
end


##########################################################################################92


"""
     print_map_ξ_GNC_multipole(
          cosmo::Cosmology,
          out::String,
          effect::Union{String,Function},
          ss = nothing;
          s1 = nothing,
          kwargs...)

Evaluate the multipole of order `L` of the chosen correlation function term, 
through the `ξ_GNC_multipole` function, for all the `s` values stored inside `ss`, and
print the results (with all the options used) in a file named `out`.
If `ss = nothing`, it is set `ss = 10 .^ range(0, 3, length = N_log)`.
If `s1 = nothing`, it is set `s1 = cosmo.s_eff`.

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
\\cdot \\hat{\\mathbf{s}}_2`` and ``\\xi`` is the chosen CF effect, for all the 
comoving distances `s` inside `ss`. 

The window functions ``F(x, \\mu)`` and ``\\phi(s)`` are calculated for the given
Cosmology `cosmo` through the functions `spline_F` and `ϕ` respectivelly.

## Inputs

- `cosmo::Cosmology` : cosmology to be used in this computation

- `effect::Union{String,Function}` : the GR effect TPCF you want to consider; you may
  specify the name of the effect as one of the following strings (recommended):

  `$(string(GaPSE.GR_EFFECTS_GNC .* " , "...))`
  
  or directly the name of the function among the following: 
  
  `$(string(string.(GaPSE.VEC_ξs_GNC) .* " , "...))`

- `out::String` : name of the file where the results must be stored.

- ``ss` : vector/range of `s` values where the function must be evaluated; if `ss = nothing`, 
  it is set `ss = 10 .^ range(0, 3, length = N_log)`.


## Optional arguments

- `s1 = nothing` : comoving distance from the observer where the TPCF should be evaluated;
  if `s1 = nothing`, it is automatically set `s1 = cosmo.s_eff` from the given input `Cosmology`.

- `kwargs...` : other keyword arguments that will be passed to `map_ξ_GNC_multipole`

See also: [`integrand_ξ_GNC_multipole`](@ref), [`ξ_GNC_multipole`](@ref),
[`map_ξ_GNC_multipole`](@ref), 
[`spline_F`](@ref), [`ϕ`](@ref), [`Cosmology`](@ref), 
[`y`](@ref), [`s2`](@ref), [`GR_EFFECTS_GNC`](@ref)
"""
function print_map_ξ_GNC_multipole(
          cosmo::Cosmology, out::String, 
          effect::Union{String,Function}, ss = nothing;
          s1 = nothing, L::Int = 0, kwargs...)

     check_parent_directory(out)
     check_namefile(out)

     s_1 = isnothing(s1) ? cosmo.s_eff : s1
     t1 = time()
     vec = map_ξ_GNC_multipole(cosmo, effect, ss; s1 = s_1, L = L, kwargs...)
     t2 = time()

     isfile(out) && run(`rm $out`)
     open(out, "w") do io
          println(io, BRAND)

          println(io, "#\n# This is an integration map on mu of the ξ L=$L multipole $effect GR effect")
          println(io, "# concerning the relativistic galaxy number counts.")
          parameters_used(io, cosmo; logo = false)
          println(io, "# computational time needed (in s) : $(@sprintf("%.4f", t2-t1))")
          print(io, "# kwards passed: ")

          if isempty(kwargs)
               println(io, "none")
          else
               print(io, "\n")
               println(io, "# \t\tL = $L")
               for key in keys(kwargs)
                    println(io, "# \t\t$(key) = $(kwargs[key])")
               end
          end
          isnothing(s1) || println(io, "#\n# NOTE: the computation is done not in " *
                                       "s1 = s_eff, because you specified in input s1 = $s1 !")
          println(io, "# ")
          println(io, "# s [Mpc/h_0] \t \t xi")
          for (s, xi) in zip(vec[1], vec[2])
               println(io, "$s \t $xi")
          end
     end
end

