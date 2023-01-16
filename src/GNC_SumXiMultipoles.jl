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



"""
     sum_ξ_GNC_multipole(s1, s, cosmo::Cosmology;
          s1 = nothing, L::Int = 0, alg::Symbol = :lobatto,
          obs::Union{Bool,Symbol} = :noobsvel,
          N_lob::Int = 100, N_trap::Int = 50,
          atol_quad::Float64 = 0.0, rtol_quad::Float64 = 1e-2,
          enhancer::Float64=1e6, N_log::Int = 1000, 
          pr::Bool = true,
          kwargs...) ::Tuple{Float64, Vector{Float64}}

Evaluate the multipole of order `L` of all the Galaxy Number Counts (GNC) 
Two-Point Correlation Function (TPCF) multipoles and their sum
in the comoving distance `s1` and a comoving distance `s` from it 
for the input `cosmo::Cosmology`.
We remember that all the distances are measured in ``h_0^{-1}\\mathrm{Mpc}``.

This function makes a for-loop on the `GaPSE.GR_EFFECTS_GNC` strings, calling 
`ξ_GNC_multipole` for each of them. To each string corresponds pretty intuitively one of the 
25 GNC effects. They are currently, in order:

`$(string(GR_EFFECTS_GNC .* " , "...))`

## Inputs

- `s1`: the comoving distance where must be evaluated the integral

- `s`: the comoving distance from `s1` where must be evaluated the integral

- `cosmo::Cosmology`: cosmology to be used in this computation

## Optional arguments

This function recall internally `ξ_GNC_multipole`, so the kwargs are the same; we report them
for comfortness:

- `L::Int = 0`: order of the Legendre polynomial to be used

- `alg::Symbol = :trap` : algorithm to be used for the integration; the valid options 
  are (other values will lead to `AssertionError`):
  - `:quad` -> the integration over ``\\mu`` will be preformed through the Julia function `quadgk` 
  from the [`QuadGK.jl`](https://github.com/JuliaMath/QuadGK.jl) Julia package, that uses an adaptive 
  Gauss-Kronrod quadrature.
  - `:trap` -> the integration over ``\\mu`` will be preformed through the Julia function `trapz` 
  from the [`Trapz.jl`](https://github.com/francescoalemanno/Trapz.jl) Julia package, that uses the
  simple trapezoidal rulae.
  - `:lobatto` -> the integration over ``\\mu`` will be preformed through the Julia function `gausslobatto` 
  from the [`FastGaussQuadrature.jl`](https://github.com/JuliaApproximation/FastGaussQuadrature.jl) Julia package, 
  that uses the Gauss-Lobatto quadrature. 
  WE RECOMMEND TO USE `:quad` FOR MONOPOLES AND `:lobatto` FOR HIGHER ORDER MULTIPOLES!
  
- `use_windows::Bool = false`: tells if the integrand must consider the two
   window function ``\\phi`` and ``\\mathcal{F}``

- `obs::Union{Bool,Symbol} = :noobsvel` : do you want to consider the observer terms in the computation of the 
  chosen GNC TPCF effect?
  - `:yes` or `true` -> all the observer effects will be considered
  - `:no` or `false` -> no observer term will be taken into account
  - `:noobsvel` -> the observer terms related to the observer velocity (that you can find in the CF concerning Doppler)
    will be neglected, the other ones will be taken into account

- `N_lob::Int = 100` : number of points to be used in the sampling made by the function `trapz`.
  Note that these options will have an effect only if you se `alg = :quad`.

- `N_trap::Int = 200` : number of points to be used in the sampling made by the function `trapz`.
  Note that these options will have an effect only if you se `alg = :quad`.

- `atol_quad::Float64 = 0.0` and `rtol_quad::Float64 = 1e-2`: absolute and relative tolerance
  to be passed to the function `quadgk`; it's recommended not to set `rtol_quad < 1e-2` 
  because the time for evaluation increase quickly.
  Note that these options will have an effect only if you se `alg = :quad`.

- `enhancer::Float64 = 1e6`: just a float number used in order to deal better with small numbers; 
  the returned value is NOT modified by this value, because after a multiplication
  the internal result is divided by `enhancer`.

- `kwargs...` : other keyword arguments that will be passed to ALL the 
  GNC TPCF effect (`ξ_GNC_Doppler`, `ξ_GNC_Lensing`, ...); if one of them has that keyword argument,
  it will use the given value, otherwise it will be unaffected.

## Returns

A tuple containing:
- the sum of all the ξ multipoles as first element
- a `Vector{Float64}` with all the values of each ξ; they are ordered
  following `GR_EFFECTS_GNC`


See also: [`integrand_ξ_GNC_multipole`](@ref), [`ξ_GNC_multipole`](@ref),
[`map_sum_ξ_GNC_multipole`](@ref), [`print_map_sum_ξ_GNC_multipole`](@ref),
[`Cosmology`](@ref), [`GR_EFFECTS_GNC`](@ref)
"""
function sum_ξ_GNC_multipole(s1, s, cosmo::Cosmology; kwargs...)
     ALL = [ξ_GNC_multipole(s1, s, effect, cosmo; specif_kwargs_GNC(effect, kwargs)...)
            for effect in GaPSE.GR_EFFECTS_GNC]

     return sum(ALL), ALL
end



##########################################################################################92



"""
     map_sum_ξ_GNC_multipole(
          effect::Union{String,Function}, ss = nothing;
          s1 = nothing, L::Int = 0, alg::Symbol = :lobatto,
          obs::Union{Bool,Symbol} = :noobsvel,
          N_lob::Int = 100, N_trap::Int = 50,
          atol_quad::Float64 = 0.0, rtol_quad::Float64 = 1e-2,
          enhancer::Float64=1e6, N_log::Int = 1000, 
          pr::Bool = true,
          kwargs...) ::Tuple{Vector{Float64}, Vector{Float64}, Vector{Vector{Float64}}}

Evaluate the multipole of order `L` of all the Galaxy Number Counts (GNC) 
Two-Point Correlation Function (TPCF) multipoles and their sum
in the comoving distance `s1`,  for all the comoving distances stored inside `ss` (representing 
the comoving distance from `s1`) for the input `cosmo::Cosmology`.
If `ss = nothing`, it is set `ss = 10 .^ range(0, log10(2 * cosmo.s_max), length=N_log)`.
If `s1 = nothing`, it is set `s1 = cosmo.s_eff`.
We remember that all the distances are measured in ``h_0^{-1}\\mathrm{Mpc}``.

This function makes a for-loop on the `GaPSE.GR_EFFECTS_GNC` strings, calling 
`map_ξ_GNC_multipole` for each of them. To each string corresponds pretty intuitively one of the 
25 GNC effects. They are currently, in order:

`$(string(GR_EFFECTS_GNC .* " , "...))`

## Inputs

- `cosmo::Cosmology` : cosmology to be used in this computation

- `ss` : vector/range of `s` values where the function must be evaluated; if `ss = nothing`, 
  it is set `ss = 10 .^ range(0, log10(2 * cosmo.s_max), length=N_log)`. This is why it is returned 
  also the vector of the "input" values.

## Optional arguments

This function recall internally `map_ξ_GNC_multipole`, so the kwargs are the same; we report them
for comfortness:

- `s1 = nothing` : comoving distance from the observer where the TPCF should be evaluated;
  if `s1 = nothing`, it is automatically set `s1 = cosmo.s_eff` from the given input `cosmo::Cosmology`.

- `L::Int = 0`: order of the Legendre polynomial to be used

- `alg::Symbol = :trap` : algorithm to be used for the integration; the valid options 
  are (other values will lead to `AssertionError`):
  - `:quad` -> the integration over ``\\mu`` will be preformed through the Julia function `quadgk` 
  from the [`QuadGK.jl`](https://github.com/JuliaMath/QuadGK.jl) Julia package, that uses an adaptive 
  Gauss-Kronrod quadrature.
  - `:trap` -> the integration over ``\\mu`` will be preformed through the Julia function `trapz` 
  from the [`Trapz.jl`](https://github.com/francescoalemanno/Trapz.jl) Julia package, that uses the
  simple trapezoidal rulae.
  - `:lobatto` -> the integration over ``\\mu`` will be preformed through the Julia function `gausslobatto` 
  from the [`FastGaussQuadrature.jl`](https://github.com/JuliaApproximation/FastGaussQuadrature.jl) Julia package, 
  that uses the Gauss-Lobatto quadrature. 
  WE RECOMMEND TO USE `:quad` FOR MONOPOLES AND `:lobatto` FOR HIGHER ORDER MULTIPOLES!
  
- `use_windows::Bool = false`: tells if the integrand must consider the two
   window function ``\\phi`` and ``\\mathcal{F}``

- `obs::Union{Bool,Symbol} = :noobsvel` : do you want to consider the observer terms in the computation of the 
  chosen GNC TPCF effect?
  - `:yes` or `true` -> all the observer effects will be considered
  - `:no` or `false` -> no observer term will be taken into account
  - `:noobsvel` -> the observer terms related to the observer velocity (that you can find in the CF concerning Doppler)
    will be neglected, the other ones will be taken into account

- `N_lob::Int = 100` : number of points to be used in the sampling made by the function `trapz`.
  Note that these options will have an effect only if you se `alg = :quad`.

- `N_trap::Int = 200` : number of points to be used in the sampling made by the function `trapz`.
  Note that these options will have an effect only if you se `alg = :quad`.

- `atol_quad::Float64 = 0.0` and `rtol_quad::Float64 = 1e-2`: absolute and relative tolerance
  to be passed to the function `quadgk`; it's recommended not to set `rtol_quad < 1e-2` 
  because the time for evaluation increase quickly.
  Note that these options will have an effect only if you se `alg = :quad`.

- `enhancer::Float64 = 1e6`: just a float number used in order to deal better with small numbers; 
  the returned value is NOT modified by this value, because after a multiplication
  the internal result is divided by `enhancer`.

- `N_log::Int = 1000` : number of points to be used in the default logaritmically-spaced 
  range for `ss`, i.e. `range(0, log10(2 * cosmo.s_max), length=N_log)`; it is ignored if `ss ≠ nothing` 

- `pr::Bool = true` : do you want the progress bar showed on screen, in order to 
  check the time needed for the computation? (`true` recommended)

- `kwargs...` : other keyword arguments that will be passed to ALL the 
  GNC TPCF effect (`ξ_GNC_Doppler`, `ξ_GNC_Lensing`, ...); if one of them has that keyword argument,
  it will use the given value, otherwise it will be unaffected.

## Returns

A tuple containing:
- as first element, the vector `ss` itself;
- as second one, the  `Vector{Float64}` of the sum of all the ξ multipoles;
- as third one, a `Vector{Vector{Float64}}` with all the values of each ξ; they are ordered
  following `GR_EFFECTS_GNC`


See also: [`map_ξ_GNC_multipole`](@ref),
[`sum_ξ_GNC_multipole`](@ref), [`print_map_sum_ξ_GNC_multipole`](@ref),
[`Cosmology`](@ref), [`GR_EFFECTS_GNC`](@ref)
"""
function map_sum_ξ_GNC_multipole(
     cosmo::Cosmology,
     ss=nothing;
     N_log::Int=1000,
     sum_xi::Bool=true,
     kwargs...)

     v_ss = isnothing(ss) ? 10 .^ range(0, log10(2 * cosmo.s_max), length=N_log) : ss

     ALL = [
          begin
               _, xis = map_ξ_GNC_multipole(cosmo, effect, v_ss;
                    sum_xi=sum_xi, specif_kwargs_GNC(effect, kwargs)...)
               xis
          end for effect in GaPSE.GR_EFFECTS_GNC
     ]

     return v_ss, sum(ALL), ALL
end



##########################################################################################92



"""
     print_map_sum_ξ_GNC_multipole(
          cosmo::Cosmology, out::String, ss = nothing;
          s1 = nothing, L::Int = 0, alg::Symbol = :lobatto,
          obs::Union{Bool,Symbol} = :noobsvel,
          N_lob::Int = 100, N_trap::Int = 50,
          atol_quad::Float64 = 0.0, rtol_quad::Float64 = 1e-2,
          enhancer::Float64=1e6, N_log::Int = 1000, 
          pr::Bool = true,
          single::Bool = true,
          kwargs...) 

Evaluate the multipole of order `L` of all the Galaxy Number Counts (GNC) 
Two-Point Correlation Function (TPCF) multipoles and their sum
in the comoving distance `s1`,  for all the comoving distances stored inside `ss` (representing 
the comoving distance from `s1`) for the input `cosmo::Cosmology`; finally,
it saves the results inside the file `out`.
If `ss = nothing`, it is set `ss = 10 .^ range(0, log10(2 * cosmo.s_max), length=N_log)`.
If `s1 = nothing`, it is set `s1 = cosmo.s_eff`.
We remember that all the distances are measured in ``h_0^{-1}\\mathrm{Mpc}``.

This function makes a for-loop on the `GaPSE.GR_EFFECTS_GNC` strings, calling 
`map_ξ_GNC_multipole` for each of them. To each string corresponds pretty intuitively one of the 
25 GNC effects. They are currently, in order:

`$(string(GR_EFFECTS_GNC .* " , "...))`

## Inputs

- `cosmo::Cosmology` : cosmology to be used in this computation

- `out::String` : name of the file where the results must be stored.

- `ss` : vector/range of `s` values where the function must be evaluated; if `ss = nothing`, 
  it is set `ss = 10 .^ range(0, log10(2 * cosmo.s_max), length=N_log)`.

## Optional arguments

This function recall internally `map_ξ_GNC_multipole`, so the kwargs are the same; we report them
for comfortness:

- `s1 = nothing` : comoving distance from the observer where the TPCF should be evaluated;
  if `s1 = nothing`, it is automatically set `s1 = cosmo.s_eff` from the given input `cosmo::Cosmology`.

- `L::Int = 0`: order of the Legendre polynomial to be used

- `alg::Symbol = :trap` : algorithm to be used for the integration; the valid options 
  are (other values will lead to `AssertionError`):
  - `:quad` -> the integration over ``\\mu`` will be preformed through the Julia function `quadgk` 
  from the [`QuadGK.jl`](https://github.com/JuliaMath/QuadGK.jl) Julia package, that uses an adaptive 
  Gauss-Kronrod quadrature.
  - `:trap` -> the integration over ``\\mu`` will be preformed through the Julia function `trapz` 
  from the [`Trapz.jl`](https://github.com/francescoalemanno/Trapz.jl) Julia package, that uses the
  simple trapezoidal rulae.
  - `:lobatto` -> the integration over ``\\mu`` will be preformed through the Julia function `gausslobatto` 
  from the [`FastGaussQuadrature.jl`](https://github.com/JuliaApproximation/FastGaussQuadrature.jl) Julia package, 
  that uses the Gauss-Lobatto quadrature. 
  WE RECOMMEND TO USE `:quad` FOR MONOPOLES AND `:lobatto` FOR HIGHER ORDER MULTIPOLES!
  
- `use_windows::Bool = false`: tells if the integrand must consider the two
   window function ``\\phi`` and ``\\mathcal{F}``

- `obs::Union{Bool,Symbol} = :noobsvel` : do you want to consider the observer terms in the computation of the 
  chosen GNC TPCF effect?
  - `:yes` or `true` -> all the observer effects will be considered
  - `:no` or `false` -> no observer term will be taken into account
  - `:noobsvel` -> the observer terms related to the observer velocity (that you can find in the CF concerning Doppler)
    will be neglected, the other ones will be taken into account

- `N_lob::Int = 100` : number of points to be used in the sampling made by the function `trapz`.
  Note that these options will have an effect only if you se `alg = :quad`.

- `N_trap::Int = 200` : number of points to be used in the sampling made by the function `trapz`.
  Note that these options will have an effect only if you se `alg = :quad`.

- `atol_quad::Float64 = 0.0` and `rtol_quad::Float64 = 1e-2`: absolute and relative tolerance
  to be passed to the function `quadgk`; it's recommended not to set `rtol_quad < 1e-2` 
  because the time for evaluation increase quickly.
  Note that these options will have an effect only if you se `alg = :quad`.

- `enhancer::Float64 = 1e6`: just a float number used in order to deal better with small numbers; 
  the returned value is NOT modified by this value, because after a multiplication
  the internal result is divided by `enhancer`.

- `N_log::Int = 1000` : number of points to be used in the default logaritmically-spaced 
  range for `ss`, i.e. `range(0, log10(2 * cosmo.s_max), length=N_log)`; it is ignored if `ss ≠ nothing` 

- `pr::Bool = true` : do you want the progress bar showed on screen, in order to 
  check the time needed for the computation? (`true` recommended)

- `single::Bool = true` : if `true`, all the CFs are printed inside the file of the sum, in a 
  table with 18 columns (first one for `ss`, second for their sum an the next 16 for each effect).
  Otherwise, a new directory "all_standalones_CFs" is created (in the same path given in `out`) and 
  they are separately saved in files there placed.

- `kwargs...` : other keyword arguments that will be passed to ALL the 
  GNC TPCF effect (`ξ_GNC_Doppler`, `ξ_GNC_Lensing`, ...); if one of them has that keyword argument,
  it will use the given value, otherwise it will be unaffected.


See also: [`map_ξ_GNC_multipole`](@ref),
[`sum_ξ_GNC_multipole`](@ref), [`map_sum_ξ_GNC_multipole`](@ref),
[`Cosmology`](@ref), [`GR_EFFECTS_GNC`](@ref)
"""
function print_map_sum_ξ_GNC_multipole(
     cosmo::Cosmology, out::String, v_ss = nothing;
     s1 = nothing, L::Int = 0,
     single::Bool = true,
     pr::Bool = true,
     sum_xi::Bool = true,
     kwargs...)

     check_parent_directory(out)
     check_namefile(out)

     dir = length(split(out, "/")) == 1 ? "all_GNC_standalones_CFs/" :
           join(split(out, "/")[begin:end-1] .* "/") * "all_GNC_standalones_CFs/"

     if single == false
          isdir(dir) || mkdir(dir)
     end

     s_1 = isnothing(s1) ? cosmo.s_eff : s1
     t1 = time()
     ss, xis, ALL = map_sum_ξ_GNC_multipole(cosmo, v_ss;
          L=L, s1=s_1, pr=pr, sum_xi=sum_xi, kwargs...)
     t2 = time()

     (pr) && println("\ntime needed for map_sum_ξ_GNC_multipole L=$L" *
                     "[in s] = $(@sprintf("%.5f", t2-t1)) \n")

     isfile(out) && run(`rm $out`)
     open(out, "w") do io
          println(io, BRAND)

          println(io, "#\n# This is an integration map on mu of the sum" *
                      " of all the ξ_GNC L=$L multipole GR effects")
          println(io, "# concerning the relativistic galaxy number counts.")
          !(single == true) ||
               println(io, "# In input was set \"single = $single\", " *
                           "so, together with their sum, all the CFs are here reported.\n#")
          !(single == false) ||
               println(io, "# In input was set \"single = $single\", so here is showed only \n" *
                           "# the sum of all the CFs, and them are saved separately in the \n" *
                           "# following directory: $dir \n#")

          parameters_used(io, cosmo; logo=false)
          println(io, "# computational time needed (in s) : $(@sprintf("%.4f", t2-t1))")
          print(io, "# kwards passed: ")

          println(io, "\n# \t\tL = $L")
          if !isempty(kwargs)
               for key in keys(kwargs)
                    println(io, "# \t\t$(key) = $(kwargs[key])")
               end
          end
          isnothing(s1) ||
               println(io, "#\n# NOTE: the computation is done not in s1 = s_eff, \n" *
                           "#\t because you specified in input s1 = $s1 Mpc/h_0!")
          println(io, "# ")

          if (single == false)
               println(io, "# s [Mpc/h_0] \t \t xi_SUM")
               for (s, xi) in zip(ss, xis)
                    println(io, "$s \t $xi")
               end

          elseif single == true
               println(io, "# 1: s [Mpc/h_0] \t 2: xi_SUM \t " *
                           join([string(i) for i in 3:length(GR_EFFECTS_GNC)+2] .*
                                ": xi_" .* GR_EFFECTS_GNC .* " \t "))
               for (i, (s, xi)) in enumerate(zip(ss, xis))
                    println(io, "$s \t $xi \t " *
                                join(["$(v[i]) \t " for v in ALL]))
               end
          else
               throw(ErrorException("how did you arrive here?"))
          end
     end

     if single == false
          for (effect, vec) in zip(GaPSE.GR_EFFECTS_GNC, ALL)
               open(dir * "xi_GNC_" * effect * "_L$L" * ".txt", "w") do io
                    println(io, "# This is an integration map on mu of the ξ multipole $effect GR effect")
                    println(io, "# concerning the relativistic galaxy number counts.")
                    parameters_used(io, cosmo)
                    #println(io, "# computational time needed (in s) : $(@sprintf("%.4f", t2-t1))")
                    print(io, "# kwards passed: ")

                    println(io, "\n# \t\tL = $L")
                    if !isempty(kwargs)
                         for key in keys(kwargs)
                              println(io, "# \t\t$(key) = $(kwargs[key])")
                         end
                    end
                    isnothing(s1) || println(io, "#\n# NOTE: the computation is done not in " *
                                                 "s1 = s_eff, because you specified in input s1 = $s1 !")
                    println(io, "# ")
                    println(io, "# s [Mpc/h_0] \t \t xi")
                    for (s, xi) in zip(ss, vec)
                         println(io, "$s \t $xi")
                    end

               end
          end
     end
end

