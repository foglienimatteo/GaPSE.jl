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

function integral_on_mu(s1, s, integrand::Function, cosmo::Cosmology;
     L::Integer = 0,
     enhancer::Float64 = 1e6,
     use_windows::Bool = true,
     μ_atol::Float64 = 1e-4,
     μ_rtol::Float64 = 1e-1,
     kwargs...
)

     f(μ) = integrand(s1, s, μ, cosmo; enhancer = enhancer, L = L,
          use_windows = use_windows, kwargs...)[1]

     #println("s1 = $s1 \t s = $s")
     int = quadgk(μ -> f(μ), -1.0, 1.0; rtol = μ_rtol, atol = μ_atol)
     #println("s1 = $s1 \t s2 = $s \t int = $int")
     return int ./ enhancer
end



function integral_on_mu(s1, s, effect::String, cosmo::Cosmology; kwargs...)
     error = "$effect is not a valid GR effect name.\n" *
             "Valid GR effect names are the following:\n" *
             string(keys(dict_gr_mu) .* " , "...)

     @assert (effect ∈ keys(dict_gr_mu)) error
     integral_on_mu(s1, s, dict_gr_mu[effect], cosmo; kwargs...)
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



function my_integral_on_mu(s1, s, integrand, cosmo::GaPSE.Cosmology;
     L::Integer = 0,
     enhancer::Float64 = 1e6,
     use_windows::Bool = true,
     μ_steps = 30,
     kwargs...
)

     μs1 = range(-1.0, -0.90, length = μ_steps)
     μs2 = range(-0.90, 0.90, length = μ_steps)
     μs3 = range(0.90, 1.0, length = μ_steps)
     μs = unique(vcat(μs1, μs2, μs3))
     fs = [integrand(s1, s, μ, cosmo; enhancer = enhancer, L = L,
          use_windows = use_windows, kwargs...) for μ in μs]

     #println("s1 = $s1 \t s = $s")
     #int = quadgk(μ -> f(μ), -1.0, 1.0; rtol = μ_rtol, atol = μ_atol)
     #println("s1 = $s1 \t s2 = $s \t int = $int")
     return trapz(μs, fs) / enhancer
end


function my_integral_on_mu(s1, s, effect::String, cosmo::GaPSE.Cosmology; kwargs...)
     error = "$effect is not a valid GR effect name.\n" *
             "Valid GR effect names are the following:\n" *
             string(keys(dict_gr_mu) .* " , "...)

     @assert (effect ∈ keys(dict_gr_mu)) error
     my_integral_on_mu(s1, s, dict_gr_mu[effect], cosmo; kwargs...)
end

##########################################################################################92


function map_integral_on_mu(
     cosmo::Cosmology,
     effect::Union{String,Function},
     v_ss::Union{Vector{Float64},Nothing} = nothing;
     s_1::Union{Float64,Nothing} = nothing,
     use_my::Bool = true,
     pr::Bool = true, enhancer = 1e6, kwargs...)

     s1 = isnothing(s_1) ? cosmo.s_eff : s_1

     t1 = time()
     ss = isnothing(v_ss) ? 10 .^ range(-1, 3, length = 100) : v_ss
     xis = if use_my == true
          println("I will use trapz.")
          my_f(s) = my_integral_on_mu(s1, s, effect, cosmo; enhancer = enhancer, kwargs...)
          omg = @showprogress [my_f(s) for s in ss]
          omg
     else
          println("I will use quadgk.")
          f(s) = integral_on_mu(s1, s, effect, cosmo; enhancer = enhancer, kwargs...)
          vec = @showprogress [f(s) for s in ss]
          omg, xis_err = [x[1] for x in vec], [x[2] for x in vec]
          omg
     end

     t2 = time()
     pr && println("\ntime needed for map_integral_on_mu for $effect [in s] = $(t2-t1)")
     return (ss, xis)
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


##########################################################################################92



function ξ_multipole(s1, s, effect::Function, cosmo::Cosmology; L::Integer = 0, kwargs...)
     error = "$(string(effect)) is not a valid GR effect function.\n" *
             "Valid GR effect functions are the following:\n" *
             string(values(dict_gr_mu) .* " , "...)
     @assert (effect ∈ values(dict_gr_mu)) error

     return (2.0 * L + 1.0) / 2.0 .* integral_on_mu(s1, s, effect, cosmo; L = L, kwargs...)
end


function ξ_multipole(s1, s, effect::String, cosmo::Cosmology; L::Integer = 0, kwargs...)
     error = "$effect is not a valid GR effect name.\n" *
             "Valid GR effect names are the following:\n" *
             string(keys(dict_gr_mu) .* " , "...)
     @assert (effect ∈ keys(dict_gr_mu)) error

     return (2.0 * L + 1.0) / 2.0 .* integral_on_mu(s1, s, dict_gr_mu[effect], cosmo; L = L, kwargs...)
end



##########################################################################################92


#=
function my_integral_on_mu(s1, s, integrand, cosmo::GaPSE.Cosmology;
        L::Integer = 0,
        enhancer::Float64 = 1e6,
        use_windows::Bool = true,
        #μ_atol::Float64 = 1e-3,
        #μ_rtol::Float64 = 1e-3,
        Numb=10,
        kwargs...
        )
    
    

    μs1 = range(-1.0,-0.9, length=Numb)
    μs2 = range(-0.9, 0.9, length=Numb)
    μs3 = range(0.9, 1.0, length=Numb)
    μs = unique(vcat(μs1, μs2, μs3))
    fs = [integrand(s1, s, μ, cosmo; enhancer = enhancer, L = L,
        use_windows = use_windows, kwargs...) for μ in μs ]
    
    

    #println("s1 = $s1 \t s = $s")
    #int = quadgk(μ -> f(μ), -1.0, 1.0; rtol = μ_rtol, atol = μ_atol)
    #println("s1 = $s1 \t s2 = $s \t int = $int")
    return trapz(μs,fs) / enhancer
end
=#

