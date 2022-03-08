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

function PS_multipole(
     f_in::Union{Function,Dierckx.Spline1D};
     int_s_min::Float64 = 1e-1, int_s_max::Float64 = 1e3,
     L::Integer = 0, N::Integer = 1024,
     pr::Bool = true, kwargs...)

     t1 = time()
     ks, pks = xicalc(s -> 2 * π^2 * f_in(s), L, 0;
          N = N, kmin = int_s_min, kmax = int_s_max, r0 = 1 / int_s_max)
     t2 = time()
     pr && println("\ntime needed for Power Spectrum  computation [in s] = $(t2-t1)\n")

     if iseven(L)
          return ks, (1 / π * (-1)^(L / 2)) .* pks #(1 / A_prime * (-1)^(L / 2)) .* pks
     else
          return ks, (1 / π * (-im)^L) .* pks #(1 / A_prime * (-im)^L) .* pks
     end
end


function PS_multipole(
     input::String;
     L::Integer = 0, N::Integer = 1024,
     pr::Bool = true,
     int_s_min::Union{Float64,Nothing} = nothing,
     int_s_max::Union{Float64,Nothing} = nothing)

     xi_table = readdlm(input, comments = true)
     ss = convert(Vector{Float64}, xi_table[:, 1])
     fs = convert(Vector{Float64}, xi_table[:, 2])

     f_in = Spline1D(ss, fs; bc = "error")

     intsmin = isnothing(int_s_min) ? min(ss...) : int_s_min
     intsmax = isnothing(int_s_max) ? max(ss...) : int_s_max
     return PS_multipole(f_in; int_s_min = intsmin, int_s_max = intsmax, N = N, L = L, pr = pr)
end



function PS_multipole(
     effect::String, cosmo::Cosmology;
     int_s_min::Float64 = 1e-1, int_s_max::Float64 = 1e3,
     L::Integer = 0, N::Integer = 100,
     pr::Bool = true, kwargs...)

     error = "$effect is not a valid GR effect name.\n" *
             "Valid GR effect names are the following:\n" *
             string(IMPLEMENTED_GR_EFFECTS .* " , "...)
     @assert (effect ∈ IMPLEMENTED_GR_EFFECTS) error

     xs, ys = map_ξ_multipole(cosmo, effect; L = L, pr = pr, kwargs...)
     f_in = Spline1D(xs, ys; bc = "error")

     t1 = time()
     ks, pks = xicalc(s -> 2 * π^2 * f_in(s), L, 0;
          N = N, kmin = int_s_min, kmax = int_s_max, r0 = 1 / int_s_max)
     t2 = time()
     pr && println("\ntime needed for Power Spectrum  computation [in s] = $(t2-t1)\n")

     if iseven(L)
          return ks, (1 / π * (-1)^(L / 2)) .* pks #(1 / A_prime * (-1)^(L / 2)) .* pks
     else
          return ks, (1 / π * (-im)^L) .* pks #(1 / A_prime * (-im)^L) .* pks
     end
end


@doc raw"""
     PS_multipole(in::String, int_s_min = 0.0, int_s_max = 1000.0; 
          L::Integer = 0, N::Integer = 1024, pr::Bool = true)
     PS_multipole(
          f_in::Union{Function,Dierckx.Spline1D}, int_s_min = 0.0, int_s_max = 1000.0; 
          L::Integer = 0, N::Integer = 1024, pr::Bool = true
          ) :: Tuple{Vector{Float64}, Vector{Float64}}

Return the `L`-order multipole from the input function `f_in`, through the
following Fast Fourier Transform and the effective redshift approximation:

```math
P_L(k) = \frac{2 L + 1}{A^{'}} (-i)^L \, \phi(s_\mathrm{eff}) \int_0^\infty 
        \mathrm{d} s \; s^2 \, j_L(ks) \, f_\mathrm{in}(s) \; ,
        \quad \; A^{'} = 8 \, \pi^2
```

This expression can be easily obtained from the standard one:
```math
\begin{split}
    P_L(k) = &\frac{2 L + 1}{A} (-i)^L \, 
        \int_0^\infty \mathrm{d} s_1 \; s_1^2 
        \int_0^\infty \mathrm{d} s \; s^2 
        \int_{-1}^{+1} \mathrm{d} \mu \;
        j_L(ks) \, \xi(s_1, s, \mu) \, \phi(s_1) \, \phi(s_2) \,
        \mathcal{L}_L(\mu) F\left(\frac{s}{s_1}, \mu \right) \\
        &\mathrm{with} \; \;s_2 = s_2(s_1, s, μ) = \sqrt{s_1^2 + s^2 + 2s_1s\mu}
        \; 
        , \quad A(s_\mathrm{max}, s_\mathrm{min}, \theta_\mathrm{max}) 
        = 2 \, \pi \, V
\end{split}
```
with the definition
```math
f_\mathrm{in}(s_1, s) =  \int_{-1}^{+1} \mathrm{d} \mu \;
        \xi(s_1, s, \mu) \, \phi(s_2) \,
        \mathcal{L}_L(\mu) \, F\left(\frac{s}{s_1}, \mu \right)
```
and the application of the effective redshift approximation.

See also: [`V_survey`](@ref), [`A`](@ref), [`A_prime`](@ref)
"""
PS_multipole



##########################################################################################92



function print_PS_multipole(in::String, out::String;
     L::Integer = 0, N::Integer = 1024,
     pr::Bool = true, kwargs...)

     pr && println("\nI'm computiong the PS_multipole from the file $in")

     time_1 = time()
     vec = PS_multipole(in; N = N, L = L, pr = pr, kwargs...)
     time_2 = time()

     pr && println("\ntime needed for Power Spectrum  computation [in s] = $(time_2-time_1)\n")

     isfile(out) && run(`rm $out`)
     open(out, "w") do io
          println(io, "# Power Spectrum Multipole computation of the file: $in")
          println(io, "#\n# For this PS_multipole computation we set: ")
          println(io, "# \t #points used in Fourier transform N = $N")
          println(io, "# \t multipole degree in consideration L = $L")
          println(io, "# computational time needed (in s) : $(@sprintf("%.4f", time_2-time_1))")
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
          println(io, "# k [h_0/Mpc] \t \t  P [(Mpc/h_0)^3]")
          for (k, pk) in zip(vec[1], vec[2])
               println(io, "$k \t $pk")
          end
     end
end


function print_PS_multipole(
     effect::String, out::String,
     cosmo::Cosmology;
     L::Integer = 0, N::Integer = 1024,
     pr::Bool = true, kwargs...)

     error = "$effect is not a valid GR effect name.\n" *
             "Valid GR effect names are the following:\n" *
             string(IMPLEMENTED_GR_EFFECTS .* " , "...)
     @assert (effect ∈ IMPLEMENTED_GR_EFFECTS) error

     pr && println("\nI'm computiong the PS_multipole for the $effect GR effect.")

     time_1 = time()
     vec = PS_multipole(effect, cosmo; N = N, L = L, pr = pr, kwargs...)
     time_2 = time()

     pr && println("\ntime needed for Power Spectrum  computation [in s] = $(time_2-time_1)\n")

     isfile(out) && run(`rm $out`)
     open(out, "w") do io
          println(io, "# Power Spectrum Multipole computation for the $effect GR effect.")
          println(io, "# The following cosmology data were given:\n#")
          parameters_used(io, cosmo)
          println(io, "#\n# For this PS_multipole computation we set: ")
          println(io, "# \t #points used in Fourier transform N = $N")
          println(io, "# \t multipole degree in consideration L = $L")
          println(io, "# computational time needed (in s) : $(@sprintf("%.4f", time_2-time_1))")
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
          println(io, "# k [h_0/Mpc] \t \t  P [(Mpc/h_0)^3]")
          for (k, pk) in zip(vec[1], vec[2])
               println(io, "$k \t $pk")
          end
     end
end