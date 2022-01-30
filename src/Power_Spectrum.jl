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


IMPLEMENTED_GR_EFFECTS = ["auto_doppler", "auto_lensing"]
IMPLEMENTED_XI_FUNCS = [integral_on_mu_doppler, integral_on_mu_lensing]
dict_gr_xi = Dict([a => b for (a, b) in zip(IMPLEMENTED_GR_EFFECTS, IMPLEMENTED_XI_FUNCS)]...)

function PS_multipole(f_in::Union{Function,Dierckx.Spline1D};
     L::Integer = 0, N = 128,
     int_s_min = 1e-2, int_s_max = 2 * s_max,
     pr::Bool = true)

     t1 = time()
     ks, pks = xicalc(x -> 2 * π^2 * f_in(x), L, 0; N = N, kmin = int_s_min, kmax = int_s_max, r0 = 1 / int_s_max)
     t2 = time()
     pr && println("\ntime needed for Power Spectrum  computation [in s] = $(t2-t1)\n")

     if iseven(L)
          return ks, ((2 * L + 1) / A_prime * ϕ(s_eff) * (-1)^(L / 2)) .* pks
     else
          return ks, ((2 * L + 1) / A_prime * ϕ(s_eff) * (-im)^L) .* pks
     end
end


function PS_multipole(in::String; int_s_min = 1e-2, int_s_max = 2 * s_max, N = 128,
     L::Integer = 0, pr::Bool = true)

     xi_table = readdlm(in, comments = true)
     ss = convert(Vector{Float64}, xi_table[2:end, 1])
     fs = convert(Vector{Float64}, xi_table[2:end, 2])
     f_in = Spline1D(ss, fs)

     return PS_multipole(f_in; int_s_min = int_s_min, int_s_max = int_s_max, N = N, L = L, pr = pr)
end



@doc raw"""
     PS_multipole(in::String; int_s_min = 1e-2, int_s_max = 2 * s_max, N = 128,
          L::Integer = 0, pr::Bool = true)
     PS_multipole(f_in::Union{Function,Dierckx.Spline1D};
          L::Integer = 0,  N = 128, 
          int_s_min = 1e-2, int_s_max = 2 * s_max,
          pr::Bool = true) :: Tuple{Vector{Float64}, Vector{Float64}}

Return the `L`-order multipole from the input function `f_in`, through the
following Fast Fourier Transform and the effective redshift approximation:

```math
P_L(k) = \frac{2 L + 1}{A^{'}} (-i)^L \, \phi(s_\mathrm{eff}) \int_0^\infty 
        \mathrm{d} s \; s^2 \, j_L(ks) \, f_\mathrm{in}(s) \; ,
        \quad \; A^{'} = 4 \, \pi^2
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



function print_PS_multipole(out::String, in::String;
     int_s_min = 1e-2, int_s_max = 2 * s_max, s1 = s_eff,
     N = 128, L::Integer = 0, pr::Bool = true, kwargs...)

     time_1 = time()

     vec = if isfile(in)
          pr && println("\nI'm computiong the PS_multipole from the file $in")
          PS_multipole(in, int_s_min = int_s_min, int_s_max = int_s_max,
               N = N, L = L, pr = pr, kwargs...)
     else
          if in ∈ IMPLEMENTED_GR_EFFECTS
               pr && println("\nI'm computiong the PS_multipole for the $in GR effect.")
               t1 = time()
               ss = 10 .^ range(-1, 3, length = 100)
               v = [dict_gr_xi[in](s1, s; L = L, kwargs...) for s in ss]
               xis, xis_err = [x[1] for x in v], [x[2] for x in v]
               t2 = time()
               pr && println("\ntime needed to create the xi map [in s] = $(t2-t1)\n")
               f_in = Spline1D(ss, xis)
               PS_multipole(f_in, int_s_min = int_s_min, int_s_max = int_s_max,
                    N = N, L = L, pr = pr, kwargs...)
          else
               throw(ErrorException(
                    "$in is neither a GR impemented effect or a file.\n" *
                    "\t The implemented GR effects are currently: \n" *
                    "\t $(IMPLEMENTED_GR_EFFECTS)"
               ))
          end
     end

     time_2 = time()

     isfile(out) && run(`rm $out`)
     open(out, "w") do io
          println(io, "# This is the Power Spectrum computation of ")
          parameters_used(io)
          println(io, "#\n# For this PS_multipole computation we set: ")
          println(io, "# \t int_s_min = $int_s_min \t int_s_max = $int_s_max ")
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

          println(io, "\nk \t P")
          for (k, pk) in zip(vec[1], vec[2])
               println(io, "$k \t $pk")
          end
     end
end
