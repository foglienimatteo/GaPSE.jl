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



function integrand_ξ_PPMatter_multipole(s, μ, cosmo::Cosmology;
     L::Int=0, use_windows::Bool=true)

     res = if use_windows == true
          #println("s1 = $s1 \\t s2 = $(s2(s1, s, μ)) \\t  y=$(y(s1, s, μ))")
          int = cosmo.ξ_matter(s)
          #println("int = $int")
          coef = (cosmo.params.b + cosmo.f_of_s(s) * μ^2)^2 * cosmo.D_of_s(s)
          #println("coef = $coef")
          int .* (coef * spline_integrF(s, μ, cosmo.windowFint) / cosmo.WFI_norm * Pl(μ, L))
     else
          int = cosmo.ξ_matter(s)
          coef = (cosmo.params.b + cosmo.f_of_s(s) * μ^2)^2
          int .* coef .* Pl(μ, L)
     end

     #println("res = $res")
     return (2.0 * L + 1.0) / 2.0 * res
end



##########################################################################################92



function ξ_PPMatter_multipole(
     s, cosmo::Cosmology;
     L::Int=0,
     use_windows::Bool=true,
     enhancer::Float64=1e6,
     N_μs::Int=200,
     μ_atol::Float64=0.0,
     μ_rtol::Float64=1e-2,
     trap::Bool=false)

     orig_f(μ) = enhancer * integrand_ξ_PPMatter_multipole(s, μ, cosmo;
          L=L, use_windows=use_windows)

     int =
          if s > 1.0 && trap == false
               quadgk(μ -> orig_f(μ), -1.0, 1.0; atol=μ_atol, rtol=μ_rtol)[1]
          else
               #=
               μs = union(
                    range(-1.0, -0.90, length=N_μs),
                    range(-0.90, 0.90, length=N_μs),
                    range(0.90, 1.0, length=N_μs)
                    )
               =#
               μs = range(-1.0 + 1e-6, 1.0 - 1e-6, length=N_μs)
               orig_fs = orig_f.(μs)
               #orig_fs
               trapz(μs, orig_fs)

          end

     return int / enhancer
end


##########################################################################################92



function map_ξ_PPMatter_multipole(cosmo::Cosmology,
     v_ss=nothing;
     pr::Bool=true,
     N_log::Int=1000,
     L::Int=0,
     kwargs...)

     t1 = time()
     ss = isnothing(v_ss) ? 10 .^ range(0, 3, length=N_log) : v_ss
     xis = pr ? begin
          @showprogress "PP Matter, L=$L: " [
               ξ_PPMatter_multipole(s, cosmo; L=L, kwargs...) for s in ss
          ]
     end : [
          ξ_PPMatter_multipole(s, cosmo; L=L, kwargs...) for s in ss
     ]

     t2 = time()
     pr && println("\ntime needed for map_ξ_PPMatter_multipole " *
                   "[in s] = $(@sprintf("%.5f", t2-t1)) ")
     return (ss, xis)
end


##########################################################################################92



function print_map_ξ_PPMatter_multipole(
     cosmo::Cosmology,
     out::String,
     v_ss=nothing;
     kwargs...)

     t1 = time()
     vec = map_ξ_PPMatter_multipole(cosmo, v_ss; kwargs...)
     t2 = time()

     isfile(out) && run(`rm $out`)
     open(out, "w") do io
          println(io, "# This is an integration map on mu of the ξ multipole concerning the PP Matter.")
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

          println(io, "# ")
          println(io, "# s [Mpc/h_0] \t \t xi")
          for (s, xi) in zip(vec[1], vec[2])
               println(io, "$s \t $xi")
          end
     end
end




