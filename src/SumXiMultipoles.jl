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



function sum_ξ_multipole(s1, s, cosmo::Cosmology;
     all::Bool = false,
     kwargs...)

     ALL = [ξ_multipole(s1, s, effect, cosmo; kwargs...)
            for effect in IMPLEMENTED_GR_EFFECTS]

     if all
          return sum(ALL), ALL
     else
          return sum(ALL)
     end
end



function map_sum_ξ_multipole(
     cosmo::Cosmology,
     v_ss::Union{Vector{Float64},Nothing} = nothing;
     N_log::Integer = 1000,
     all::Bool = false,
     kwargs...)

     ss = isnothing(v_ss) ? 10 .^ range(-1, 3, length = N_log) : v_ss


     ALL = [
          begin
               _, xis = map_ξ_multipole(cosmo, effect, ss; kwargs...)
               xis
          end for effect in IMPLEMENTED_GR_EFFECTS
     ]

     if all
          return ss, sum(ALL), ALL
     else
          return ss, sum(ALL)
     end
end


function print_map_sum_ξ_multipole(
     cosmo::Cosmology,
     out::String,
     v_ss::Union{Vector{Float64},Nothing} = nothing;
     s_1::Union{Float64,Nothing} = nothing,
     L::Integer = 0,
     all::Bool = true,
     single::Bool = true,
     kwargs...)

     dir = length(split(out, "/")) == 1 ? "all_standalones_CFs/" :
           join(split(out, "/")[begin:end-1] .* "/") * "all_standalones_CFs/"

     if single == false
          isdir(dir) || mkdir(dir)
     end

     s1 = isnothing(s_1) ? cosmo.s_eff : s_1
     t1 = time()
     ss, xis, ALL = map_sum_ξ_multipole(cosmo, v_ss;
          L = L, s_1 = s1, all = true, kwargs...)
     t2 = time()

     isfile(out) && run(`rm $out`)
     open(out, "w") do io
          println(io, "# This is an integration map on mu of the sum" *
                      " of all the ξ_L=$L multipole GR effects.")
          !(all == true && single == true) ||
               println(io, "# In input was set \"all = $all\" and \"single = $single\" \n" *
                           "# so, together with their sum, all the CFs are here reported.\n#")
          !(all == true && single == false) ||
               println(io, "# In input was set \"all = $all\" and \"single = $single\" \n" *
                           "# so here is showed only the sum of all the CFs, and them are\n" *
                           "# saved separately in the directory: $dir \n#")
          !(all == false) ||
               println(io, "# In input was set \"all = $all\"" *
                           "# so only the sum of all the CFs was saved in a file (this one).\n#")

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

          if (all == false) || (all == true && single == false)
               println(io, "# s [Mpc/h_0] \t \t xi_SUM")
               for (s, xi) in zip(ss, xis)
                    println(io, "$s \t $xi")
               end

          elseif all == true && single == true
               println(io, "# s [Mpc/h_0] \t xi_SUM \t " *
                           join("xi_" .* GaPSE.IMPLEMENTED_GR_EFFECTS .* " \t "))
               for (i, (s, xi)) in enumerate(zip(ss, xis))
                    println(io, "$s \t $xi \t " *
                                join(["$(v[i]) \t " for v in ALL]))
               end
          else
               throw(ErrorException("how did you arrive here?"))
          end
     end

     if all == true && single == false
          for (effect, vec) in zip(IMPLEMENTED_GR_EFFECTS, ALL)
               open(dir * "xi_" * effect * "_L$L" * ".txt", "w") do io
                    println(io, "# This is an integration map on mu of the ξ multipole $effect GR effect.")
                    parameters_used(io, cosmo)
                    #println(io, "# computational time needed (in s) : $(@sprintf("%.4f", t2-t1))")
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
                    for (s, xi) in zip(ss, vec)
                         println(io, "$s \t $xi")
                    end

               end
          end
     end
end

