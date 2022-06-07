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


function FFTLog_PS_multipole(ss, xis;
     pr::Bool=true,
     L::Int=0, ν::Float64=1.5, n_extrap_low::Int=500,
     n_extrap_high::Int=500, n_pad::Int=500)

     XIS = xis .* ss .^ 3

     t1 = time()
     plan = FFTLog.SingleBesselPlan(x=ss,
          n_extrap_low=n_extrap_low, ν=ν,
          n_extrap_high=n_extrap_high, n_pad=n_pad)
     FFTLog.prepare_FFTLog!(plan, [L])
     pks = reshape(FFTLog.evaluate_FFTLog(plan, XIS), (:,))
     ks = reshape(FFTLog.get_y(plan), (:,))
     t2 = time()
     pr && println("\ntime needed for this Power Spectrum computation [in s] = $(t2-t1)\n")


     if iseven(L)
          return ks, (1 / A_prime * (-1)^(L / 2)) .* pks
     else
          return ks, (1 / A_prime * (-im)^L) .* pks
     end
end



function FFTLog_all_PS_multipole(input::String,
     group::String=VALID_GROUPS[end];
     L::Int=0, pr::Bool=true, ν::Float64=1.5, n_extrap_low::Int=500,
     n_extrap_high::Int=500, n_pad::Int=500)

     check_group(group; valid_groups=VALID_GROUPS)
     check_fileisingroup(input, group)

     pr && begin
          print("\nI'm computing the PS_multipole from the file \"$input\"")
          if group == "GNC"
               println("for the Galaxy Number Counts.")
          elseif group == "LD"
               println("for the Luminosity Distance perturbations.")
          elseif group == "GNCxLD"
               println("for the cross correlations between " *
                       "Galaxy Number Counts and Luminosity Distance perturbations.")
          elseif group == "LDxGNC"
               println("for the cross correlations between " *
                       "Luminosity Distance perturbations and Galaxy Number Counts.")
          else
               println("(no specific group considered).")
          end
     end

     sps = group == "GNC" ? "GNC GR effects" :
           group == "LD" ? "LD GR effects" :
           group == "GNCxLD" ? "GNCxLD GR effects" :
           group == "LDxGNC" ? "LDxGNC GR effects" :
           "generic file"

     table = readdlm(input; comments=true)
     ss = convert(Vector{Float64}, table[:, 1])
     all_xis = [convert(Vector{Float64}, col)
                for col in eachcol(table[:, 2:end])]

     ALL_XIS = [xis .* ss .^ 3 for xis in all_xis]

     plan = FFTLog.SingleBesselPlan(x=ss,
          n_extrap_low=n_extrap_low, ν=ν,
          n_extrap_high=n_extrap_high, n_pad=n_pad)
     FFTLog.prepare_FFTLog!(plan, [L])
     ALL_pks = pr ? begin
          @showprogress sps * ", L=$L: " [reshape(FFTLog.evaluate_FFTLog(plan, XIS), (:,))
               for XIS in ALL_XIS] 
          end : begin
               [reshape(FFTLog.evaluate_FFTLog(plan, XIS), (:,))
               for XIS in ALL_XIS] 
          end

     ks = reshape(FFTLog.get_y(plan), (:,))

     if iseven(L)
          return ks, [(1 / A_prime * (-1)^(L / 2)) .* pks for pks in ALL_pks]
     else
          return ks, [(1 / A_prime * (-im)^L) .* pks for pks in ALL_pks]
     end
end

