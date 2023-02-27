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
     GenericWindow(
          comdist::Vector{Float64}
          multipoles::Vector{Vector{Float64}}
          splines::Vector{Dierckx.Spline1D}
     )

Stores the multipoles of a generic window function, computed as:
```math
     Q_{\\ell_1} = \\int_0^{\\infty} \\mathrm{d}s_1 \\, s_1^2 \\, \\phi(s_1) \\, F_{\\ell_1}(s_1,s)
```

with some FFT algorithm. See Eq. (2.13) of Castorina, Di Dio (2021) for more details.

## Constructors

     GenericWindow(file::String; comments::Bool=true, xdt::DataType=Float64, ydt::DataType=Float64)

Read the file `file` and create a `GenericWindow` struct having
- as first element the comoving distances stored in the first column (with the input type `xdt`)
- as second element a vector that contains all the following columns 
  (with the input type `ydt`), which are the multipoles `L=0,1,2,...` of the Window FUnction considered
- as second element a vector that contains all the splines of that multipoles

If the file start with comments (lines starting with #), set `comments = true`.
"""
struct GenericWindow
     comdist::Vector{Float64}
     multipoles::Vector{Vector{Float64}}
     splines::Vector{Dierckx.Spline1D}

     function GenericWindow(file::String; comments::Bool=true, xdt::DataType=Float64, ydt::DataType=Float64)

          table = readdlm(file, comments=comments)
          ss = vecstring_to_vecnumbers(table[:, 1]; dt=xdt)
          all = [vecstring_to_vecnumbers(col; dt=ydt)
                 for col in eachcol(table[:, 2:end])]

          for col in all
               @assert length(ss) == length(col) "the columns must have all the same length"
          end

          splines = [Spline1D(ss, col; bc="error") for col in all]

          new(ss, all, splines)
     end
end



"""
     create_file_for_XiMultipoles(out::String, names::Vector{String}, 
          effect::Union{String, Integer}, group::String="GNC"; 
          comments::Bool=true, xdt::DataType=Float64, ydt::DataType=Float64)

Read the column number `effect` (if is an Integer) or the one corresponding to the GR effect `effect`
for the input group `group` (if is a String) from all the filenames stored in the Vector `names`, 
and save them in a file named `out`.

The first column of `out` will be the same as the first column of the first filename in `names`; it
is however checked internally if the first column of all the other files coincides with this one.
The following columns of `out` follow the order in `names`. Note that `effect`, if passed as Integer,
must be > 1 (because 1 is the index of the first column, used as x-axis).

`group` must be one among the following: $(GaPSE.VALID_GROUPS)
If `group=$(GaPSE.VALID_GROUPS[end])`, then `effect` must be an integer (because 
you are not selecting a specific effect in one of the native GaPSE groups).

`xdt` and `ydt` are the data types to be used for respectively the first column and the 2-3-4-... columns 
of `out`.
Set `comments=true` if the files in `names` start with a header that must be skipped (its lines must start with
#, otherwise they will not be recognised as comments).

## Example

julia> run(`cat file_1.txt`)
# Generic comment line
# of the file_1.txt
1.0  0.999999  0.34545   0.00991
...  ...       ...       ...

julia> run(`cat file_2.txt`)
# same, for file_2.txt
1.0  0.58244  0.12123    0.000154
...  ...       ...       ...

julia> create_file_for_XiMultipoles("mix.txt", ["file_1.txt", "file_2.txt"], 3, "generic");
julia> run(`cat mix.txt`)
###############
#    GaPSE    #
############### 
#
#
# This is a table containing the multipoles of the Two-Point Correlation Function (TPCF) 
# for a generic group effect [not given, provied only the index 2] taken from the files:
#   - L = 0 : file_1.txt
#   - L = 1 : file_2.txt
#
# s [Mpc/h_0] 	xi_{L=0} 	 xi_{L=1} 		 
1.0   0.34545   0.12123 
...   ...       ...
"""
function create_file_for_XiMultipoles(out::String, names::Vector{String},
     effect::Union{String,Integer}, group::String="GNC";
     comments::Bool=true, xdt::DataType=Float64, ydt::DataType=Float64)

     @assert length(names) > 0 "at least one file name must be given!"
     @assert group ∈ VALID_GROUPS "group must be one among $(VALID_GROUPS); $group is not valid!"

     if typeof(effect) <: Integer
          @assert effect > 1 "if you pass effect as number, it must" *
                             " be an integer > 1, not $effect !"

     elseif group == "GNC"
          error = "$effect is not a valid GR effect name.\n" *
                  "Valid GR effect names are the following:\n" *
                  string(GaPSE.GR_EFFECTS_GNC .* " , "...)
          @assert (effect ∈ GaPSE.GR_EFFECTS_GNC) error

     elseif group == "LD"
          error = "$effect is not a valid GR effect name.\n" *
                  "Valid GR effect names are the following:\n" *
                  string(GaPSE.GR_EFFECTS_LD .* " , "...)
          @assert (effect ∈ GaPSE.GR_EFFECTS_LD) error

     elseif group == "GNCxLD"
          error = "$effect is not a valid GR effect name.\n" *
                  "Valid GR effect names are the following:\n" *
                  string(GaPSE.GR_EFFECTS_GNCxLD .* " , "...)
          @assert (effect ∈ GaPSE.GR_EFFECTS_GNCxLD) error

     elseif group == "LDxGNC"
          error = "$effect is not a valid GR effect name.\n" *
                  "Valid GR effect names are the following:\n" *
                  string(GaPSE.GR_EFFECTS_LDxGNC .* " , "...)
          @assert (effect ∈ GaPSE.GR_EFFECTS_LDxGNC) error

     elseif group == "generic"
          throw(AssertionError(
               """if you set group="generic", you should't give in input a string (because 
               you are not selecting a specific effect in a group) but an integer for effect."""
          ))
     else
          throw(AssertionError("how did you arrive here?"))
     end


     index = (typeof(effect) <: Integer || group == "generic") ? begin
          effect
     end : (group == "GNC") ? begin
          GaPSE.INDEX_GR_EFFECT_GNC[effect] + 2
     end : (group == "LD") ? begin
          GaPSE.INDEX_GR_EFFECT_LD[effect] + 2
     end : (group == "GNCxLD") ? begin
          GaPSE.INDEX_GR_EFFECT_GNCxLD[effect] + 2
     end : (group == "LDxGNC") ? begin
          GaPSE.INDEX_GR_EFFECT_LDxGNC[effect] + 2
     end : throw(AssertionError("how did you arrive here?"))

     table = readdlm(names[1], comments=comments)
     ss = GaPSE.vecstring_to_vecnumbers(table[:, 1]; dt=xdt)

     ALL = [
          begin
               ops = readdlm(name, comments=comments)
               here_ss = GaPSE.vecstring_to_vecnumbers(ops[:, 1]; dt=xdt)
               col = GaPSE.vecstring_to_vecnumbers(ops[:, index]; dt=ydt)
               @assert length(ss) == length(col) "the columns must have all the same length, $name differs!"
               @assert all(ss .≈ here_ss) "the ss must have all the same values, $name differs!"
               col
          end for name in names]

     open(out, "w") do io
          println(io, GaPSE.BRAND)
          println(io, "#\n# This is a table containing the multipoles of the Two-Point Correlation Function (TPCF) ")

          EFFECT = typeof(effect) == String ? effect : "[not given, provied only the index $effect]"
          if group == "GNC"
               println(io, "# for the Galaxy Number Counts GR effect $EFFECT" *
                           "\n#  taken from the files:")
          elseif group == "LD"
               println(io, "# for the Luminosity Distance perturbations GR effect $EFFECT" *
                           "\n# taken from the files:")
          elseif group == "GNCxLD"
               println(io, "# for the cross correlations between " *
                           "Galaxy Number Counts and Luminosity Distance perturbations \n#" *
                           "effect $EFFECT taken from the files:")
          elseif group == "LDxGNC"
               println(io, "# for the cross correlations between " *
                           "Luminosity Distance perturbations and Galaxy Number Counts \n#" *
                           "effect $EFFECT taken from the files:")
          elseif group == "generic"
               println(io, "# for a generic group " *
                           "effect $EFFECT taken from the files:")
          else
               throw(AssertionError("how did you arrive here?"))
          end

          for (i, name) in enumerate(names)
               println(io, "#   - L = $(i-1) : $name")
          end

          println(io, "#\n# s [Mpc/h_0] \t" .* join(["xi_{L=$(i-1)} \t " for i in 1:length(names)]))
          for (i, s) in enumerate(ss)
               println(io, "$s \t " *
                           join([" $(v[i]) \t " for v in ALL]))
          end
     end
end



"""
     XiMultipoles(
          comdist::Vector{Float64}
          multipoles::Vector{Vector{Float64}}
     )

Stores the multipoles of a generic Two-Point Correlation Function.

## Constructors

     XiMultipoles(file::String; comments::Bool=true, xdt::DataType=Float64, ydt::DataType=Float64)

Read the file `file` and create a `XiMultipoles` struct having
- as first element the comoving distances stored in the first column (with the input type `xdt`)
- as second element a vector that contains all the following columns 
  (with the input type `ydt`), which are the multipoles `L=0,1,2,...` of the TPCF considered

If the file start with comments (lines starting with #), set `comments = true`.
"""
struct XiMultipoles
     comdist::Vector{Float64}
     multipoles::Vector{Vector{Float64}}

     function XiMultipoles(file::String; comments::Bool=true, xdt::DataType=Float64, ydt::DataType=Float64)

          table = readdlm(file, comments=comments)
          ss = vecstring_to_vecnumbers(table[:, 1]; dt=xdt)
          all = [vecstring_to_vecnumbers(col; dt=ydt)
                 for col in eachcol(table[:, 2:end])]

          for col in all
               @assert length(ss) == length(col) "the columns must have all the same length"
          end

          new(ss, all)
     end
end





##########################################################################################92





function PS_multipole_GenWin(
     ximult::Union{XiMultipoles,String}, genwin::Union{GenericWindow,String};
     alg::Symbol=:fftlog, L::Int=0,
     cut_first_n::Int=0, cut_last_n::Int=0,
     kwargs...)

     @assert cut_first_n ≥ 0 "cut_first_n ≥ 0 must hold!"
     @assert cut_last_n ≥ 0 "cut_last_n ≥ 0 must hold!"
     @assert cut_first_n + cut_last_n < length(ss) "cut_first_n + cut_last_n < length(ss) must hold!"

     XIMULT = typeof(ximult) == String ? XiMultipoles(ximult) : ximult
     GENWIN = typeof(genwin) == String ? GenericWindow(genwin) : genwin

     a, b = 1 + cut_first_n, length(ss) - cut_last_n
     SS, VEC_FS = XIMULT.comdist[a:b], [vec[a:b] for vec in XIMULT.multipoles]

     if alg == :twofast

          ks, _ = GaPSE.TwoFAST_PS_multipole(SS, VEC_FS[1]; cut_first_n=0, cut_last_n=0, kwargs...)
          res = zeros(length(ks))

          for l in 0:length(XIMULT.multipoles)-1, l1 in 0:length(GENWIN.splines)-1

               w3j = convert(Float64, wigner3j(L, l, l1, 0, 0, 0))

               if (w3j ≈ 0.0) == false
                    qls = [GENWIN.splines[l1+1](s) for s in SS]
                    _, pks = GaPSE.TwoFAST_PS_multipole(SS, VEC_FS[l+1] .* qls; L=L, cut_first_n=0, cut_last_n=0, kwargs...)
                    res += w3j .* pks
               end
          end

          return ks, res

     elseif alg == :fftlog

          ks, _ = GaPSE.FFTLog_PS_multipole(SS, VEC_FS[1]; cut_first_n=0, cut_last_n=0, kwargs...)
          res = zeros(length(ks))

          for l in 0:length(XIMULT.multipoles)-1, l1 in 0:length(GENWIN.splines)-1

               w3j = convert(Float64, wigner3j(L, l, l1, 0, 0, 0))

               if (w3j ≈ 0.0) == false
                    qls = [GENWIN.splines[l1+1](s) for s in SS]
                    _, pks = GaPSE.FFTLog_PS_multipole(SS, VEC_FS[l+1] .* qls; L=L, cut_first_n=0, cut_last_n=0, kwargs...)
                    res += w3j .* pks
               end
          end

          return ks, res

     else
          throw(AssertionError(
               "The algorithm ':$alg' does not exist! The available ones are:\n" *
               "\t ':fftlog' (default), ':twofast' ."
          ))
     end
end
