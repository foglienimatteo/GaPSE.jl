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


##########################################################################################92


function warning(io::IO, msg::String)
    red = "\033[1m\033[31m" 
    printstyled(io, "WARNING: " * msg * "\n"; color=:red, bold=true)
    return
end
warning(msg::String) = warning(stdout, msg)



function check_compatible_dicts(ref::Dict, b::Dict, name::String = "NO-NAME")
    for k in keys(b)
        @assert haskey(ref, k) ((typeof(k) == Symbol ? ":$k" : "$k" ) * 
                " is not a valid key for the Dict $(name)!\n"*
                "the allowed ones are the following:\n $(keys(ref))")

        if (typeof(ref[k]) <: Real) && !(typeof(ref[k]) <: Union{Bool, Integer})
            @assert typeof(b[k]) <: Real && typeof(b[k]) ≠ Bool
                "the value associated with the (valid) key $k the Dict $(name) "*
                    "must be of type T<:Real (not Bool) !"*
                "\n The type of the input value $(b[k]) is instead "* 
                "$(typeof(b[k]))"
        else
            @assert typeof(b[k]) == typeof(ref[k]) "the value associated"*
            " with the (valid) key $k for the Dict $(name) "*
                "must be of type $(typeof(ref[k]))!"*
            "\n The type of the input value $(b[k]) is instead "* 
            "$(typeof(b[k]))"
        end
    end

    return nothing
end


function my_println_vec(io::IO, vec::Vector{T}, name::String; N::Integer = 5) where {T}
     @assert N > 1 "N must be an integer >1, not $N !"

     println(io, name * " = [")
     for (i, el) in enumerate(vec)
          print(io, string(el) * " , ")
          (i % N ≠ 0) || print(io, "\n")
     end
     (length(vec) % N ≠ 0) && print(io, "\n")
     print(io, "];\n")

     return nothing
end

function my_println_vec(vec::Vector{T}, name::String; N::Integer = 5) where {T}
     my_println_vec(stdout, vec, name; N = N)
end
