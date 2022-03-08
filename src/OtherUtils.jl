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



function check_compatible_dicts(ref::Dict, b::Dict, name::String)
    for k in keys(b)
    @assert haskey(ref, k) ((typeof(k) == Symbol ? ":$k" : "$k" ) * 
                " is not a valid key for the Dict $(name)!\n"*
            "the allowed ones are the following:\n $(keys(ref))")
        if (typeof(ref[k]) <: Real) && (typeof(ref[k]) ≠ Bool)
            
        @assert typeof(b[k]) <: Real "the value associated"*
        " with the (valid) key $k the Dict $(name) "*
            "must be of type T<:Real !"*
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
