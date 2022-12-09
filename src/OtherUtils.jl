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


function warning(io::IO, msg::String)
    red = "\033[1m\033[31m"
    printstyled(io, "WARNING: " * msg * "\n"; color=:red, bold=true)
    return
end
warning(msg::String) = warning(stdout, msg)



"""
    check_compatible_dicts(ref::Dict, b::Dict, name::String = "NO-NAME")

Compare the field of two `Dict` and check if the second one (`b`) is "compatible"
with the first one (`ref`), i.e.:
- checks if each of the key in `b` is also a key of `ref`;
- for each `key` of `b`:
    - if `typeof(ref[key]) <: Real` and `!(typeof(ref[k]) <: Union{Bool, Int})`,
      checks that `typeof(b[k]) <: Real && typeof(b[k]) ≠ Bool`
    - otherwise, checks that ` typeof(b[k]) == typeof(ref[k])`

If someone of the check mentioned is `false`, raise an `AssertionError`, otherwise
return `nothing`. 
The string `name` is only used inside the AssertionError messages for the correct
name of the input `b` dictionary.
"""
function check_compatible_dicts(ref::Dict, b::Dict, name::String="NO-NAME")
    for k in keys(b)
        @assert haskey(ref, k) ((typeof(k) == Symbol ? ":$k" : "$k") *
                                " is not a valid key for the Dict $(name)!\n" *
                                "the allowed ones are the following:\n $(keys(ref))")

        if (typeof(ref[k]) <: Real) && !(typeof(ref[k]) <: Union{Bool,Int})
            @assert typeof(b[k]) <: Real && typeof(b[k]) ≠ Bool
            "the value associated with the (valid) key $k the Dict $(name) " *
            "must be of type T<:Real (not Bool) !" *
            "\n The type of the input value $(b[k]) is instead " *
            "$(typeof(b[k]))"
        else
            @assert typeof(b[k]) == typeof(ref[k]) "the value associated" *
                                                   " with the (valid) key $k for the Dict $(name) " *
                                                   "must be of type $(typeof(ref[k]))!" *
                                                   "\n The type of the input value $(b[k]) is instead " *
                                                   "$(typeof(b[k]))"
        end
    end

    return nothing
end



##########################################################################################92



function my_println_vec(io::IO, vec::Vector{T}, name::String; N::Int=5) where {T}
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

function my_println_vec(vec::Vector{T}, name::String; N::Int=5) where {T}
    my_println_vec(stdout, vec, name; N=N)
end



##########################################################################################92



function my_println_dict(io::IO, dict::Dict; pref::String="", N::Int=3)
    @assert N > 1 "N must be an integer >1, not $N !"

    print(io, pref)
    for (i, key) in enumerate(keys(dict))
        print(io, "$key = $(dict[key]) \t ")
        (i % N ≠ 0) || (length(keys(dict)) ≠ i && print(io, "\n" * pref))
    end
    (length(keys(dict)) % N ≠ 0) || print(io, "\n")

    return nothing
end

function my_println_dict(dict::Dict; pref::String="", N::Int=3)
    my_println_dict(stdout, dict; pref=pref, N=N)
end;



##########################################################################################92



function parent_directory(s::String)
    splitted = split(s, "/")

    if length(splitted) == 1 || (splitted[end] == "" && length(splitted) == 2)
        return "./"
    end

    ss = string(split(s, "/")[begin:end-1] .* "/"...)

    if ss == s
        return string(split(ss, "/")[begin:end-2] .* "/"...)
    else
        return ss
    end
end

function check_parent_directory(s::String)
    pd = parent_directory(s)

    pd == "./" && (return nothing)

    @assert isdir(pd) "$pd is not an existing directory!"
end


function return_namefile(s::String)
    splitted = split(s, "/")
    @assert splitted[end] ≠ "" "$s is a valid name for a directory, not for a file!"
    name = string(splitted[end])
    return name
end

function check_namefile(s::String)
    name = return_namefile(s)
    splitted_2 = split(name, ".")
    @assert splitted_2[end] ≠ name "$name has no extension (like .txt, .dat, ...)!"
    ve = ["txt", "dat"]
    @assert splitted_2[end] ∈ ve "$(splitted_2[end]) is not a valid extrension; they are: \n $(ve)"
end

function check_group(s::String; valid_groups::Vector{String}=VALID_GROUPS)
    error = """ "$s" is not a valid group name; they are the following:\n$(valid_groups)"""
    @assert s ∈ valid_groups error
end

function check_fileisingroup(input::String, group::String; 
        valid_groups::Vector{String}=VALID_GROUPS, comments::Bool=true)

    check_group(group; valid_groups=valid_groups)
    l = LENGTH_VALID_GROUPS[findfirst(x -> x == group, valid_groups)]
    table = readdlm(input, comments=comments)
    xs = convert(Vector{Float64}, table[:, 1])
    all_YS = [begin
        @assert length(col) == length(xs) "mismatch in length of the columns." 
        @assert col[end] ≠ "" "mismatch in length of the columns." 
        convert(Vector{Float64}, col)
        end
        for col in eachcol(table[:, 2:end])]

    if !isnothing(l)
        er = "the column numbers inside $input, which is $(length(all_YS) + 1), " *
                "does not match with the length of the " *
                "chosen group $group, which is $l."
        @assert (length(all_YS) + 1 == l) er
    end
    for (i, ys) in enumerate(all_YS)
        err = "the column number $i has length = $(length(ys)) instead of "
        "$(length(xs)) (which is the length of the first column, the x-axis one). "
        @assert length(ys) == length(xs) err
    end
end


##########################################################################################92

function number_to_string(x::Number)
    if iszero(imag(x))
        return "$(real(x))"
    elseif iszero(real(x))
        return "$(imag(x))im"
    else
        if imag(x) > 0
            return "$(real(x))+$(imag(x))im"
        else
            return "$(real(x))-$(abs(imag(x)))im"
        end
    end
end


function vecstring_to_vecnumbers(v; dt::DataType = Float64)
    try 
        return convert(Vector{dt}, v) 
    catch e 
        [parse(dt, el) for el in v]
    end
end


##########################################################################################92


function readxy(input::String; comments::Bool=true, 
        xdt::DataType = Float64, ydt::DataType = Float64)
    table = readdlm(input, comments=comments)
    xs = vecstring_to_vecnumbers(table[:, 1]; dt = xdt)
    ys = vecstring_to_vecnumbers(table[:, 2]; dt = ydt)    
    return (xs, ys)
end;

function readxall(input::String; comments::Bool=true, 
        xdt::DataType = Float64, ydt::DataType = Float64)
    table = readdlm(input, comments=comments)
    xs =  vecstring_to_vecnumbers(table[:, 1]; dt = xdt)
    all = [vecstring_to_vecnumbers(col; dt = ydt)
           for col in eachcol(table[:, 2:end])]
    return (xs, all)
end;

function readxyall(input::String; comments::Bool=true, 
        xdt::DataType = Float64, ydt::DataType = Float64, zdt::DataType = Float64)
    table = readdlm(input, comments=comments)
    xs =  vecstring_to_vecnumbers(table[:, 1]; dt = xdt)
    ys = vecstring_to_vecnumbers(table[:, 2]; dt = ydt)
    all = [vecstring_to_vecnumbers(col; dt = zdt)
           for col in eachcol(table[:, 3:end])]
    return (xs, ys, all)
end;



##########################################################################################92



"""
    sample_subdivision_begin(x_min, x_stop, x_end; 
        frac_begin::Float64 = 0.5, N::Int = 100, ass::Bool = true)

Return a vector of `N+2` points inside the interval `x_min ≤ x ≤ x_max` linearly distributed
with two different sampling:
- in `x_min ≤ x ≤ x_stop` there are around `frac_begin * N` points;
- in `x_stop ≤ x ≤ x_max` there are around `(1.0 - frac_begin) * N` points.

`frac_begin` is then the fraction of the `N` points that is inside the LEFT INTERVAL.
If `ass::Bool` is set to `false` the assert checks on the input data will not be performed. 
"""
function sample_subdivision_begin(x_min, x_stop, x_max; frac_begin::Float64 = 0.5, N::Int = 100, ass::Bool = true)
    if ass == true
        @assert 0.0 < frac_begin < 1.0 "frac_begin must be in 0.0 < frac_begin < 1.0, frac_begin = $frac_begin is not valid!" 
        @assert x_min < x_stop "x_min < x_stop must hold, x_min = $x_min and x_stop = $x_stop do not!"
        @assert x_stop < x_max "x_stop < x_max must hold, x_stop = $x_stop and x_max = $x_max do not!"
        @assert N > 3 "N > 3 must hold, N = $N do not!"
    end

    return unique(
        vcat(
            range(x_min, x_stop; length = Int( ceil( frac_begin * N)) + 1) ,
            range(x_stop, x_max; length = Int( ceil( (1.0 - frac_begin) * N)) + 1)
        )
    )
end


"""
    sample_subdivision_middle(x_min, x_start, x_stop, x_max; 
        frac_middle::Float64 = 0.5, rel_frac_begin::Union{Float64, Nothing} = nothing, 
        N::Int = 100, ass::Bool = true)

Return a vector of `N+3` points inside the interval `x_min ≤ x ≤ x_max` linearly distributed
with three different sampling, depending on the values of `frac_middle` and `rel_frac_begin`.

If `rel_frac_begin == nothing`, defining the relative size of the segment `x_min ≤ x ≤ x_start`
compared to the "masked one" `x_min ≤ x ≤ x_start || x_stop ≤ x ≤ x_max` as 
`rel_prop_begin = (x_start - x_min) / (x_max - x_min - x_stop + x_start)`:
- in `x_min ≤ x ≤ x_start` there are around `(1.0 - frac_middle) * rel_prop_begin * N` points;
- in `x_start ≤ x ≤ x_stop` there are around `frac_middle * N` points;
- in `x_stop ≤ x ≤ x_max` there are around `(1.0 - frac_middle) * (1.0 - rel_prop_begin) * N` points.

If `rel_frac_begin` is instead a float inside the interval `0.0 < rel_frac_begin < 1.0`: 
- in `x_min ≤ x ≤ x_start` there are around `(1.0 - frac_middle) * rel_frac_begin * N` points;
- in `x_start ≤ x ≤ x_stop` there are around `frac_middle * N` points;
- in `x_stop ≤ x ≤ x_max` there are around `(1.0 - frac_middle) * (1.0 - rel_frac_begin) * N` points.

`frac_middle` is then the fraction of the `N` points that is inside the MIDDLE INTERVAL, while 
`rel_frac_begin` is the one inside the LEFT INTERVAL COMPARED TO THE MASKED TOTAL ONE.
If `ass::Bool` is set to `false` the assert checks on the input data will not be performed. 
"""
function sample_subdivision_middle(x_min, x_start, x_stop, x_max; 
        frac_middle::Float64 = 0.5, rel_frac_begin::Union{Float64, Nothing} = nothing, 
        N::Int = 100, ass::Bool = true)
    
    if ass == true
        @assert 0.0 < frac_middle < 1.0 "frac_middle must be in 0.0 < frac_middle < 1.0, frac_middle = $frac_middle is not valid!" 
        @assert isnothing(rel_frac_begin) || 0.0 < rel_frac_begin < 1.0 "rel_frac_begin must be in "*
            "0.0 < rel_frac_begin < 1.0, rel_frac_begin = $rel_frac_begin not valid!" 
        @assert x_min < x_start "x_min < x_start must hold, x_min = $x_min and x_start = $x_start do not!"
        @assert x_start < x_stop "x_start < x_stop must hold, x_start = $x_start and x_stop = $x_stop do not!"
        @assert x_stop < x_max "x_stop < x_max must hold, x_stop = $x_stop and x_max = $x_max do not!"
        @assert N > 5 "N > 5 must hold, N = $N do not!"
    end
    
    
    if isnothing(rel_frac_begin)
        rel_prop_begin = (x_start - x_min) / (x_max - x_min - x_stop + x_start)
        
        return unique(
            vcat(
                range(x_min, x_start; length = Int( ceil( (1.0 - frac_middle) * rel_prop_begin * N) ) + 1 ),
                range(x_start, x_stop; length = Int( ceil( frac_middle * N ) + 1 ) ),
                range(x_stop, x_max; length = Int( ceil( (1.0-frac_middle) * (1.0 - rel_prop_begin) * N) ) + 1 ),
            )
        )
    else
        return unique(
            vcat(
                range(x_min, x_start; length = Int( ceil( (1.0 - frac_middle) * rel_frac_begin * N) ) + 1 ),
                range(x_start, x_stop; length = Int( ceil( frac_middle * N ) + 1 ) ),
                range(x_stop, x_max; length = Int( ceil( (1.0 - frac_middle) * (1.0 - rel_frac_begin) * N) ) + 1 ),
            )
        )
    end
end