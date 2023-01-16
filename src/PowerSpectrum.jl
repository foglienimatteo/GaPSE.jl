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

function PS_multipole(ss, fs;
     alg::Symbol=:fftlog,
     cut_first_n::Int=0, cut_last_n::Int=0,
     kwargs...)

     @assert length(ss) == length(fs) "length(ss) == length(fs) must hold!"
     @assert cut_first_n ≥ 0 "cut_first_n ≥ 0 must hold!"
     @assert cut_last_n ≥ 0 "cut_last_n ≥ 0 must hold!"
     @assert cut_first_n + cut_last_n < length(ss) "cut_first_n + cut_last_n < length(ss) must hold!"

     a, b = 1 + cut_first_n, length(ss) - cut_last_n

     if alg == :twofast
          return TwoFAST_PS_multipole(ss[a:b], fs[a:b]; cut_first_n=0, cut_last_n=0, kwargs...)

     elseif alg == :fftlog
          return FFTLog_PS_multipole(ss[a:b], fs[a:b]; cut_first_n=0, cut_last_n=0, kwargs...)

     else
          throw(AssertionError(
               "The algorithm ':$alg' does not exist! The available ones are:\n" *
               "\t ':fftlog' (default), ':twofast' ."
          ))
     end
end


function PS_multipole(input::String; comments=true, kwargs...)
     ss, fs = GaPSE.readxy(input; comments=comments)

     return PS_multipole(ss, fs; kwargs...)
end




"""
     PS_multipole(ss, fs; 
          pr::Bool = true, L::Int = 0, 
          alg::Symbol = :fftlog, 
          cut_first_n::Int = 0, cut_last_n::Int = 0, 
          kwargs...
          ) ::Tuple{Vector{Float64}, Vector{Float64}}

     PS_multipole(input::String; 
          kwargs...)

Return the `L`-order PS multipole through the
following Fast Fourier Transform and the effective redshift approximation:

```math
P_L(k) = \\frac{2 L + 1}{A^{'}} (-i)^L \\, \\phi(s_\\mathrm{eff}) \\int_0^\\infty 
        \\mathrm{d} s \\; s^2 \\, j_L(ks) \\, f_\\mathrm{in}(s) \\; ,
        \\quad \\; A^{'} = \\frac{1}{4\\,\\pi}
```


The second method reads the input file, takes the first column as `ss` and the second as `fs`
and recalls the first method.
     
Currenlty, there are two algorithms you can choose in order to perform the computation; you can choose 
which one to use through the keyword value `alg`:
- `alg = :fftlog` (default and recommended option) will employ the [FFTLog](https://github.com/marcobonici/FFTLog.jl) 
  algorithm.
- `alg = :twofast` will employ the TwoFAST `xicalc` function of the [TwoFAST](https://github.com/hsgg/TwoFAST.jl) 
  Julia package. Note that in the computation the integration range ``0\\leq s \\leq \\infty`` 
  is reduced to `int_s_min ≤ s ≤ int_s_max`. This alogrithm is not the ideal choise, because TwoFAST is conceived
  for the direction PS -> TPCF, while is not 100% trustworthy for the other way round.

IMPORTANT: no matter which algorithm you choose, you will need to give the input data in a
LOGARITHMICALLY DISTRIBUTED scale. A linear distribution does not fit for the algorithms to apply.

##  Optional arguments

Depending on the algorithm you choose, the options would change. The options in common are:
- `pr::Bool=true` : want to print the automatic messages to the screen?
- `L::Int=0` : which multipole order should I use for this computation? IT MUST MATCH 
  THE MULTIPOLE ORDER OF THE INPUT TPCF!
- `cut_first_n::Int=0` and `cut_last_n::Int=0` : you can cout the first and/or last n elements
  of the input data, if they are highly irregular.

The specific ones for `alg = :fftlog` are:
- `ν::Union{Float64,Nothing} = nothing` : bias parameter, i.e. exponent used to "balance" the curve;
  if `nothing`, will be set automatically to `1.5`
- `n_extrap_low::Int = 500` and `n_extrap_high::Int = 500` : number of points to concatenate on the left/right
  of the input x-axis `ss` vector, logarithmically distributed with the same ratio of the left/right-edge
  elements of `ss`.
- `n_pad::Int = 500` : number of zeros to be concatenated both on the left and
  on the right of the input function. They stabilize a lot the algorithm.

The specific ones for `alg = :twofast` are:
- `epl::Bool=true` : do you want to extend the edges of the input vectors using two fitted
  power-laws (obtained from `EPLs`)
- `N_left::Int = 12` and `N_right::Int = 12` : number of points from left
  right edges to be used for the power law fitting in `EPLs`. They matters only
  if in the given input file ξ is not defined until the extremes of integration
  `int_s_min` and `int_s_max`.
- `int_s_min::Float64 = 1e-1` and `int_s_max::Float64 = 1e3`: extremes of integration; if `epl`
  is set to `false`, their values will be automatically set to `min(ss...)` and `max(ss...)`
  respectively. Their values do matter only if `epl=true`. 
- `p0_left=[-2.0, 1.0]` and `p0_right=[-2.0, 1.0]`: vectors with the initial values for the left/right 
  power-law fitting of `EPLs`; the power-law is in the form ``y = f(x) = b * x^s``, so the first vector 
  value is the initial value of ``s`` (and of course the second is the one of ``b``).
- `k0::Union{Nothing,Float64} = nothing` : starting point for the `xicalc` function; if `nothing`, 
  it will be set `k0 = 1.0 / max(ss...)`
- `right::Union{Float64,Nothing} = nothing` : do you want to cut the output elements with 
  `ks .> right`? if set to `nothing`, no cut will be done.
- `N::Int = 1024` : number of points to be used in Fourier transform 


## Analytical derivation 

The analytical expression previously showed can be easily obtained from the 
standard one:
```math
\\begin{split}
    P_L(k) = &\\frac{2 L + 1}{A} (-i)^L \\, 
        \\int_0^\\infty \\mathrm{d} s_1 \\; s_1^2 
        \\int_0^\\infty \\mathrm{d} s \\; s^2 
        \\int_{-1}^{+1} \\mathrm{d} \\mu \\;
        j_L(ks) \\, \\xi(s_1, s, \\mu) \\, \\phi(s_1) \\, \\phi(s_2) \\,
        \\mathcal{L}_L(\\mu) F\\left(\\frac{s}{s_1}, \\mu \\right) \\\\
        &\\mathrm{with} \\; \\;s_2 = s_2(s_1, s, μ) = \\sqrt{s_1^2 + s^2 + 2s_1s\\mu}
        \\; 
        , \\quad A(s_\\mathrm{max}, s_\\mathrm{min}, \\theta_\\mathrm{max}) 
        \\frac{
          V(s_\\mathrm{max}, s_\\mathrm{min}, \\theta_\\mathrm{max})
          }{4 \\, \\pi^2}
\\end{split}
```

with the definition

```math
f_\\mathrm{in}(s_1, s) =  \\int_{-1}^{+1} \\mathrm{d} \\mu \\;
        \\xi(s_1, s, \\mu) \\, \\phi(s_2) \\,
        \\mathcal{L}_L(\\mu) \\, F\\left(\\frac{s}{s_1}, \\mu \\right)
```
and the application of the effective redshift approximation.

## Returns

A `Tuple{Vector{Float64}, Vector{Float64}}` with:
- the `k` values vector as first element;
- the correspoding PS `pk` values vector as second one.


See also: [`V_survey`](@ref), [`A`](@ref), [`A_prime`](@ref),
[`EPLs`](@ref),  [`print_PS_multipole`](@ref)
"""
PS_multipole



##########################################################################################92



function print_PS_multipole(ss, fs, out::String;
     L::Int=0, pr::Bool=true, alg::Symbol=:fftlog, kwargs...)

     check_parent_directory(out)
     check_namefile(out)

     pr && println("""\nI'm computing the PS_multipole from the two input vectors.""")

     time_1 = time()
     vec = PS_multipole(ss, fs; L=L, pr=pr, alg=alg, kwargs...)
     time_2 = time()

     N = length(vec[1])

     isfile(out) && run(`rm $out`)
     open(out, "w") do io
          println(io, BRAND)

          println(io, "# Power Spectrum Multipole computation from two input vectors.")
          println(io, "#\n# For this PS_multipole computation we set: ")
          println(io, "# \t multipole degree in consideration L = $L")
          println(io, "# \t algorithm chosen for the computation: :$alg")
          println(io, "# \t #points used in Fourier transform N = $N")
          println(io, "# computational time needed (in s) : $(@sprintf("%.4f", time_2-time_1))")
          print(io, "# kwards passed to \"print_PS_multipole\": ")

          if isempty(kwargs)
               println(io, "none")
          else
               print(io, "\n")
               for key in keys(kwargs)
                    println(io, "# \t\t$(key) = $(kwargs[key])")
               end
          end
          println(io, "# ")
          println(io, "# k [h_0/Mpc] \t \t  P [(Mpc/h_0)^3]")
          for (k, pk) in zip(vec[1], vec[2])
               println(io, "$k \t " * GaPSE.number_to_string(pk))
          end
     end
end


function print_PS_multipole(input::String, out::String;
     comments=true, kwargs...)
     ss, fs = GaPSE.readxy(input; comments=comments)

     print_PS_multipole(ss, fs, out; kwargs...)
end



"""
     print_PS_multipole(ss, fs, out::String;
          L::Int=0, pr::Bool=true, alg::Symbol=:fftlog, kwargs...)
     print_PS_multipole(input::String, out::String;
          kwargs...)

Takes in input a filename `input` where is stored a TPCF multipole,
calculate the `L`-order PS multipole through the
following Fast Fourier Transform and the effective redshift approximation

```math
P_L(k) = \\frac{2 L + 1}{A^{'}} (-i)^L \\, \\phi(s_\\mathrm{eff}) \\int_0^\\infty 
        \\mathrm{d} s \\; s^2 \\, j_L(ks) \\, f_\\mathrm{in}(s) \\; ,
        \\quad \\; A^{'} = \\frac{1}{4\\,\\pi}
```

where ``f_\\mathrm{in}`` is the function samples by `ss` and `xis`,
and save it in the file `out`, together with the options used for the computation.

The second method reads the input file, takes the first column as `ss` and the second as `fs`
and recalls the first method.


## Optional arguments

Depending on the algorithm you choose, the options would change. The options in common are:
- `pr::Bool=true` : want to print the automatic messages to the screen?
- `L::Int=0` : which multipole order should I use for this computation? IT MUST MATCH 
  THE MULTIPOLE ORDER OF THE INPUT TPCF!
- `cut_first_n::Int=0` and `cut_last_n::Int=0` : you can cout the first and/or last n elements
  of the input data, if they are highly irregular.
- `alg::Symbol = :fftlog` : algorithm to be used for the computation. Currenlty, there are two algorithms 
  you can coose in order to perform the computation:
  - `alg = :fftlog` (default and recommended option) will employ the [FFTLog](https://github.com/marcobonici/FFTLog.jl) 
    algorithm.
  - `alg = :twofast` will employ the TwoFAST `xicalc` function of the [TwoFAST](https://github.com/hsgg/TwoFAST.jl) 
    Julia package. Note that in the computation the integration range ``0\\leq s \\leq \\infty`` 
    is reduced to `int_s_min ≤ s ≤ int_s_max`. This alogrithm is not the ideal choise, because TwoFAST is conceived
    for the direction PS -> TPCF, while is not 100% trustworthy for the other way round.


The specific ones for `alg = :fftlog` are:
- `ν::Union{Float64,Nothing} = nothing` : bias parameter, i.e. exponent used to "balance" the curve;
  if `nothing`, will be set automatically to `1.5`
- `n_extrap_low::Int = 500` and `n_extrap_high::Int = 500` : number of points to concatenate on the left/right
  of the input x-axis `ss` vector, logarithmically distributed with the same ratio of the left/right-edge
  elements of `ss`.
- `n_pad::Int = 500` : number of zeros to be concatenated both on the left and
  on the right of the input function. They stabilize a lot the algorithm.

The specific ones for `alg = :twofast` are:
- `epl::Bool=true` : do you want to extend the edges of the input vectors using two fitted
  power-laws (obtained from `EPLs`)
- `N_left::Int = 12` and `N_right::Int = 12` : number of points from left
  right edges to be used for the power law fitting in `EPLs`. They matters only
  if in the given input file ξ is not defined until the extremes of integration
  `int_s_min` and `int_s_max`.
- `int_s_min::Float64 = 1e-1` and `int_s_max::Float64 = 1e3`: extremes of integration; if `epl`
  is set to `false`, their values will be automatically set to `min(ss...)` and `max(ss...)`
  respectively. Their values do matter only if `epl=true`. 
- `p0_left=[-2.0, 1.0]` and `p0_right=[-2.0, 1.0]`: vectors with the initial values for the left/right 
  power-law fitting of `EPLs`; the power-law is in the form ``y = f(x) = b * x^s``, so the first vector 
  value is the initial value of ``s`` (and of course the second is the one of ``b``).
- `k0::Union{Nothing,Float64} = nothing` : starting point for the `xicalc` function; if `nothing`, 
  it will be set `k0 = 1.0 / max(ss...)`
- `right::Union{Float64,Nothing} = nothing` : do you want to cut the output elements with 
  `ks .> right`? if set to `nothing`, no cut will be done.
- `N::Int = 1024` : number of points to be used in Fourier transform 

See also: [`V_survey`](@ref), [`A`](@ref), [`A_prime`](@ref),
[`EPLs`](@ref), [`PS_multipole`](@ref)
"""
print_PS_multipole



##########################################################################################92



"""
     function all_PS_multipole(input::String,
          group::String=VALID_GROUPS[end];
          L::Int = 0, pr::Bool = true, 
          alg::Symbol=:fftlog, kwargs...
          ) ::Tuple{Vector{Float64}, Vector{Vector{Float64}}}

Given an input file where the first column is the x-axis data one and all the 
following columns are the corresponding y-data ones, this function computes all the 
Power Spectra of each y-data column and return a Tuple containing
- as first element, the `ks` values, common to all the PS
- as second element, a vector where in each position there is the Power Spectra corresponding 
  to the associated inputy y-data. 

The `group::String=VALID_GROUPS[end]` argument allow you to specify the group of the input TPCF, 
if they were computed through GAPSE. The allowed values for this argument are:
`$(string(GaPSE.VALID_GROUPS .* " , "...))`

If you choose a group pay attention that the number of input TPCF must match the group number 
(16, 25, 20 and 20 respectively). The last group name (which is also the default value) is used in 
case the input xis do not belog to a specific group (and so no predefined number is expected).


## Optional arguments

Depending on the algorithm you choose, the options would change. The options in common are:
- `pr::Bool=true` : want to print the automatic messages to the screen?
- `L::Int=0` : which multipole order should I use for this computation? IT MUST MATCH 
  THE MULTIPOLE ORDER OF THE INPUT TPCF!
- `cut_first_n::Int=0` and `cut_last_n::Int=0` : you can cout the first and/or last n elements
  of the input data, if they are highly irregular.
- `alg::Symbol = :fftlog` : algorithm to be used for the computation. Currenlty, there are two algorithms 
  you can coose in order to perform the computation:
  - `alg = :fftlog` (default and recommended option) will employ the [FFTLog](https://github.com/marcobonici/FFTLog.jl) 
    algorithm.
  - `alg = :twofast` will employ the TwoFAST `xicalc` function of the [TwoFAST](https://github.com/hsgg/TwoFAST.jl) 
    Julia package. Note that in the computation the integration range ``0\\leq s \\leq \\infty`` 
    is reduced to `int_s_min ≤ s ≤ int_s_max`. This alogrithm is not the ideal choise, because TwoFAST is conceived
    for the direction PS -> TPCF, while is not 100% trustworthy for the other way round.

The specific ones for `alg = :fftlog` are:
- `ν::Union{Float64,Nothing} = nothing` : bias parameter, i.e. exponent used to "balance" the curve;
  if `nothing`, will be set automatically to `1.5`
- `n_extrap_low::Int = 500` and `n_extrap_high::Int = 500` : number of points to concatenate on the left/right
  of the input x-axis `ss` vector, logarithmically distributed with the same ratio of the left/right-edge
  elements of `ss`.
- `n_pad::Int = 500` : number of zeros to be concatenated both on the left and
  on the right of the input function. They stabilize a lot the algorithm.

The specific ones for `alg = :twofast` are:
- `epl::Bool=true` : do you want to extend the edges of the input vectors using two fitted
  power-laws (obtained from `EPLs`)
- `N_left::Int = 12` and `N_right::Int = 12` : number of points from left
  right edges to be used for the power law fitting in `EPLs`. They matters only
  if in the given input file ξ is not defined until the extremes of integration
  `int_s_min` and `int_s_max`.
- `int_s_min::Float64 = 1e-1` and `int_s_max::Float64 = 1e3`: extremes of integration; if `epl`
  is set to `false`, their values will be automatically set to `min(ss...)` and `max(ss...)`
  respectively. Their values do matter only if `epl=true`. 
- `p0_left=[-2.0, 1.0]` and `p0_right=[-2.0, 1.0]`: vectors with the initial values for the left/right 
  power-law fitting of `EPLs`; the power-law is in the form ``y = f(x) = b * x^s``, so the first vector 
  value is the initial value of ``s`` (and of course the second is the one of ``b``).
- `k0::Union{Nothing,Float64} = nothing` : starting point for the `xicalc` function; if `nothing`, 
  it will be set `k0 = 1.0 / max(ss...)`
- `right::Union{Float64,Nothing} = nothing` : do you want to cut the output elements with 
  `ks .> right`? if set to `nothing`, no cut will be done.
- `N::Int = 1024` : number of points to be used in Fourier transform 



See also: [`EPLs`](@ref), [`PS_multipole`](@ref)
"""
function all_PS_multipole(input::String,
     group::String=VALID_GROUPS[end];
     L::Int=0, pr::Bool=true,
     alg::Symbol=:fftlog, kwargs...)

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

     time_1 = time()

     ks, VEC = if alg == :twofast
          TwoFAST_all_PS_multipole(input,
               group; L=L, pr=false, kwargs...)

     elseif alg == :fftlog
          FFTLog_all_PS_multipole(input,
               group; L=L, pr=false, kwargs...)

     else
          throw(AssertionError(
               "The algorithm ':$alg' does not exist! The available ones are:\n" *
               "\t ':fftlog' (default), ':twofast' ."
          ))
     end


     time_2 = time()

     pr && println("\ntime needed for all the Power Spectra L = $L computation [in s] = $(time_2-time_1)\n")

     return ks, VEC
end



"""
     print_all_PS_multipole(input::String, out::String,
          group::String = VALID_GROUPS[end]; 
          L::Int = 0, pr::Bool = true, 
          alg::Symbol = :fftlog,
          kwargs...)

Given an `input`` file where the first column is the x-axis data one and all the 
following columns are the corresponding y-data ones, this function computes all the 
Power Spectra of each y-data column and print in a file named `out` 
- as first column, the `ks` values, common to all the PS
- as folowing columns, a vector where in each position there is the Power Spectra corresponding 
  to the associated inputy y-data. 

The `group::String=VALID_GROUPS[end]` argument allow you to specify the group of the input TPCF, 
if they were computed through GAPSE. The allowed values for this argument are:
`$(string(GaPSE.VALID_GROUPS .* " , "...))`

If you choose a group pay attention that the number of input TPCF must match the group number 
(16, 25, 20 and 20 respectively). The last group name (which is also the default value) is used in 
case the input xis do not belog to a specific group (and so no predefined number is expected).

## Optional arguments


Depending on the algorithm you choose, the options would change. The options in common are:
- `pr::Bool=true` : want to print the automatic messages to the screen?
- `L::Int=0` : which multipole order should I use for this computation? IT MUST MATCH 
  THE MULTIPOLE ORDER OF THE INPUT TPCF!
- `cut_first_n::Int=0` and `cut_last_n::Int=0` : you can cout the first and/or last n elements
  of the input data, if they are highly irregular.
- `alg::Symbol = :fftlog` : algorithm to be used for the computation. Currenlty, there are two algorithms 
  you can coose in order to perform the computation:
  - `alg = :fftlog` (default and recommended option) will employ the [FFTLog](https://github.com/marcobonici/FFTLog.jl) 
    algorithm.
  - `alg = :twofast` will employ the TwoFAST `xicalc` function of the [TwoFAST](https://github.com/hsgg/TwoFAST.jl) 
    Julia package. Note that in the computation the integration range ``0\\leq s \\leq \\infty`` 
    is reduced to `int_s_min ≤ s ≤ int_s_max`. This alogrithm is not the ideal choise, because TwoFAST is conceived
    for the direction PS -> TPCF, while is not 100% trustworthy for the other way round.

The specific ones for `alg = :fftlog` are:
- `ν::Union{Float64,Nothing} = nothing` : bias parameter, i.e. exponent used to "balance" the curve;
  if `nothing`, will be set automatically to `1.5`
- `n_extrap_low::Int = 500` and `n_extrap_high::Int = 500` : number of points to concatenate on the left/right
  of the input x-axis `ss` vector, logarithmically distributed with the same ratio of the left/right-edge
  elements of `ss`.
- `n_pad::Int = 500` : number of zeros to be concatenated both on the left and
  on the right of the input function. They stabilize a lot the algorithm.

The specific ones for `alg = :twofast` are:
- `epl::Bool=true` : do you want to extend the edges of the input vectors using two fitted
  power-laws (obtained from `EPLs`)
- `N_left::Int = 12` and `N_right::Int = 12` : number of points from left
  right edges to be used for the power law fitting in `EPLs`. They matters only
  if in the given input file ξ is not defined until the extremes of integration
  `int_s_min` and `int_s_max`.
- `int_s_min::Float64 = 1e-1` and `int_s_max::Float64 = 1e3`: extremes of integration; if `epl`
  is set to `false`, their values will be automatically set to `min(ss...)` and `max(ss...)`
  respectively. Their values do matter only if `epl=true`. 
- `p0_left=[-2.0, 1.0]` and `p0_right=[-2.0, 1.0]`: vectors with the initial values for the left/right 
  power-law fitting of `EPLs`; the power-law is in the form ``y = f(x) = b * x^s``, so the first vector 
  value is the initial value of ``s`` (and of course the second is the one of ``b``).
- `k0::Union{Nothing,Float64} = nothing` : starting point for the `xicalc` function; if `nothing`, 
  it will be set `k0 = 1.0 / max(ss...)`
- `right::Union{Float64,Nothing} = nothing` : do you want to cut the output elements with 
  `ks .> right`? if set to `nothing`, no cut will be done.
- `N::Int = 1024` : number of points to be used in Fourier transform 


See also: [`EPLs`](@ref), [`PS_multipole`](@ref)
"""
function print_all_PS_multipole(input::String, out::String,
     group::String=VALID_GROUPS[end];
     L::Int=0, pr::Bool=true,
     alg::Symbol=:fftlog,
     kwargs...)

     check_parent_directory(out)
     check_namefile(out)
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

     time_1 = time()

     ks, VEC = all_PS_multipole(input, group;
          alg=alg, L=L, pr=false, kwargs...)

     time_2 = time()

     N = length(ks)

     pr && println("\ntime needed for all the Power Spectra computation [in s] = $(time_2-time_1)\n")

     isfile(out) && run(`rm $out`)
     open(out, "w") do io
          println(io, BRAND)

          print(io, "#\n# Power Spectra Multipole computation ")
          if group == "GNC"
               println(io, "for the Galaxy Number Counts GR effect" *
                           "\n# from the file: $input")
          elseif group == "LD"
               println(io, "for the Luminosity Distance perturbations GR effect" *
                           "\n# from the file: $input")
          elseif group == "GNCxLD"
               println(io, "for the cross correlations between \n#" *
                           "Galaxy Number Counts and Luminosity Distance perturbations " *
                           "from the file:\n# $input")
          elseif group == "LDxGNC"
               println(io, "for the cross correlations between \n#" *
                           "Luminosity Distance perturbations and Galaxy Number Counts " *
                           "from the file:\n# $input")
          else
               println(io, "without any specific group considered" *
                           "\n# from the file: $input")
          end

          println(io, "#\n# For this PS_multipole computation we set: ")
          println(io, "# \t multipole degree in consideration L = $L")
          println(io, "# \t algorithm chosen for the computation: :$alg")
          println(io, "# \t #points used in Fourier transform N = $N")
          println(io, "# overall computational time needed (in s) : $(@sprintf("%.4f", time_2-time_1))")
          print(io, "# kwards passed to \"print_all_PS_multipole\": ")

          if isempty(kwargs)
               println(io, "none")
          else
               print(io, "\n")
               for key in keys(kwargs)
                    println(io, "# \t\t$(key) = $(kwargs[key])")
               end
          end
          println(io, "# ")
          println(io, "# (all the following Power Spectra are measured in (Mpc/h_0)^3)")

          if group ≠ VALID_GROUPS[end]
               effs = group == "GNC" ? GR_EFFECTS_GNC :
                      group == "LD" ? GR_EFFECTS_LD :
                      group == "GNCxLD" ? GR_EFFECTS_GNCxLD :
                      group == "LDxGNC" ? GR_EFFECTS_LDxGNC :
                      throw(ErrorException("how did you arrive here???"))

               println(io, "# 1: k [h_0/Mpc] \t 2: P_SUM \t " *
                           join([string(i) for i in 3:length(effs)+2] .*
                                ": P_" .* effs .* " \t "))
          else
               println(io, "# 1: k [h_0/Mpc] \t" *
                           join([string(i) for i in 2:length(VEC)+1] .*
                                ": P_" .*
                                [string(i) for i in 2:length(VEC)+1] .* " \t ")
               )
          end

          for (i, k) in enumerate(ks)
               println(io, "$k \t" * join([GaPSE.number_to_string(v[i]) * " \t " for v in VEC]))
          end
     end
end
