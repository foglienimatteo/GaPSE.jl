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
# MERCHANTABILITY or FITNEss FOR A PARTICULAR PURPOSE. See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with GaPSE. If not, see <http://www.gnu.org/licenses/>.
#


function integrand_ξ_LD_multipole(s1, s, μ, effect::Function, cosmo::Cosmology;
    L::Int=0, use_windows::Bool=true, kwargs...)

    s2_value = s2(s1, s, μ)
    y_value = y(s1, s, μ)

    res = if use_windows == true
        #println("s1 = $s1 \\t s2 = $(s2(s1, s, μ)) \\t  y=$(y(s1, s, μ))")
        int = effect(s1, s2_value, y_value, cosmo; kwargs...)
        #print(" F_val = ", spline_integrF(s, μ, cosmo.windowFint))
        #isapprox(cost, 3.1749350588720085e10; rtol=1e-5) || print(" \t cost = " , cosmo.WFI_norm)
        #isapprox(Pl(μ, L), 1.0; rotl=1e-4) || println(" \t Lp = " , Pl(μ, L) )
        #println("IF = $val")

        int .* (spline_integrF(s, μ, cosmo.windowFint) / cosmo.WFI_norm * Pl(μ, L))
        #=
        ϕ_s2 = ϕ(s2_value, cosmo.s_min, cosmo.s_max)
        (ϕ_s2 > 0.0) || (return 0.0)
        int = effect(s1, s2_value, y_value, cosmo; kwargs...)
        int .* (ϕ_s2 * spline_F(s / s1, μ, cosmo.windowF) * Pl(μ, L))
        =#
    else
        int = effect(s1, s2_value, y_value, cosmo; kwargs...)
        int .* Pl(μ, L)
    end

    return (2.0 * L + 1.0) / 2.0 * res
end


function integrand_ξ_LD_multipole(s1, s, μ, effect::String, cosmo::Cosmology; kwargs...)
    error = "$effect is not a valid GR effect name.\n" *
            "Valid GR effect names are the following:\n" *
            string(GaPSE.GR_EFFECTS_LD .* " , "...)
    @assert (effect ∈ GaPSE.GR_EFFECTS_LD) error

    return integrand_ξ_LD_multipole(s1, s, μ, DICT_GR_ξs_LD[effect], cosmo; kwargs...)
end



"""
    integrand_ξ_LD_multipole(s1, s, μ, effect::Function, cosmo::Cosmology;
        L::Int = 0, use_windows::Bool = true, kwargs...)

    integrand_ξ_LD_multipole(s1, s, μ, effect::String, cosmo::Cosmology; kwargs...)

Return the integrand on ``\\mu = \\hat{\\mathbf{s}}_1 \\cdot \\hat{\\mathbf{s}}`` 
of the chosen perturbation of the Luminosity Distance (LD) Two-Point Correlation Function (TPCF)
term, i.e. the following function ``f(s_1, s, \\mu)``:

```math
     f_L(s_1, s, \\mu) = \\xi \\left(s_1, s_2, y\\right) 
          \\, \\mathcal{L}_L(\\mu) \\, \\times 
     \\begin{cases} 
          \\frac{1}{\\mathcal{N}}\\mathcal{F}(s, \\mu) \\quad \\mathrm{use\\_windows == true} \\\\
          1 \\quad\\quad \\mathrm{use\\_windows == false}
     \\end{cases}
```

where:
- ``s_2 = s_2(s_1, s, \\mu) = \\sqrt{s_1^2 + s^2 + 2 \\, s_1 \\, s \\, \\mu}`` 
- ``y = \\cos{\\theta} = \\hat{\\mathbf{s}}_1 \\cdot \\hat{\\mathbf{s}}_2 = \\frac{\\mu \\, s + s_1}{s_2(s_1, s, \\mu)}``
- ``\\xi`` is the corresponding LD TPCF effect
- ``\\mathcal{L}_L(\\mu)`` is the Legendre polynomial of order ``L``
- ``\\mathcal{F}(s, \\mu)`` is the integrated window function stored in `cosmo::Cosmology` (check the documentation of `WindowFIntegrated`)
- ``\\mathcal{N}`` is the integrated window function norm (check the documentation of `WindowFIntegrated`)

In the former method you have to pass as an input the LD `effect` function you want 
to integrate, while in the (recommended) latter one it's necessary to specify the
name of the LD TPCF term among the following:

`$(string(GR_EFFECTS_LD .* " , "...))`

to which correspond the following functions, respectively:

`$(string(string.(VEC_ξs_LD) .* " , "...))`

Internally, the second method recall the first, so all the keyword arguments are in common.

Note that these is an internal conversion of coordiate sistems
from `(s1, s, μ)` to `(s1, s2, y)` thorugh the functions `y` and `s2`. The inverse relations are:
- ``s = s(s_1, s_2, y) = \\sqrt{s_1^2 + s_2^2 - 2 \\, s_1 \\, s_2 \\, y}``
- ``\\mu = \\cos{\\theta} = \\hat{\\mathbf{s}}_1 \\cdot \\hat{\\mathbf{s}} = \\frac{y \\, s_2 - s_1}{s(s_1, s_2, y)}``.


## Inputs

- `s1`: the comoving distance where must be evaluated the integral

- `s`: the comoving distance from `s1` where must be evaluated the integral

- `μ`: the cosine between `s1` and `s` where must be evaluated the integral

- `effect::Function` or `effect::String` : name/function you want to consider, among the ones
  inside `GaPSE.GR_EFFECTS_LD`/`GaPSE.VEC_ξs_LD` respectively.

- `cosmo::Cosmology`: cosmology to be used in this computation


## Optional arguments 

- `L::Int = 0`: order of the Legendre polynomial to be used

- `use_windows::Bool = false`: tells if the integrand must consider the two
   window function ``\\phi`` and ``\\mathcal{F}``

- `kwargs...` : other keyword arguments that will be passed to the selected 
  GR TPCF effect (`ξ_LD_Doppler`, `ξ_LD_Lensing`, ...)


See also: [`ξ_LD_multipole`](@ref), [`map_ξ_LD_multipole`](@ref),
[`print_map_ξ_LD_multipole`](@ref),
[`WindowFIntegrated`](@ref), [`ϕ`](@ref), [`Cosmology`](@ref), 
[`y`](@ref), [`s2`](@ref), [`GR_EFFECTS_LD`](@ref), [`GaPSE.VEC_ξs_LD`](@ref)
"""
integrand_ξ_LD_multipole


##########################################################################################92




function ξ_LD_multipole(
    s1, s, effect::Function, cosmo::Cosmology;
    L::Int=0, alg::Symbol=:lobatto,
    use_windows::Bool=true,
    N_lob::Int=100, N_trap::Int=200,
    atol_quad::Float64=0.0, rtol_quad::Float64=1e-2,
    enhancer::Float64=1e6,
    kwargs...)

    error = "$(string(effect)) is not a valid GR effect function for galaxy number counts.\n" *
            "Valid GR effect functions for galaxy number counts are the following:\n" *
            string(string.(VEC_ξs_LD) .* " , "...)
    @assert (effect ∈ VEC_ξs_LD) error
    @assert alg ∈ VALID_INTEGRATION_ALGORITHM ":$alg is not a valid Symbol for \"alg\"; they are: \n\t" *
                                            "$(":".*string.(VALID_INTEGRATION_ALGORITHM) .* vcat([" , " for i in 1:length(VALID_INTEGRATION_ALGORITHM)-1], " .")... )"

    @assert N_trap > 2 "N_trap must be >2,  N_trap = $N_trap is not!"
    @assert N_lob > 2 "N_lob must be >2,  N_lob = $N_lob is not!"
    @assert atol_quad ≥ 0.0 "atol_quad must be ≥ 0.0,  atol_quad = $atol_quad is not!"
    @assert rtol_quad ≥ 0.0 "rtol_trap must be ≥ 0.0,  rtol_quad = $rtol_quad is not!"
    @assert L ≥ 0 "L must be ≥ 0, L = $L is not!"


    orig_f(μ) = enhancer * integrand_ξ_LD_multipole(s1, s, μ, effect, cosmo;
        L=L, use_windows=use_windows, kwargs...)

    int = if alg == :lobatto
        xs, ws = gausslobatto(N_lob)
        dot(ws, orig_f.(xs))

    elseif alg == :quad
        quadgk(μ -> orig_f(μ), -1.0, 1.0; atol=atol_quad, rtol=rtol_quad)[1]

    elseif alg == :trap
        μs = union(
            range(-1.0, -0.98, length=Int(ceil(N_trap / 3) + 1)),
            range(-0.98, 0.98, length=Int(ceil(N_trap / 3) + 1)),
            range(0.98, 1.0, length=Int(ceil(N_trap / 3) + 1))
        )
        #μs = range(-1.0 + 1e-6, 1.0 - 1e-6, length=N_trap)
        orig_fs = orig_f.(μs)
        trapz(μs, orig_fs)

    else
        throw(AssertionError("how the hell did you arrive here?"))
    end

    return int / enhancer
end



function ξ_LD_multipole(s1, s, effect::String, cosmo::Cosmology; kwargs...)
    error = "$effect is not a valid GR effect name for galaxy number counts.\n" *
            "Valid GR effect names for galaxy number counts are the following:\n" *
            string(GaPSE.GR_EFFECTS_LD .* " , "...)

    @assert (effect ∈ GaPSE.GR_EFFECTS_LD) error
    ξ_LD_multipole(s1, s, DICT_GR_ξs_LD[effect], cosmo; kwargs...)
end



"""
    ξ_LD_multipole(
        s1, s, effect::Function, cosmo::Cosmology;
        L::Int = 0, alg::Symbol = :lobatto, 
        use_windows::Bool = true, 
        N_lob::Int = 100, N_trap::Int = 200, 
        atol_quad::Float64 = 0.0, rtol_quad::Float64 = 1e-2,
        enhancer::Float64 = 1e6, 
        kwargs...) ::Float64

    ξ_LD_multipole(s1, s, effect::String, cosmo::Cosmology; 
        kwargs...) ::Float64

Evaluate the multipole of order `L` of the chosen perturbation of the Luminosity 
Distance (LD) Two-Point Correlation Function (TPCF) 
term i.e. the following function ``\\xi_L(s_1, s)``:

```math
    \\xi_L(s_1, s) = \\frac{2 L + 1}{2} \\int_{-1}^{+1} \\mathrm{d}\\mu \\; 
        \\xi (s_1, s_2, \\cos{\\theta}) \\, \\mathcal{L}_L(\\mu) \\,  \\times
    \\begin{cases} 
        \\frac{1}{\\mathcal{N}}\\mathcal{F}(s, \\mu) \\quad \\mathrm{use\\_windows == true} \\\\
        1 \\quad\\quad \\mathrm{use\\_windows == false}
    \\end{cases}
```

where:
- ``s_2 = s_2(s_1, s, \\mu) = \\sqrt{s_1^2 + s^2 + 2 \\, s_1 \\, s \\, \\mu}`` 
- ``y = \\cos{\\theta} = \\hat{\\mathbf{s}}_1 \\cdot \\hat{\\mathbf{s}}_2 = \\frac{\\mu \\, s + s_1}{s_2(s_1, s, \\mu)}``
- ``\\xi`` is the corresponding LD TPCF effect
- ``\\mathcal{L}_L(\\mu)`` is the Legendre polynomial of order ``L``
- ``\\mathcal{F}(s, \\mu)`` is the integrated window function stored in `cosmo::Cosmology` (check the documentation of `WindowFIntegrated`)
- ``\\mathcal{N}`` is the integrated window function norm (check the documentation of `WindowFIntegrated`) 

In the former method you have to pass as an input the LD `effect` function you want 
to integrate, while in the (recommended) latter one it's necessary to specify the
name of the LD TPCF term among the following:

`$(string(GR_EFFECTS_LD .* " , "...))`

to which correspond the following functions:

`$(string(string.(VEC_ξs_LD) .* " , "...))`

Internally, the second method recall the first, so all the keyword arguments are in common.

Note that these is an internal conversion of coordiate sistems
from `(s1, s, μ)` to `(s1, s2, y)` thorugh the functions `y` and `s2`. The inverse relations are:
- ``s = s(s_1, s_2, y) = \\sqrt{s_1^2 + s_2^2 - 2 \\, s_1 \\, s_2 \\, y}``
- ``\\mu = \\cos{\\theta} = \\hat{\\mathbf{s}}_1 \\cdot \\hat{\\mathbf{s}} = \\frac{y \\, s_2 - s_1}{s(s_1, s_2, y)}``.


## Inputs

- `s1`: the comoving distance where must be evaluated the integral

- `s`: the comoving distance from `s1` where must be evaluated the integral

- `effect::Function` or `effect::String` : name/function you want to consider, among the ones
  inside `GaPSE.GR_EFFECTS_LD`/`GaPSE.VEC_ξs_LD` respectively.

- `cosmo::Cosmology`: cosmology to be used in this computation


## Optional arguments

- `L::Int = 0`: order of the Legendre polynomial to be used

- `alg::Symbol = :trap` : algorithm to be used for the integration; the valid options 
  are (other values will lead to `AssertionError`):
  - `:quad` -> the integration over ``\\mu`` will be preformed through the Julia function `quadgk` 
  from the [`QuadGK.jl`](https://github.com/JuliaMath/QuadGK.jl) Julia package, that uses an adaptive 
  Gauss-Kronrod quadrature.
  - `:trap` -> the integration over ``\\mu`` will be preformed through the Julia function `trapz` 
  from the [`Trapz.jl`](https://github.com/francescoalemanno/Trapz.jl) Julia package, that uses the
  simple trapezoidal rulae.
  - `:lobatto` -> the integration over ``\\mu`` will be preformed through the Julia function `gausslobatto` 
  from the [`FastGaussQuadrature.jl`](https://github.com/JuliaApproximation/FastGaussQuadrature.jl) Julia package, 
  that uses the Gauss-Lobatto quadrature. 
  WE RECOMMEND TO USE `:quad` FOR MONOPOLES AND `:lobatto` FOR HIGHER ORDER MULTIPOLES!
  
- `use_windows::Bool = false`: tells if the integrand must consider the two
   window function ``\\phi`` and ``\\mathcal{F}``

- `N_lob::Int = 100` : number of points to be used in the sampling made by the function `trapz`.
  Note that these options will have an effect only if you se `alg = :quad`.

- `N_trap::Int = 200` : number of points to be used in the sampling made by the function `trapz`.
  Note that these options will have an effect only if you se `alg = :quad`.

- `atol_quad::Float64 = 0.0` and `rtol_quad::Float64 = 1e-2`: absolute and relative tolerance
  to be passed to the function `quadgk`; it's recommended not to set `rtol_quad < 1e-2` 
  because the time for evaluation increase quickly.
  Note that these options will have an effect only if you se `alg = :quad`.

- `enhancer::Float64 = 1e6`: just a float number used in order to deal better with small numbers; 
  the returned value is NOT modified by this value, because after a multiplication
  the internal result is divided by `enhancer`.

- `kwargs...` : other keyword arguments that will be passed to the selected 
  LD TPCF effect (`ξ_LD_Doppler`, `ξ_LD_Lensing`, ...)

See also: [`integrand_ξ_LD_multipole`](@ref), 
[`map_ξ_LD_multipole`](@ref), [`print_map_ξ_LD_multipole`](@ref)
[`WindowFIntegrated`](@ref), [`Cosmology`](@ref), 
[`y`](@ref), [`s2`](@ref), [`GR_EFFECTS_LD`](@ref), [`VEC_ξs_LD`](@ref)
"""
ξ_LD_multipole


##########################################################################################92


"""
    map_ξ_LD_multipole(cosmo::Cosmology,
        effect::Union{String,Function}, ss = nothing;
        s1 = nothing, L::Int = 0, alg::Symbol = :lobatto,
        N_lob::Int = 100, N_trap::Int = 50,
        atol_quad::Float64 = 0.0, rtol_quad::Float64 = 1e-2,
        enhancer::Float64=1e6, N_log::Int = 1000, 
        pr::Bool = true,
        kwargs...) ::Tuple{Vector{Float64}, Vector{Float64}}

Evaluate the multipole of order `L` of the chosen perturbation of the Luminosity 
Distance (LD) Two-Point Correlation Function (TPCF) term for all the comoving distance 
values stored inside `ss`.
If `ss = nothing`, it is set `ss = 10 .^ range(0, log10(2 * cosmo.s_max), length=N_log)`.
If `s1 = nothing`, it is set `s1 = cosmo.s_eff`.

The function evaluated is then the following ``\\xi_L(s_1, s)``:

```math
    \\xi_L(s_1, s) = \\frac{2 L + 1}{2} \\int_{-1}^{+1} \\mathrm{d}\\mu \\; 
        \\xi (s_1, s_2, \\cos{\\theta}) \\, \\mathcal{L}_L(\\mu) \\,  \\times
    \\begin{cases} 
        \\frac{1}{\\mathcal{N}}\\mathcal{F}(s, \\mu) \\quad \\mathrm{use\\_windows == true} \\\\
        1 \\quad\\quad \\mathrm{use\\_windows == false}
    \\end{cases}
```

where:
- ``s_2 = s_2(s_1, s, \\mu) = \\sqrt{s_1^2 + s^2 + 2 \\, s_1 \\, s \\, \\mu}`` 
- ``y = \\cos{\\theta} = \\hat{\\mathbf{s}}_1 \\cdot \\hat{\\mathbf{s}}_2 = \\frac{\\mu \\, s + s_1}{s_2(s_1, s, \\mu)}``
- ``\\xi`` is the corresponding LD TPCF effect
- ``\\mathcal{L}_L(\\mu)`` is the Legendre polynomial of order ``L``
- ``\\mathcal{F}(s, \\mu)`` is the integrated window function stored in `cosmo::Cosmology` (check the documentation of `WindowFIntegrated`)
- ``\\mathcal{N}`` is the integrated window function norm (check the documentation of `WindowFIntegrated`) 

Note that these is an internal conversion of coordiate sistems
from `(s1, s, μ)` to `(s1, s2, y)` thorugh the functions `y` and `s2`. The inverse relations are:
- ``s = s(s_1, s_2, y) = \\sqrt{s_1^2 + s_2^2 - 2 \\, s_1 \\, s_2 \\, y}``
- ``\\mu = \\cos{\\theta} = \\hat{\\mathbf{s}}_1 \\cdot \\hat{\\mathbf{s}} = \\frac{y \\, s_2 - s_1}{s(s_1, s_2, y)}``.


## Inputs

- `cosmo::Cosmology` : cosmology to be used in this computation

- `effect::Union{String,Function}` : the GR effect TPCF you want to consider; you may
  specify the name of the effect as one of the following strings (recommended):

  `$(string(GR_EFFECTS_LD .* " , "...))`
  
  or directly the name of the function among the following: 
  
  `$(string(string.(VEC_ξs_LD) .* " , "...))`

- `ss` : vector/range of `s` values where the function must be evaluated; if `ss = nothing`, 
  it is set `ss = 10 .^ range(0, log10(2 * cosmo.s_max), length=N_log)`. This is why it is returned 
  also the vector of the "input" values.


## Optional arguments

- `s1 = nothing` : comoving distance from the observer where the TPCF should be evaluated;
  if `s1 = nothing`, it is automatically set `s1 = cosmo.s_eff` from the given input `cosmo::Cosmology`.

- `L::Int = 0`: order of the Legendre polynomial to be used

- `alg::Symbol = :trap` : algorithm to be used for the integration; the valid options 
  are (other values will lead to `AssertionError`):
  - `:quad` -> the integration over ``\\mu`` will be preformed through the Julia function `quadgk` 
  from the [`QuadGK.jl`](https://github.com/JuliaMath/QuadGK.jl) Julia package, that uses an adaptive 
  Gauss-Kronrod quadrature.
  - `:trap` -> the integration over ``\\mu`` will be preformed through the Julia function `trapz` 
  from the [`Trapz.jl`](https://github.com/francescoalemanno/Trapz.jl) Julia package, that uses the
  simple trapezoidal rulae.
  - `:lobatto` -> the integration over ``\\mu`` will be preformed through the Julia function `gausslobatto` 
  from the [`FastGaussQuadrature.jl`](https://github.com/JuliaApproximation/FastGaussQuadrature.jl) Julia package, 
  that uses the Gauss-Lobatto quadrature. 
  WE RECOMMEND TO USE `:quad` FOR MONOPOLES AND `:lobatto` FOR HIGHER ORDER MULTIPOLES!
  
- `use_windows::Bool = false`: tells if the integrand must consider the two
   window function ``\\phi`` and ``\\mathcal{F}``

- `N_lob::Int = 100` : number of points to be used in the sampling made by the function `trapz`.
  Note that these options will have an effect only if you se `alg = :quad`.

- `N_trap::Int = 200` : number of points to be used in the sampling made by the function `trapz`.
  Note that these options will have an effect only if you se `alg = :quad`.

- `atol_quad::Float64 = 0.0` and `rtol_quad::Float64 = 1e-2`: absolute and relative tolerance
  to be passed to the function `quadgk`; it's recommended not to set `rtol_quad < 1e-2` 
  because the time for evaluation increase quickly.
  Note that these options will have an effect only if you se `alg = :quad`.

- `enhancer::Float64 = 1e6`: just a float number used in order to deal better with small numbers; 
  the returned value is NOT modified by this value, because after a multiplication
  the internal result is divided by `enhancer`.

- `N_log::Int = 1000` : number of points to be used in the default logaritmically-spaced 
  range for `ss`, i.e. `range(0, log10(2 * cosmo.s_max), length=N_log)`; it is ignored if `ss ≠ nothing` 

- `pr::Bool = true` : do you want the progress bar showed on screen, in order to 
  check the time needed for the computation? (`true` recommended)

- `kwargs...` : other keyword arguments that will be passed to the selected 
  LD TPCF effect (`ξ_LD_Doppler`, `ξ_LD_Lensing`, ...)

# Returns

A `Tuple{Vector{Float64}, Vector{Float64}}`, which has as first element the `ss` vector
and as second one the corresponding ξ value evaluated.

See also: [`integrand_ξ_LD_multipole`](@ref), [`ξ_LD_multipole`](@ref),
[`print_map_ξ_LD_multipole`](@ref),
[`WindowFIntegrated`](@ref), [`Cosmology`](@ref), 
[`y`](@ref), [`s2`](@ref), [`GR_EFFECTS_LD`](@ref), [`GaPSE.VEC_ξs_LD`](@ref)
"""
function map_ξ_LD_multipole(cosmo::Cosmology,
    effect::Union{String,Function}, ss=nothing;
    s1=nothing, L::Int=0, alg::Symbol=:lobatto,
    N_lob::Int=100, N_trap::Int=50,
    atol_quad::Float64=0.0, rtol_quad::Float64=1e-2,
    enhancer::Float64=1e6, N_log::Int=1000,
    pr::Bool=true, sum_xi::Bool=false,
    kwargs...)

    @assert alg ∈ VALID_INTEGRATION_ALGORITHM ":$alg is not a valid Symbol for \"alg\"; they are: \n\t" *
                                            "$(":".*string.(VALID_INTEGRATION_ALGORITHM) .* vcat([" , " for i in 1:length(VALID_INTEGRATION_ALGORITHM)-1], " .")... )"

    @assert N_trap > 2 "N_trap must be >2,  N_trap = $N_trap is not!"
    @assert N_lob > 2 "N_lob must be >2,  N_lob = $N_lob is not!"
    @assert atol_quad ≥ 0.0 "atol_quad must be ≥ 0.0,  atol_quad = $atol_quad is not!"
    @assert rtol_quad ≥ 0.0 "rtol_trap must be ≥ 0.0,  rtol_quad = $rtol_quad is not!"
    @assert L ≥ 0 "L must be ≥ 0, L = $L is not!"

    t1 = time()


    s_1 = isnothing(s1) ? cosmo.s_eff : s1
    v_ss = isnothing(ss) ? 10 .^ range(0, log10(2 * cosmo.s_max), length=N_log) : ss
    xis = [0.0 for i in 1:length(v_ss)]

    orig_f(μ, s) = enhancer * integrand_ξ_LD_multipole(s_1, s, μ, effect, cosmo;
        L=L, kwargs...)

    if alg == :lobatto
        μs, ws = gausslobatto(N_lob)

        global xis = pr ? begin
            @showprogress "$effect, L=$L: " [
                dot(ws, [orig_f(μ, s) for μ in μs]) / enhancer for s in v_ss
            ]
        end : [
            dot(ws, [orig_f(μ, s) for μ in μs]) / enhancer for s in v_ss
        ]

    elseif alg == :quad

        global xis = pr ? begin
            @showprogress "$effect, L=$L: " [
                quadgk(μ -> orig_f(μ, s), -1.0, 1.0;
                        atol=atol_quad, rtol=rtol_quad)[1] / enhancer for s in v_ss
            ]
        end : [
            quadgk(μ -> orig_f(μ, s), -1.0, 1.0;
                atol=atol_quad, rtol=rtol_quad)[1] / enhancer for s in v_ss
        ]

    elseif alg == :trap

        μs = union(
            range(-1.0, -0.98, length=Int(ceil(N_trap / 3) + 1)),
            range(-0.98, 0.98, length=Int(ceil(N_trap / 3) + 1)),
            range(0.98, 1.0, length=Int(ceil(N_trap / 3) + 1))
        )
        #μs = range(-1.0 + 1e-6, 1.0 - 1e-6, length=N_trap)

        global xis = pr ? begin
            @showprogress "$effect, L=$L: " [
                trapz(μs, [orig_f(μ, s) for μ in μs]) / enhancer for s in v_ss
            ]
        end : [
            trapz(μs, [orig_f(μ, s) for μ in μs]) / enhancer for s in v_ss
        ]

    else
        throw(AssertionError("how the hell did you arrive here?"))
    end

    #=
    xis = pr ? begin
        @showprogress "$effect, L=$L: " [
            ξ_LD_multipole(s_1, s, effect, cosmo; L=L, kwargs...) for s in v_ss
        ]
    end : [
        ξ_LD_multipole(s_1, s, effect, cosmo; L=L, kwargs...) for s in v_ss
    ]
    =#

    t2 = time()
    (!sum_xi) && (pr) && println("\ntime needed for map_ξ_LD_multipole for $effect, L=$L " *
                                "[in s] = $(@sprintf("%.5f", t2-t1)) \n")
    return (v_ss, xis)
end


##########################################################################################92


"""
    print_map_ξ_LD_multipole(
        cosmo::Cosmology, out::String,
        effect::Union{String,Function},
        ss = nothing;
        s1 = nothing, L::Int = 0, alg::Symbol = :lobatto,
        N_lob::Int = 100, N_trap::Int = 50,
        atol_quad::Float64 = 0.0, rtol_quad::Float64 = 1e-2,
        enhancer::Float64=1e6, N_log::Int = 1000, 
        pr::Bool = true,
        kwargs...)

Evaluate the multipole of order `L` of the of the chosen perturbation of the Luminosity 
Distance (LD) Two-Point Correlation Function (TPCF) term for all the comoving distance 
values stored inside `ss`, and print the results (with all the options used) 
in a file named `out`.
If `ss = nothing`, it is set `ss = 10 .^ range(0, log10(2 * cosmo.s_max), length=N_log)`.
If `s1 = nothing`, it is set `s1 = cosmo.s_eff`.

The function evaluated is then the following ``\\xi_L(s_1, s)``:

```math
    \\xi_L(s_1, s) = \\frac{2 L + 1}{2} \\int_{-1}^{+1} \\mathrm{d}\\mu \\; 
        \\xi (s_1, s_2, \\cos{\\theta}) \\, \\mathcal{L}_L(\\mu) \\,  \\times
    \\begin{cases} 
        \\frac{1}{\\mathcal{N}}\\mathcal{F}(s, \\mu) \\quad \\mathrm{use\\_windows == true} \\\\
        1 \\quad\\quad \\mathrm{use\\_windows == false}
    \\end{cases}
```

where:
- ``s_2 = s_2(s_1, s, \\mu) = \\sqrt{s_1^2 + s^2 + 2 \\, s_1 \\, s \\, \\mu}`` 
- ``y = \\cos{\\theta} = \\hat{\\mathbf{s}}_1 \\cdot \\hat{\\mathbf{s}}_2 = \\frac{\\mu \\, s + s_1}{s_2(s_1, s, \\mu)}``
- ``\\xi`` is the corresponding LD TPCF effect
- ``\\mathcal{L}_L(\\mu)`` is the Legendre polynomial of order ``L``
- ``\\mathcal{F}(s, \\mu)`` is the integrated window function stored in `cosmo::Cosmology` (check the documentation of `WindowFIntegrated`)
- ``\\mathcal{N}`` is the integrated window function norm (check the documentation of `WindowFIntegrated`) 

Note that these is an internal conversion of coordiate sistems
from `(s1, s, μ)` to `(s1, s2, y)` thorugh the functions `y` and `s2`. The inverse relations are:
- ``s = s(s_1, s_2, y) = \\sqrt{s_1^2 + s_2^2 - 2 \\, s_1 \\, s_2 \\, y}``
- ``\\mu = \\cos{\\theta} = \\hat{\\mathbf{s}}_1 \\cdot \\hat{\\mathbf{s}} = \\frac{y \\, s_2 - s_1}{s(s_1, s_2, y)}``.

## Inputs

- `cosmo::Cosmology` : cosmology to be used in this computation

- `effect::Union{String,Function}` : the GR effect TPCF you want to consider; you may
  specify the name of the effect as one of the following strings (recommended):

  `$(string(GR_EFFECTS_LD .* " , "...))`
  
  or directly the name of the function among the following: 
  
  `$(string(string.(VEC_ξs_LD) .* " , "...))`

- `out::String` : name of the file where the results must be stored.

- `ss` : vector/range of `s` values where the function must be evaluated; if `ss = nothing`, 
  it is set `ss = 10 .^ range(0, log10(2 * cosmo.s_max), length=N_log)`.


## Optional arguments

This function recall internally `map_ξ_LD_multipole`, so the kwargs are the same; we report them
for comfortness:

- `s1 = nothing` : comoving distance from the observer where the TPCF should be evaluated;
  if `s1 = nothing`, it is automatically set `s1 = cosmo.s_eff` from the given input `cosmo::Cosmology`.

- `L::Int = 0`: order of the Legendre polynomial to be used

- `alg::Symbol = :trap` : algorithm to be used for the integration; the valid options 
  are (other values will lead to `AssertionError`):
  - `:quad` -> the integration over ``\\mu`` will be preformed through the Julia function `quadgk` 
  from the [`QuadGK.jl`](https://github.com/JuliaMath/QuadGK.jl) Julia package, that uses an adaptive 
  Gauss-Kronrod quadrature.
  - `:trap` -> the integration over ``\\mu`` will be preformed through the Julia function `trapz` 
  from the [`Trapz.jl`](https://github.com/francescoalemanno/Trapz.jl) Julia package, that uses the
  simple trapezoidal rulae.
  - `:lobatto` -> the integration over ``\\mu`` will be preformed through the Julia function `gausslobatto` 
  from the [`FastGaussQuadrature.jl`](https://github.com/JuliaApproximation/FastGaussQuadrature.jl) Julia package, 
  that uses the Gauss-Lobatto quadrature. 
  WE RECOMMEND TO USE `:quad` FOR MONOPOLES AND `:lobatto` FOR HIGHER ORDER MULTIPOLES!
  
- `use_windows::Bool = false`: tells if the integrand must consider the two
   window function ``\\phi`` and ``\\mathcal{F}``

- `N_lob::Int = 100` : number of points to be used in the sampling made by the function `trapz`.
  Note that these options will have an effect only if you se `alg = :quad`.

- `N_trap::Int = 200` : number of points to be used in the sampling made by the function `trapz`.
  Note that these options will have an effect only if you se `alg = :quad`.

- `atol_quad::Float64 = 0.0` and `rtol_quad::Float64 = 1e-2`: absolute and relative tolerance
  to be passed to the function `quadgk`; it's recommended not to set `rtol_quad < 1e-2` 
  because the time for evaluation increase quickly.
  Note that these options will have an effect only if you se `alg = :quad`.

- `enhancer::Float64 = 1e6`: just a float number used in order to deal better with small numbers; 
  the returned value is NOT modified by this value, because after a multiplication
  the internal result is divided by `enhancer`.

- `N_log::Int = 1000` : number of points to be used in the default logaritmically-spaced 
  range for `ss`, i.e. `range(0, log10(2 * cosmo.s_max), length=N_log)`; it is ignored if `ss ≠ nothing` 

- `pr::Bool = true` : do you want the progress bar showed on screen, in order to 
  check the time needed for the computation? (`true` recommended)

- `kwargs...` : other keyword arguments that will be passed to the selected 
  LD TPCF effect (`ξ_LD_Doppler`, `ξ_LD_Lensing`, ...)

See also: [`integrand_ξ_LD_multipole`](@ref), [`ξ_LD_multipole`](@ref),
[`print_map_ξ_LD_multipole`](@ref),
[`WindowFIntegrated`](@ref), [`Cosmology`](@ref), 
[`y`](@ref), [`s2`](@ref), [`GR_EFFECTS_LD`](@ref), [`GaPSE.VEC_ξs_LD`](@ref)
"""
function print_map_ξ_LD_multipole(
    cosmo::Cosmology, out::String,
    effect::Union{String,Function}, ss=nothing;
    s1=nothing, L::Int=0, kwargs...)

    check_parent_directory(out)
    check_namefile(out)

    s_1 = isnothing(s1) ? cosmo.s_eff : s1
    t1 = time()
    vec = map_ξ_LD_multipole(cosmo, effect, ss; s1=s_1, L=L, kwargs...)
    t2 = time()

    isfile(out) && run(`rm $out`)
    open(out, "w") do io
        println(io, BRAND)

        println(io, "#\n# This is an integration map on mu of the ξ L=$L multipole $effect GR effect")
        println(io, "# concerning the relativistic galaxy number counts.")
        parameters_used(io, cosmo; logo=false)
        println(io, "# computational time needed (in s) : $(@sprintf("%.4f", t2-t1))")
        print(io, "# kwards passed: ")

        println(io, "\n# \t\tL = $L")
        if !isempty(kwargs)
            for key in keys(kwargs)
                println(io, "# \t\t$(key) = $(kwargs[key])")
            end
        end

        isnothing(s1) || println(io, "#\n# NOTE: the computation is done not in " *
                                    "s1 = s_eff, because you specified in input s1 = $s1 !")
        println(io, "# ")
        println(io, "# s [Mpc/h_0] \t \t xi")
        for (s, xi) in zip(vec[1], vec[2])
            println(io, "$s \t $xi")
        end
    end
end

