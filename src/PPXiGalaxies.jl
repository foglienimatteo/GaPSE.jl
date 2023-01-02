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


function ξ_PPGalaxies_L0(P::Point, cosmo::Cosmology)
     b = cosmo.params.b
     s = P.comdist

     Peff = Point(cosmo.s_eff, cosmo)
     D, f = Peff.D, Peff.f

     (b^2 + 2 * b * f / 3 + f^2 / 5) * D^2 * cosmo.tools.I00(s)
end

function ξ_PPGalaxies_L0(s, cosmo::Cosmology)
     P = Point(s, cosmo)
     return ξ_PPGalaxies_L0(P, cosmo)
end


"""
     ξ_PPGalaxies_L0(P::Point, cosmo::Cosmology)
     ξ_PPGalaxies_L0(s, cosmo::Cosmology)

Return the value of the Two-Point Correlation Function (TPCF) monopole of the Galaxies
in the Plane-Parallel approximation.
In the first method, you should pass the `Point` where to evaluate that function,
while in the second (that internally recalls the first) you must provide the 
comoving distance `s`.

The analytical expression of such TPCF monopole is the following:
```math
\\xi^{\\mathrm{pp}}_0(s,z) = D^2(z) \\, I_0^0(s)\\, \\left(b^2 + 
\\frac{2}{3} \\, b \\, f(z) + \\frac{1}{5} \\, f^2(z)\\right)
```
where ``b`` is the galaxy bias, ``z`` the redshift, ``D`` the linear growth factor, 
``f`` the linear growth rate and ``I_\\ell^n`` is defined as
```math
I_\\ell^n(s) = \\int_0^{+\\infty} \\frac{\\mathrm{d}q}{2\\pi^2} 
\\, q^2 \\, P(q) \\, \\frac{j_\\ell(qs)}{(qs)^n}
```
with ``P(q)`` as the matter Power Spectrum at ``z=0`` and ``j_\\ell`` as spherical
Bessel function of order ``\\ell``.

All the cosmological data needed for this computation are taken from the input struct `cosmo::Cosmology`.

See also: [`Point`](@ref), [`Cosmology`](@ref)
"""
ξ_PPGalaxies_L0

function ξ_PPGalaxies_L2(P::Point, cosmo::Cosmology)
     b = cosmo.params.b
     s = P.comdist

     Peff = Point(cosmo.s_eff, cosmo)
     D, f = Peff.D, Peff.f

     -(4 * b * f / 3 + 4 * f^2 / 7) * D^2 * cosmo.tools.I20(s)
end

function ξ_PPGalaxies_L2(s, cosmo::Cosmology)
     P = Point(s, cosmo)
     return ξ_PPGalaxies_L2(P, cosmo)
end

"""
     ξ_PPGalaxies_L2(P::Point, cosmo::Cosmology)
     ξ_PPGalaxies_L2(s, cosmo::Cosmology)

Return the value of the Two-Point Correlation Function (TPCF) quadrupole of the Galaxies
in the Plane-Parallel approximation.
In the first method, you should pass the `Point` where to evaluate that function,
while in the second (that internally recalls the first) you must provide the 
comoving distance `s`.

The analytical expression of such TPCF monopole is the following:
```math
\\xi^{\\mathrm{pp}}_2(s,z) = - D^2(z) \\, I_2^0(s)\\, \\left(\\frac{4}{3} 
\\, b \\, f(z) + \\frac{4}{7} \\, f^2(z)\\right)
```
where ``b`` is the galaxy bias, ``z`` the redshift, ``D`` the linear growth factor, 
``f`` the linear growth rate and ``I_\\ell^n`` is defined as
```math
I_\\ell^n(s) = \\int_0^{+\\infty} \\frac{\\mathrm{d}q}{2\\pi^2} 
\\, q^2 \\, P(q) \\, \\frac{j_\\ell(qs)}{(qs)^n}
```
with ``P(q)`` as the matter Power Spectrum at ``z=0`` and ``j_\\ell`` as spherical
Bessel function of order ``\\ell``.

All the cosmological data needed for this computation are taken from the input struct `cosmo::Cosmology`.

See also: [`Point`](@ref), [`Cosmology`](@ref)
"""
ξ_PPGalaxies_L2

function ξ_PPGalaxies_L4(P::Point, cosmo::Cosmology)
     s = P.comdist

     Peff = Point(cosmo.s_eff, cosmo)
     D, f = Peff.D, Peff.f

     8 * f^2 / 35 * D^2 * cosmo.tools.I40(s)
end

function ξ_PPGalaxies_L4(s, cosmo::Cosmology)
     P = Point(s, cosmo)
     return ξ_PPGalaxies_L4(P, cosmo)
end


"""
     ξ_PPGalaxies_L4(P::Point, cosmo::Cosmology)
     ξ_PPGalaxies_L4(s, cosmo::Cosmology)

Return the value of the Two-Point Correlation Function (TPCF) hexadecapole of the Galaxies
in the Plane-Parallel approximation.
In the first method, you should pass the `Point` where to evaluate that function,
while in the second (that internally recalls the first) you must provide the 
comoving distance `s`.

The analytical expression of such TPCF monopole is the following:
```math
\\xi^{\\mathrm{pp}}_4(s,z) = D^2(z) \\, I_4^0(s)\\, 
\\left(\\frac{8}{35} \\, f^2(z)\\right)
```
where ``b`` is the galaxy bias, ``z`` the redshift, ``D`` the linear growth factor, 
``f`` the linear growth rate and ``I_\\ell^n`` is defined as
```math
I_\\ell^n(s) = \\int_0^{+\\infty} \\frac{\\mathrm{d}q}{2\\pi^2} 
\\, q^2 \\, P(q) \\, \\frac{j_\\ell(qs)}{(qs)^n}
```
with ``P(q)`` as the matter Power Spectrum at ``z=0`` and ``j_\\ell`` as spherical
Bessel function of order ``\\ell``.

All the cosmological data needed for this computation are taken from the input struct `cosmo::Cosmology`.

See also: [`Point`](@ref), [`Cosmology`](@ref)
"""
ξ_PPGalaxies_L4




"""
     ξ_PPGalaxies(s1, μ, cosmo::Cosmology)

Return the value of the Two-Point Correlation Function (TPCF) of the Galaxies
in the Plane-Parallel approximation in the given comoving distance `s` and cosine
value for the Legendre polynomials `μ` (for the given `cosmo::Cosmology`).

The analytical expression of such TPCF monopole is the following:
```math
\\begin{split}
\\xi^{\\mathrm{pp}}(s,&z,\\mu) = \\xi^{\\mathrm{pp}}_0(s,z) + 
    \\xi^{\\mathrm{pp}}_2(s,z) \\mathcal{L}_2(\\mu) + 
    \\xi^{\\mathrm{pp}}_4(s,z) \\mathcal{L}_4(\\mu)  \\\\
&\\xi^{\\mathrm{pp}}_0(s,z) = D^2(z) \\, I_0^0(s)\\, \\left(b^2 + 
    \\frac{2}{3} \\, b \\, f(z) + \\frac{1}{5} \\, f^2(z)\\right)  \\\\
&\\xi^{\\mathrm{pp}}_2(s,z) = - D^2(z) \\, I_2^0(s)\\, \\left(\\frac{4}{3} 
    \\, b \\, f(z) + \\frac{4}{7} \\, f^2(z)\\right)  \\\\
&\\xi^{\\mathrm{pp}}_4(s,z) = D^2(z) \\, I_4^0(s)\\, 
    \\left(\\frac{8}{35} \\, f^2(z)\\right)
\\end{split}
```
where ``b`` is the galaxy bias, ``z`` the redshift, ``D`` the linear growth factor, 
``f`` the linear growth rate, ``\\mathcal{L}_\\ell`` the Legendre polynomial of order ``\\ell`` 
and ``I_\\ell^n`` is defined as
```math
I_\\ell^n(s) = \\int_0^{+\\infty} \\frac{\\mathrm{d}q}{2\\pi^2} 
\\, q^2 \\, P(q) \\, \\frac{j_\\ell(qs)}{(qs)^n}
```
with ``P(q)`` as the matter Power Spectrum at ``z=0`` and ``j_\\ell`` as spherical
Bessel function of order ``\\ell``.

All the cosmological data needed for this computation are taken from the input struct `cosmo::Cosmology`.

See also: [`Point`](@ref), [`Cosmology`](@ref),
[`integrand_ξ_PPGalaxies_multipole`](@ref), [`ξ_PPGalaxies_multipole`](@ref) 
[`map_ξ_PPGalaxies_multipole`](@ref), [`print_map_ξ_PPGalaxies_multipole`](@ref)
"""
function ξ_PPGalaxies(s1, y, cosmo::Cosmology)
     return ξ_PPGalaxies_L0(s1, cosmo) + ξ_PPGalaxies_L2(s1, cosmo) * Pl(y, 2) +
            ξ_PPGalaxies_L4(s1, cosmo) * Pl(y, 4)
end



##########################################################################################92


"""
     integrand_ξ_PPGalaxies_multipole(s, μ, cosmo::Cosmology;
          L::Int=0, use_windows::Bool=true)

Return the integrand on ``\\mu = \\hat{\\mathbf{s}}_1 \\cdot \\hat{\\mathbf{s}}`` 
of the of the Galaxies Two-Point Correlation Function (TPCF) in the Plane Parallel (PP) 
approximation, i.e. the following function ``f(s, \\mu)``:

```math
     f_L(s, \\mu) = \\xi^{\\mathrm{pp}} \\left(s, \\mu\\right) 
          \\, \\mathcal{L}_L(\\mu) \\, \\times 
    \\begin{cases} 
        \\frac{1}{\\mathcal{N}}\\mathcal{F}(s, \\mu) \\quad \\mathrm{use\\_windows == true} \\\\
        1 \\quad\\quad \\mathrm{use\\_windows == false}
    \\end{cases}
```

where:
- ``\\xi`` is the Two-Point Correlation Function (TPCF) of the Galaxies in the Plane-Parallel approximation,
  computed from `ξ_PPGalaxies`.
- ``\\mathcal{L}_L(\\mu)`` is the Legendre polynomial of order ``L``
- ``\\mathcal{F}(s, \\mu)`` is the integrated window function stored in `cosmo::Cosmology` (check the documentation of `WindowFIntegrated`)
- ``\\mathcal{N}`` is the integrated window function norm (check the documentation of `WindowFIntegrated`)

## Inputs

- `s`: the comoving distance  where must be evaluated the integral

- `μ`: the cosine between `s1` and `s` where must be evaluated the integral

- `cosmo::Cosmology`: cosmology to be used in this computation

## Optional arguments 

- `L::Int = 0`: order of the Legendre polynomial to be used

- `use_windows::Bool = false`: tells if the integrand must consider the two
   window function ``\\phi`` and ``\\mathcal{F}``

See also:[`ξ_PPGalaxies`](@ref), [`ξ_PPGalaxies_multipole`](@ref), 
[`map_ξ_PPGalaxies_multipole`](@ref), [`print_map_ξ_PPGalaxies_multipole`](@ref)
[`WindowFIntegrated`](@ref), [`Cosmology`](@ref), 
"""
function integrand_ξ_PPGalaxies_multipole(s, μ, cosmo::Cosmology;
     L::Int=0, use_windows::Bool=true)

     res = if use_windows == true
          #println("s1 = $s1 \\t s2 = $(s2(s1, s, μ)) \\t  y=$(y(s1, s, μ))")
          #int = cosmo.ξ_matter(s)
          int = ξ_PPGalaxies(s, μ, cosmo)
          #println("int = $int")
          #coef = (cosmo.params.b + cosmo.f_of_s(s) * μ^2)^2 * cosmo.D_of_s(s)
          #println("coef = $coef")
          int .* (spline_integrF(s, μ, cosmo.windowFint) / cosmo.WFI_norm * Pl(μ, L))
     else
          #int = cosmo.ξ_matter(s)
          int = ξ_PPGalaxies(s, μ, cosmo)
          #coef = (cosmo.params.b + cosmo.f_of_s(s) * μ^2)^2
          int .* Pl(μ, L)
     end

     #println("res = $res")
     return (2.0 * L + 1.0) / 2.0 * res
end



##########################################################################################92



"""
     ξ_PPGalaxies_multipole(
          s, cosmo::Cosmology;
          L::Int = 0, use_windows::Bool = true,
          atol_quad::Float64 = 0.0,
          rtol_quad::Float64 = 1e-2
          enhancer::Float64 = 1e6 ) ::Float64


Evaluate the multipole of order `L` of the Galaxies Two-Point Correlation Function (TPCF) in the Plane 
Parallel (PP)  term i.e. the following function ``\\xi^{\\mathrm{pp}} (s)``:

```math
     \\xi^{\\mathrm{pp}} (s) = \\frac{2 L + 1}{2} \\int_{-1}^{+1} \\mathrm{d}\\mu \\; 
    \\xi^{\\mathrm{pp}} \\left(s, \\mu\\right) 
          \\, \\mathcal{L}_L(\\mu) \\, \\times 
    \\begin{cases} 
        \\frac{1}{\\mathcal{N}}\\mathcal{F}(s, \\mu) \\quad \\mathrm{use\\_windows == true} \\\\
        1 \\quad\\quad \\mathrm{use\\_windows == false}
    \\end{cases}
```

where:
- ``\\xi`` is the Two-Point Correlation Function (TPCF) of the Galaxies in the Plane-Parallel approximation,
  computed from `ξ_PPGalaxies`.
- ``\\mathcal{L}_L(\\mu)`` is the Legendre polynomial of order ``L``
- ``\\mathcal{F}(s, \\mu)`` is the integrated window function stored in `cosmo::Cosmology` (check the documentation of `WindowFIntegrated`)
- ``\\mathcal{N}`` is the integrated window function norm (check the documentation of `WindowFIntegrated`)

The integration over ``\\mu`` is preformed through the Julia function `quadgk` 
from the [`QuadGK.jl`](https://github.com/JuliaMath/QuadGK.jl) Julia package, that uses an adaptive 
Gauss-Kronrod quadrature.

## Inputs

- `s`: the comoving distance  where must be evaluated the integral

- `cosmo::Cosmology`: cosmology to be used in this computation

## Optional arguments 

- `L::Int = 0`: order of the Legendre polynomial to be used

- `use_windows::Bool = false`: tells if the integrand must consider the two
   window function ``\\phi`` and ``\\mathcal{F}``

- `atol_quad::Float64 = 0.0` and `rtol_quad::Float64 = 1e-2`: absolute and relative tolerance
  to be passed to the function `quadgk`; it's recommended not to set `rtol_quad < 1e-2` 
  because the time for evaluation increase quickly.

- `enhancer::Float64 = 1e6`: just a float number used in order to deal better with small numbers; 
  the returned value is NOT modified by this value, because after a multiplication
  the internal result is divided by `enhancer`.

See also: [`ξ_PPGalaxies`](@ref), [`integrand_ξ_PPGalaxies_multipole`](@ref), 
[`map_ξ_PPGalaxies_multipole`](@ref), [`print_map_ξ_PPGalaxies_multipole`](@ref)
[`WindowFIntegrated`](@ref), [`Cosmology`](@ref), 
"""
function ξ_PPGalaxies_multipole(
     s, cosmo::Cosmology;
     L::Int=0, use_windows::Bool=true,
     enhancer::Float64=1e6,
     atol_quad::Float64=0.0,
     rtol_quad::Float64=1e-2)

     orig_f(μ) = enhancer * integrand_ξ_PPGalaxies_multipole(s, μ, cosmo;
          L=L, use_windows=use_windows)

     int = quadgk(μ -> orig_f(μ), -1.0, 1.0; atol=atol_quad, rtol=rtol_quad)[1]

     return int / enhancer
end


##########################################################################################92


"""
     map_ξ_PPGalaxies_multipole(
          cosmo::Cosmology, ss = nothing;
          L::Int = 0, use_windows::Bool = true,
          atol_quad::Float64 = 0.0,
          rtol_quad::Float64 = 1e-2
          enhancer::Float64 = 1e6
          pr::Bool = true,
          N_log::Int = 1000,
          kwargs...) ::Tuple{Vector{Float64}, Vector{Float64}}


Evaluate the multipole of order `L` of the Galaxies Two-Point Correlation Function (TPCF) in the Plane 
Parallel (PP) term for all the comoving distance 
values stored inside `ss`.
If `ss = nothing`, it is set `ss = 10 .^ range(0, log10(2 * cosmo.s_max), length=N_log)`.

The function evaluated is then the following ``\\xi^{\\mathrm{pp}} (s)``:

```math
     \\xi^{\\mathrm{pp}} (s) = \\frac{2 L + 1}{2} \\int_{-1}^{+1} \\mathrm{d}\\mu \\; 
    \\xi^{\\mathrm{pp}} \\left(s, \\mu\\right) 
          \\, \\mathcal{L}_L(\\mu) \\, \\times 
    \\begin{cases} 
        \\frac{1}{\\mathcal{N}}\\mathcal{F}(s, \\mu) \\quad \\mathrm{use\\_windows == true} \\\\
        1 \\quad\\quad \\mathrm{use\\_windows == false}
    \\end{cases}
```

where:
- ``\\xi`` is the Two-Point Correlation Function (TPCF) of the Galaxies in the Plane-Parallel approximation,
  computed from `ξ_PPGalaxies`.
- ``\\mathcal{L}_L(\\mu)`` is the Legendre polynomial of order ``L``
- ``\\mathcal{F}(s, \\mu)`` is the integrated window function stored in `cosmo::Cosmology` (check the documentation of `WindowFIntegrated`)
- ``\\mathcal{N}`` is the integrated window function norm (check the documentation of `WindowFIntegrated`)

The integration over ``\\mu`` is preformed through the Julia function `quadgk` 
from the [`QuadGK.jl`](https://github.com/JuliaMath/QuadGK.jl) Julia package, that uses an adaptive 
Gauss-Kronrod quadrature.

## Inputs

- `cosmo::Cosmology`: cosmology to be used in this computation

- `ss` : vector/range of `s` values where the function must be evaluated; if `ss = nothing`, 
  it is set `ss = 10 .^ range(0, log10(2 * cosmo.s_max), length=N_log)`. This is why it is returned 
  also the vector of the "input" values.

## Optional arguments 

This function recall internally `ξ_PPGalaxies_multipole`, so the kwargs of the latter are valid also for the former; 
we report them for comfortness:

- `L::Int = 0`: order of the Legendre polynomial to be used

- `use_windows::Bool = false`: tells if the integrand must consider the two
   window function ``\\phi`` and ``\\mathcal{F}``

- `atol_quad::Float64 = 0.0` and `rtol_quad::Float64 = 1e-2`: absolute and relative tolerance
  to be passed to the function `quadgk`; it's recommended not to set `rtol_quad < 1e-2` 
  because the time for evaluation increase quickly.

- `enhancer::Float64 = 1e6`: just a float number used in order to deal better with small numbers; 
  the returned value is NOT modified by this value, because after a multiplication
  the internal result is divided by `enhancer`.

- `N_log::Int = 1000` : number of points to be used in the default logaritmically-spaced 
  range for `ss`, i.e. `range(0, log10(2 * cosmo.s_max), length=N_log)`; it is ignored if `ss ≠ nothing` 

- `pr::Bool = true` : do you want the progress bar showed on screen, in order to 
  check the time needed for the computation? (`true` recommended)

# Returns

A `Tuple{Vector{Float64}, Vector{Float64}}`, which has as first element the `ss` vector
and as second one the corresponding ξ value evaluated.

See also: [`ξ_PPGalaxies`](@ref), [`integrand_ξ_PPGalaxies_multipole`](@ref), 
[`ξ_PPGalaxies_multipole`](@ref), [`print_map_ξ_PPGalaxies_multipole`](@ref)
[`WindowFIntegrated`](@ref), [`Cosmology`](@ref), 
"""
function map_ξ_PPGalaxies_multipole(
     cosmo::Cosmology, ss=nothing;
     L::Int=0, pr::Bool=true,
     N_log::Int=1000,
     kwargs...)

     t1 = time()
     v_ss = isnothing(ss) ? 10 .^ range(0, log10(2 * cosmo.s_max), length=N_log) : ss
     xis = pr ? begin
          @showprogress "PP Galaxies, L=$L: " [
               ξ_PPGalaxies_multipole(s, cosmo; L=L, kwargs...) for s in v_ss
          ]
     end : [
          ξ_PPGalaxies_multipole(s, cosmo; L=L, kwargs...) for s in v_ss
     ]

     t2 = time()
     pr && println("\ntime needed for map_ξ_PPGalaxies_multipole " *
                   "[in s] = $(@sprintf("%.5f", t2-t1)) ")
     return (v_ss, xis)
end


##########################################################################################92


"""
     print_map_ξ_PPGalaxies_multipole(
         cosmo::Cosmology, out::String,
            ss = nothing;
         L::Int = 0, use_windows::Bool = true,
          atol_quad::Float64 = 0.0,
          rtol_quad::Float64 = 1e-2
          enhancer::Float64 = 1e6
        pr::Bool = true,
         N_log::Int = 1000,
         kwargs...)


Evaluate the multipole of order `L` of the Galaxies Two-Point Correlation Function (TPCF) in the Plane 
Parallel (PP) term for all the comoving distance values stored inside `ss`, 
and print the results (with all the options used) in a file named `out`.
If `ss = nothing`, it is set `ss = 10 .^ range(0, log10(2 * cosmo.s_max), length=N_log)`.

The function evaluated is then the following ``\\xi^{\\mathrm{pp}} (s)``:

```math
     \\xi^{\\mathrm{pp}} (s) = \\frac{2 L + 1}{2} \\int_{-1}^{+1} \\mathrm{d}\\mu \\; 
    \\xi^{\\mathrm{pp}} \\left(s, \\mu\\right) 
          \\, \\mathcal{L}_L(\\mu) \\, \\times 
    \\begin{cases} 
        \\frac{1}{\\mathcal{N}}\\mathcal{F}(s, \\mu) \\quad \\mathrm{use\\_windows == true} \\\\
        1 \\quad\\quad \\mathrm{use\\_windows == false}
    \\end{cases}
```

where:
- ``\\xi`` is the Two-Point Correlation Function (TPCF) of the Galaxies in the Plane-Parallel approximation,
  computed from `ξ_PPGalaxies`.
- ``\\mathcal{L}_L(\\mu)`` is the Legendre polynomial of order ``L``
- ``\\mathcal{F}(s, \\mu)`` is the integrated window function stored in `cosmo::Cosmology` (check the documentation of `WindowFIntegrated`)
- ``\\mathcal{N}`` is the integrated window function norm (check the documentation of `WindowFIntegrated`)

The integration over ``\\mu`` is preformed through the Julia function `quadgk` 
from the [`QuadGK.jl`](https://github.com/JuliaMath/QuadGK.jl) Julia package, that uses an adaptive 
Gauss-Kronrod quadrature.

## Inputs

- `cosmo::Cosmology`: cosmology to be used in this computation

- `out::String` : name of the file where the results must be stored.

- `ss` : vector/range of `s` values where the function must be evaluated; if `ss = nothing`, 
  it is set `ss = 10 .^ range(0, log10(2 * cosmo.s_max), length=N_log)`. This is why it is returned 
  also the vector of the "input" values.

## Optional arguments 

This function recall internally `map_ξ_PPGalaxies_multipole`, so the kwargs of the latter are valid also for the former; 
we report them for comfortness:

- `L::Int = 0`: order of the Legendre polynomial to be used

- `use_windows::Bool = false`: tells if the integrand must consider the two
   window function ``\\phi`` and ``\\mathcal{F}``

- `atol_quad::Float64 = 0.0` and `rtol_quad::Float64 = 1e-2`: absolute and relative tolerance
  to be passed to the function `quadgk`; it's recommended not to set `rtol_quad < 1e-2` 
  because the time for evaluation increase quickly.

- `enhancer::Float64 = 1e6`: just a float number used in order to deal better with small numbers; 
  the returned value is NOT modified by this value, because after a multiplication
  the internal result is divided by `enhancer`.

- `N_log::Int = 1000` : number of points to be used in the default logaritmically-spaced 
  range for `ss`, i.e. `range(0, log10(2 * cosmo.s_max), length=N_log)`; it is ignored if `ss ≠ nothing` 

- `pr::Bool = true` : do you want the progress bar showed on screen, in order to 
  check the time needed for the computation? (`true` recommended)

See also: [`ξ_PPGalaxies`](@ref), [`integrand_ξ_PPGalaxies_multipole`](@ref), 
[`ξ_PPGalaxies_multipole`](@ref), [`map_ξ_PPGalaxies_multipole`](@ref)
[`WindowFIntegrated`](@ref), [`Cosmology`](@ref), 
"""
function print_map_ξ_PPGalaxies_multipole(
     cosmo::Cosmology,
     out::String,
     ss=nothing;
     L::Int=0,
     kwargs...)

     t1 = time()
     vec = map_ξ_PPGalaxies_multipole(cosmo, ss; L=L, kwargs...)
     t2 = time()

     isfile(out) && run(`rm $out`)
     open(out, "w") do io
          println(io, BRAND)

          println(io, "# This is an integration map on mu of the ξ L=$L multipole concerning the PP Galaxies.")
          parameters_used(io, cosmo; logo=false)
          println(io, "# computational time needed (in s) : $(@sprintf("%.4f", t2-t1))")
          print(io, "# kwards passed: ")

          println(io, "\n# \t\tL = $L")
          if !isempty(kwargs)
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




