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

function ξ_PPGalaxies(s1, y, cosmo::Cosmology)
     return ξ_PPGalaxies_L0(s1, cosmo) + ξ_PPGalaxies_L2(s1, cosmo) * Pl(y, 2) +
            ξ_PPGalaxies_L4(s1, cosmo) * Pl(y, 4)
end



##########################################################################################92



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



function ξ_PPGalaxies_multipole(
     s, cosmo::Cosmology;
     L::Int=0,
     use_windows::Bool=true,
     enhancer::Float64=1e6,
     atol_quad::Float64=0.0,
     rtol_quad::Float64=1e-2)

     orig_f(μ) = enhancer * integrand_ξ_PPGalaxies_multipole(s, μ, cosmo;
          L=L, use_windows=use_windows)

     int = quadgk(μ -> orig_f(μ), -1.0, 1.0; atol=atol_quad, rtol=rtol_quad)[1]

     return int / enhancer
end


##########################################################################################92



function map_ξ_PPGalaxies_multipole(cosmo::Cosmology,
     v_ss=nothing;
     pr::Bool=true,
     N_log::Int=1000,
     L::Int=0,
     kwargs...)

     t1 = time()
     ss = isnothing(v_ss) ? 10 .^ range(0, 3, length=N_log) : v_ss
     xis = pr ? begin
          @showprogress "PP Galaxies, L=$L: " [
               ξ_PPGalaxies_multipole(s, cosmo; L=L, kwargs...) for s in ss
          ]
     end : [
          ξ_PPGalaxies_multipole(s, cosmo; L=L, kwargs...) for s in ss
     ]

     t2 = time()
     pr && println("\ntime needed for map_ξ_PPGalaxies_multipole " *
                   "[in s] = $(@sprintf("%.5f", t2-t1)) ")
     return (ss, xis)
end


##########################################################################################92



function print_map_ξ_PPGalaxies_multipole(
     cosmo::Cosmology,
     out::String,
     v_ss=nothing;
     L::Int=0,
     kwargs...)

     t1 = time()
     vec = map_ξ_PPGalaxies_multipole(cosmo, v_ss; L=L, kwargs...)
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




