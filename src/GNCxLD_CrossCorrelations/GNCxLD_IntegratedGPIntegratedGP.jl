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
    integrand_ξ_GNCxLD_IntegratedGP_IntegratedGP(IP1::Point, IP2::Point,
        P1::Point, P2::Point, y, cosmo::Cosmology;
        b1=nothing, b2=nothing, s_b1=nothing, s_b2=nothing, 
        𝑓_evo1=nothing, 𝑓_evo2=nothing, s_lim=nothing ) ::Float64

Return the integrand of the integrated gravitational potential 
auto-correlation function ``\\xi^{\\int\\phi\\int\\phi} (s_1, s_2, \\cos{\\theta})``, 
i.e. the function ``f(s_1, s_2, y, \\chi_1, \\chi_2)`` defined as follows:  

```math
f(s_1, s_2, y, \\chi_1, \\chi_2) = J_{40}(s_1, s_2, y, \\chi_1, \\chi_2) \\tilde{I}^4_0(\\chi)
```
where ``\\chi = \\sqrt{\\chi_1^2 + \\chi_2^2 - 2 \\, \\chi_1 \\, \\chi_2 \\, y} ``,
``y = \\cos{\\theta} = \\hat{\\mathbf{s}}_1 \\cdot \\hat{\\mathbf{s}}_2`` and:
```math
\\begin{split}
    &J_{40}(s_1, s_2, y, \\chi_1, \\chi_2)  = 
        \\frac{
                9 \\mathcal{H}_0^4 \\Omega_{M0}^2 D(\\chi_1) D(\\chi_2) \\chi^4
        }{    a(\\chi_1) a(\\chi_2) s_1 s_2} 
        (s_2 \\mathcal{H}(\\chi_2) \\mathcal{R}(s_2) (f(\\chi_2)-1) - 1) 
        (s_1 \\mathcal{H}(\\chi_1) \\mathcal{R}(s_1) (f(\\chi_1)-1) - 1)\\\\[5pt]
    &\\tilde{I}^4_0 (s) = \\int_0^\\infty \\frac{\\mathrm{d}q}{2\\pi^2} 
        q^2 \\, P(q) \\, \\frac{j_0(q s) - 1}{(q s)^4}
\\end{split}
```


## Inputs

- `IP1::Point` and `IP2::Point`: `Point` inside the integration limits, placed 
  at comoving distance `χ1` and `χ2` respectively.

- `P1::Point` and `P2::Point`: extreme `Point` of the integration, placed 
  at comoving distance `s1` and `s2` respectively.

- `y`: the cosine of the angle between the two points `P1` and `P2`

- `cosmo::Cosmology`: cosmology to be used in this computation


See also: [`ξ_GNCxLD_IntegratedGP_IntegratedGP`](@ref), [`integrand_on_mu_IntegratedGP`](@ref)
[`integral_on_mu`](@ref), [`ξ_GNC_multipole`](@ref)
"""
function integrand_ξ_GNCxLD_IntegratedGP_IntegratedGP(IP1::Point, IP2::Point,
    P1::Point, P2::Point, y, cosmo::Cosmology;
    b1=nothing, b2=nothing, s_b1=nothing, s_b2=nothing,
    𝑓_evo1=nothing, 𝑓_evo2=nothing, s_lim=nothing)

    s1 = P1.comdist
    s2, ℜ_s2 = P2.comdist, P2.ℛ_LD 
    χ1, D1, a1, ℋ1, f1 = IP1.comdist, IP1.D, IP1.a, IP1.ℋ, IP1.f
    χ2, D2, a2, ℋ2, f2 = IP2.comdist, IP2.D, IP2.a, IP2.ℋ, IP2.f

    Ω_M0 = cosmo.params.Ω_M0
    s_b_s1 = isnothing(s_b1) ? cosmo.params.s_b1 : s_b1
    𝑓_evo_s1 = isnothing(𝑓_evo1) ? cosmo.params.𝑓_evo1 : 𝑓_evo1

    s_lim = isnothing(s_lim) ? cosmo.params.s_lim : s_lim
    ℛ_s1 = func_ℛ_GNC(s1, P1.ℋ, P1.ℋ_p; s_b=s_b_s1, 𝑓_evo=𝑓_evo_s1, s_lim=s_lim)

    Δχ_square = χ1^2 + χ2^2 - 2 * χ1 * χ2 * y
    Δχ = Δχ_square > 0 ? √(Δχ_square) : 0

    factor = - 9 * Δχ^4 * ℋ0^4 * Ω_M0^2 * D1 * D2 / (s1 * s2 * a1 * a2)
    parenth_1 = s1 * ℋ1 * ℛ_s1 * (f1 - 1) - 5 * s_b_s1 + 2
    parenth_2 = s2 * ℋ2 * ℜ_s2 * (f2 - 1) - 1

    I04_tilde = cosmo.tools.I04_tilde(Δχ)

    return factor * parenth_1 * parenth_2 * I04_tilde
end



function ξ_GNCxLD_IntegratedGP_IntegratedGP(P1::Point, P2::Point, y, cosmo::Cosmology;
    en::Float64 = 1e10, N_χs_2::Int = 100, kwargs...)

    χ1s = P1.comdist .* range(1e-6, 1.0, length = N_χs_2)
    χ2s = P2.comdist .* range(1e-6, 1.0, length = N_χs_2 + 7)

    IP1s = [GaPSE.Point(x, cosmo) for x in χ1s]
    IP2s = [GaPSE.Point(x, cosmo) for x in χ2s]

    int_ξ_igp = [
        en * GaPSE.integrand_ξ_GNCxLD_IntegratedGP_IntegratedGP(IP1, IP2, P1, P2, y, cosmo; kwargs...)
        for IP1 in IP1s, IP2 in IP2s
    ]

    res = trapz((χ1s, χ2s), int_ξ_igp)
    #println("res = $res")

    #=
    χ1s = [x for x in range(0, P1.comdist, length = N_χs)[begin+1:end]]
    l = Int(floor(N_χs/2))
    matrix_χ2s = [begin
        a = [x for x in range(0, P2.comdist, length=l)[begin+1:end]];
        b = [x for x in range(x1-focus, x1+focus, length=l)];
        vcat(a[a.<x1-focus], b, a[a.>x1+focus])
        end for x1 in χ1s]

    IP1s = [GaPSE.Point(x, cosmo) for x in χ1s]
    matrix_IP2s = [[GaPSE.Point(x, cosmo) for x in y] for y in matrix_χ2s]
    matrix_int_ξs = [
        [en * GaPSE.integrand_ξ_GNCxLD_IntegratedGP_IntegratedGP(IP1, IP2, P1, P2, y, cosmo) 
        for IP2 in matrix_IP2s[i]]
        for (i,IP1) in enumerate(IP1s)]
    
    vec_trapz = [trapz(χ2s,int_ξs) for (χ2s,int_ξs) in zip(matrix_χ2s, matrix_int_ξs)]
    res = trapz(χ1s, vec_trapz)
    =#
    return res / en
end


function ξ_GNCxLD_IntegratedGP_IntegratedGP(s1, s2, y, cosmo::Cosmology; kwargs...)
    P1, P2 = Point(s1, cosmo), Point(s2, cosmo)
    return ξ_GNCxLD_IntegratedGP_IntegratedGP(P1, P2, y, cosmo; kwargs...)
end



"""
    ξ_GNCxLD_IntegratedGP_IntegratedGP(P1::Point, P2::Point, y, cosmo::Cosmology; 
        en::Float64 = 1e10, N_χs::Int = 100,
        b1=nothing, b2=nothing, s_b1=nothing, s_b2=nothing, 
        𝑓_evo1=nothing, 𝑓_evo2=nothing, s_lim=nothing ) ::Float64

    ξ_GNCxLD_IntegratedGP_IntegratedGP(s1, s2, y, cosmo::Cosmology; kwargs...)

Return the integrated gravitational potential auto-correlation function 
``\\xi^{\\int\\phi\\int\\phi}(s_1, s_2, \\cos{\\theta})`` concerning the perturbed
luminosity distance, defined as follows:
    
```math
\\xi^{\\int\\phi\\int\\phi} (s_1, s_2, \\cos{\\theta}) = 
    \\int_0^{s_1} \\mathrm{d} \\chi_1 \\int_0^{s_2}\\mathrm{d} \\chi_2 \\;
    J_{40}(s_1, s_2, y, \\chi_1, \\chi_2) \\, \\tilde{I}^4_0(\\chi)
```
where ``\\chi = \\sqrt{\\chi_1^2 + \\chi_2^2 - 2 \\, \\chi_1 \\, \\chi_2 \\, y} ``,
``y = \\cos{\\theta} = \\hat{\\mathbf{s}}_1 \\cdot \\hat{\\mathbf{s}}_2`` and:
```math
\\begin{split}
    &J_{40}(s_1, s_2, y, \\chi_1, \\chi_2)  = 
        \\frac{
            9 \\mathcal{H}_0^4 \\Omega_{M0}^2 D(\\chi_1) D(\\chi_2) \\chi^4
        }{    a(\\chi_1) a(\\chi_2) s_1 s_2} 
        (s_2 \\mathcal{H}(\\chi_2) \\mathcal{R}(s_2) (f(\\chi_2)-1) - 1) 
        (s_1 \\mathcal{H}(\\chi_1) \\mathcal{R}(s_1) (f(\\chi_1)-1) - 1)\\\\[5pt]
    &\\tilde{I}^4_0 (s) = \\int_0^\\infty \\frac{\\mathrm{d}q}{2\\pi^2} 
        q^2 \\, P(q) \\, \\frac{j_0(q s) - 1}{(q s)^4}
\\end{split}
```
and ``P(q)`` is the input power spectrum.


The computation is made applying [`trapz`](@ref) (see the 
[Trapz](https://github.com/francescoalemanno/Trapz.jl) Julia package) to
the integrand function `integrand_ξ_GNC_Lensing`.


## Inputs

- `s1` and `s2`: comovign distances where the function must be evaluated

- `y`: the cosine of the angle between the two points `P1` and `P2`

- `cosmo::Cosmology`: cosmology to be used in this computation


## Optional arguments 

- `en::Float64 = 1e10`: just a float number used in order to deal better 
  with small numbers.

- `N_χs::Int = 100`: number of points to be used for sampling the integral
  along the ranges `(0, s1)` (for `χ1`) and `(0, s1)` (for `χ2`); it has been checked that
  with `N_χs ≥ 50` the result is stable.

See also: [`integrand_ξ_GNCxLD_IntegratedGP_IntegratedGP`](@ref), [`integrand_on_mu_IntegratedGP`](@ref)
[`integral_on_mu`](@ref), [`ξ_GNC_multipole`](@ref)
"""
ξ_GNCxLD_IntegratedGP_IntegratedGP


##########################################################################################92

##########################################################################################92

##########################################################################################92



function ξ_LDxGNC_IntegratedGP_IntegratedGP(s1, s2, y, cosmo::Cosmology; 
        b1=nothing, b2=nothing, s_b1=nothing, s_b2=nothing, 
        𝑓_evo1=nothing, 𝑓_evo2=nothing, s_lim=nothing, kwargs...)
    
    b1 = isnothing(b1) ? cosmo.params.b1 : b1
    b2 = isnothing(b2) ? cosmo.params.b2 : b2
    s_b1 = isnothing(s_b1) ? cosmo.params.s_b1 : s_b1
    s_b2 = isnothing(s_b2) ? cosmo.params.s_b2 : s_b2
    𝑓_evo1 = isnothing(𝑓_evo1) ? cosmo.params.𝑓_evo1 : 𝑓_evo1
    𝑓_evo2 = isnothing(𝑓_evo2) ? cosmo.params.𝑓_evo2 : 𝑓_evo2

    ξ_GNCxLD_IntegratedGP_IntegratedGP(s2, s1, y, cosmo; 
        b1=b2, b2=b1, s_b1=s_b2, s_b2=s_b1,
        𝑓_evo1=𝑓_evo2, 𝑓_evo2=𝑓_evo1, s_lim=s_lim, kwargs...)
end
