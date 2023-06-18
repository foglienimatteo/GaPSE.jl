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
	integrand_ξ_GNCxLD_Lensing_IntegratedGP(
		IP1::Point, IP2::Point, P1::Point, P2::Point,
		y, cosmo::Cosmology;
        b1=nothing, b2=nothing, s_b1=nothing, s_b2=nothing,
        𝑓_evo1=nothing, 𝑓_evo2=nothing, s_lim=nothing ) ::Float64

Return the integrand of theLensing-IntegratedGP cross-correlation function 
``\\xi^{\\kappa\\int\\phi} (s_1, s_2, \\cos{\\theta})``, i.e. the function 
``f(s_1, s_2, y, \\chi_1, \\chi_2)`` defined as follows:  

```math
f(s_1, s_2, y, \\chi_1, \\chi_2) = 
    \\frac{9}{2}\\mathcal{H}_0^4\\Omega_{M0}^2 
    \\frac{D(\\chi_1)D(\\chi_2)\\chi_2(s_1 - \\chi_1)}{s_1a(\\chi_1)a(\\chi_2)}
    \\left(\\mathcal{H}(\\chi_2)(f(\\chi_2) - 1)\\mathcal{R}(s_2) -\\frac{1}{s_2} \\right)
    \\left( J_{31} I^3_1(\\chi) + J_{22} I^2_2(\\chi) \\right)
```

where ``D_1 = D(\\chi_1)``, ``D_2 = D(\\chi_2)`` and so on, ``\\mathcal{H} = a H``, 
``\\chi = \\sqrt{\\chi_1^2 + \\chi_2^2 - 2\\chi_1\\chi_2\\cos{\\theta}}``, 
``y = \\cos{\\theta} = \\hat{\\mathbf{s}}_1 \\cdot \\hat{\\mathbf{s}}_2``) 
and the ``J`` coefficients are given by 

```math
\\begin{align*}
    J_{31} & = -2y\\chi^2 \\\\
    J_{22} & = \\chi_1\\chi_2(1-y^2)
\\end{align*}
```

## Inputs

- `IP1::Point` and `IP2::Point`: `Point` inside the integration limits, placed 
  at comoving distance `χ1` and `χ2` respectively.

- `P1::Point` and `P2::Point`: extreme `Point` of the integration, placed 
  at comoving distance `s1` and `s2` respectively.

- `y`: the cosine of the angle between the two points `P1` and `P2`

- `cosmo::Cosmology`: cosmology to be used in this computation


See also: [`ξ_GNCxLD_Lensing_IntegratedGP`](@ref), [`integrand_on_mu_Lensing_IntegratedGP`](@ref)
[`integral_on_mu`](@ref), [`ξ_LD_multipole`](@ref)
"""
function integrand_ξ_GNCxLD_Lensing_IntegratedGP(
    IP1::Point, IP2::Point, P1::Point, P2::Point, y, cosmo::Cosmology;
    b1=nothing, b2=nothing, s_b1=nothing, s_b2=nothing,
    𝑓_evo1=nothing, 𝑓_evo2=nothing, s_lim=nothing)

    s1 = P1.comdist
    s2, ℜ_s2 = P2.comdist, P2.ℛ_LD
    χ1, D1, a1 = IP1.comdist, IP1.D, IP1.a
    χ2, D2, a2, f2, ℋ2 = IP2.comdist, IP2.D, IP2.a, IP2.f, IP2.ℋ
 
    Ω_M0 = cosmo.params.Ω_M0
    s_b_s1 = isnothing(s_b1) ? cosmo.params.s_b1 : s_b1

    Δχ_square = χ1^2 + χ2^2 - 2 * χ1 * χ2 * y
    Δχ = √(Δχ_square) > 1e-8 ? √(Δχ_square) : 1e-8

    prefactor = 9 / 2 * ℋ0^4 * Ω_M0^2
    factor = D1 * D2 * χ2 * (s1 - χ1) * (5 * s_b_s1 - 2)/ (s1 * s2 * a1 * a2)
    parenth = s2 * ℋ2 * ℜ_s2 * (f2 - 1) - 1 

    new_J31 = 2 * y * Δχ^2
    new_J22 = χ1 * χ2 * (y^2 - 1)

    I13 = cosmo.tools.I13(Δχ)
    I22 = cosmo.tools.I22(Δχ)

    res = prefactor * factor * parenth * (new_J22 * I22 + new_J31 * I13)

    return res
end


function integrand_ξ_GNCxLD_Lensing_IntegratedGP(
    χ1::Float64, χ2::Float64,
    s1::Float64, s2::Float64,
    y, cosmo::Cosmology;
    kwargs...)

    P1, P2 = Point(s1, cosmo), Point(s2, cosmo)
    IP1, IP2 = Point(χ1, cosmo), Point(χ2, cosmo)
    return integrand_ξ_GNCxLD_Lensing_IntegratedGP(IP1, IP2, P1, P2, y, cosmo; kwargs...)
end


"""
    ξ_GNCxLD_Lensing_IntegratedGP(s1, s2, y, cosmo::Cosmology;
        en::Float64 = 1e6, N_χs::Int = 100,
        b1=nothing, b2=nothing, s_b1=nothing, s_b2=nothing,
        𝑓_evo1=nothing, 𝑓_evo2=nothing, s_lim=nothing ) ::Float64

Return theLensing-IntegratedGP cross-correlation function 
``\\xi^{\\kappa\\int\\phi} (s_1, s_2, \\cos{\\theta})`` concerning the perturbed
luminosity distance, defined as follows:
    
```math
\\xi^{\\kappa\\int\\phi} (s_1, s_2, \\cos{\\theta}) = 
     \\frac{9}{2}\\mathcal{H}_0^4\\Omega_{M0}^2 
     &\\mathrm{d}\\chi_1 \\int_0^{s_2} \\mathrm{d}\\chi_2 
     \\frac{D(\\chi_1)D(\\chi_2)\\chi_2(s_1 - \\chi_1)}{s_1a(\\chi_1)a(\\chi_2)} \\\\
     &\\left(\\mathcal{H}(\\chi_2)(f(\\chi_2) - 1)\\mathcal{R}(s_2) -\\frac{1}{s_2} \\right)
     \\left( J_{31} I^3_1(\\chi) + J_{22} I^2_2(\\chi) \\right)
\\end{split}
```

where ``D_1 = D(\\chi_1)``, ``D_2 = D(\\chi_2)`` and so on, ``\\mathcal{H} = a H``, 
``\\chi = \\sqrt{\\chi_1^2 + \\chi_2^2 - 2\\chi_1\\chi_2\\cos{\\theta}}``, 
``y = \\cos{\\theta} = \\hat{\\mathbf{s}}_1 \\cdot \\hat{\\mathbf{s}}_2``) 
and the ``J`` coefficients are given by 

```math
\\begin{align*}
    J_{31} & = -2y\\chi^2 \\\\
    J_{22} & = \\chi_1\\chi_2(1-y^2)
\\end{align*}
```

The computation is made applying [`trapz`](@ref) (see the 
[Trapz](https://github.com/francescoalemanno/Trapz.jl) Julia package) to
the integrand function `integrand_ξ_GNCxLD_Lensing_IntegratedGP`.



## Inputs

- `s1` and `s2`: comovign distances where the function must be evaluated

- `y`: the cosine of the angle between the two points `P1` and `P2`

- `cosmo::Cosmology`: cosmology to be used in this computation


## Optional arguments 

- `en::Float64 = 1e6`: just a float number used in order to deal better 
  with small numbers;

- `N_χs::Int = 100`: number of points to be used for sampling the integral
  along the ranges `(0, s1)` (for `χ1`) and `(0, s1)` (for `χ2`); it has been checked that
  with `N_χs ≥ 50` the result is stable.


See also: [`integrand_ξ_GNCxLD_Lensing_IntegratedGP`](@ref), [`integrand_on_mu_Lensing_IntegratedGP`](@ref)
[`integral_on_mu`](@ref), [`ξ_LD_multipole`](@ref)
"""
function ξ_GNCxLD_Lensing_IntegratedGP(
    P1::Point, P2::Point, y, cosmo::Cosmology;
    en::Float64 = 1e6, N_χs_2::Int = 100, kwargs...)

    χ1s = P1.comdist .* range(1e-6, 1.0, length = N_χs_2)
    χ2s = P2.comdist .* range(1e-6, 1.0, length = N_χs_2 + 7)

    IP1s = [GaPSE.Point(x, cosmo) for x in χ1s]
    IP2s = [GaPSE.Point(x, cosmo) for x in χ2s]

    int_ξs = [
        en * GaPSE.integrand_ξ_GNCxLD_Lensing_IntegratedGP(IP1, IP2, P1, P2, y, cosmo; kwargs...)
        for IP1 in IP1s, IP2 in IP2s
    ]

    res = trapz((χ1s, χ2s), int_ξs)
    #println("res = $res")
    return res / en
end


function ξ_GNCxLD_Lensing_IntegratedGP(s1, s2, y, cosmo::Cosmology; kwargs...)
    P1, P2 = Point(s1, cosmo), Point(s2, cosmo)
    return ξ_GNCxLD_Lensing_IntegratedGP(P1, P2, y, cosmo; kwargs...)
end


##########################################################################################92

##########################################################################################92

##########################################################################################92



function ξ_LDxGNC_IntegratedGP_Lensing(s1, s2, y, cosmo::Cosmology; 
    b1=nothing, b2=nothing, s_b1=nothing, s_b2=nothing,
    𝑓_evo1=nothing, 𝑓_evo2=nothing, s_lim=nothing, kwargs...)
    
    b1 = isnothing(b1) ? cosmo.params.b1 : b1
    b2 = isnothing(b2) ? cosmo.params.b2 : b2
    s_b1 = isnothing(s_b1) ? cosmo.params.s_b1 : s_b1
    s_b2 = isnothing(s_b2) ? cosmo.params.s_b2 : s_b2
    𝑓_evo1 = isnothing(𝑓_evo1) ? cosmo.params.𝑓_evo1 : 𝑓_evo1
    𝑓_evo2 = isnothing(𝑓_evo2) ? cosmo.params.𝑓_evo2 : 𝑓_evo2

    ξ_GNCxLD_Lensing_IntegratedGP(s2, s1, y, cosmo; 
        b1=b2, b2=b1, s_b1=s_b2, s_b2=s_b1,
        𝑓_evo1=𝑓_evo2, 𝑓_evo2=𝑓_evo1, s_lim=s_lim, kwargs...)
end
