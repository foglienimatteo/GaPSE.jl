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
    ξ_GNCxLD_Newtonian_LocalGP(P1::Point, P2::Point, y, cosmo::Cosmology;
        b1=nothing, b2=nothing, s_b1=nothing, s_b2=nothing,
        𝑓_evo1=nothing, 𝑓_evo2=nothing, s_lim=nothing ) ::Float64

Return the LocalGP-LocalGP cross-correlation function concerning the perturbed
luminosity distance, defined as follows:

```math
\\xi^{v_{\\parallel}\\phi} (s_1, s_2, \\cos{\\theta}) = 
    \\frac{3}{2 a(s_2)} \\mathcal{H}(s_1) f(s_1) D(s_1)
    \\mathcal{R}(s_1) \\mathcal{H}_0^2 \\Omega_{M0} D(s_2)
    (1 + \\mathcal{R}(s_2)) (s_2\\cos{\\theta} - s_1) s^2 I^3_1(s)
```
where ``\\mathcal{H} = a H``,
``y = \\cos{\\theta} = \\hat{\\mathbf{s}}_1 \\cdot \\hat{\\mathbf{s}}_2`` and :

```math
I^n_l(s) = \\int_0^\\infty \\frac{\\mathrm{d}q}{2\\pi^2} q^2 \\, P(q) \\, \\frac{j_l(qs)}{(q s)^n}
```

## Inputs

- `P1::Point` and `P2::Point`: `Point` where the CF has to be calculated; they contain all the 
  data of interest needed for this calculus (comoving distance, growth factor and so on)
     
- `y`: the cosine of the angle between the two points `P1` and `P2`

- `cosmo::Cosmology`: cosmology to be used in this computation


See also: [`Point`](@ref), [`Cosmology`](@ref)
"""
function ξ_GNCxLD_Newtonian_LocalGP(P1::Point, P2::Point, y, cosmo::Cosmology;
    b1=nothing, b2=nothing, s_b1=nothing, s_b2=nothing,
    𝑓_evo1=nothing, 𝑓_evo2=nothing, s_lim=nothing)

    s1, D1, f1 = P1.comdist, P1.D, P1.f
    s2, D2, a2, ℜ2 = P2.comdist, P2.D, P2.a, P2.ℛ_LD

    Ω_M0 = cosmo.params.Ω_M0
    b1 = isnothing(b1) ? cosmo.params.b1 : b1

    Δs = s(s1, s2, y)

    common = - ℋ0^2 * Ω_M0 * (1 + ℜ2) / a2
    factor = f1 * ((3 * y^2 - 1) * s2^2 - 4 * y * s1 * s2 + 2 * s1^2)

    J20 = - 1 / 2 * (3 * b1 + f1) * (- 2 * y * s1 * s2 + s1^2 + s2^2)

    I00 = cosmo.tools.I00(Δs)
    I20 = cosmo.tools.I20(Δs)
    I40 = cosmo.tools.I40(Δs)
    I02 = cosmo.tools.I02(Δs)

    return D1 * D2 * common * (
        factor * (1 / 30 * I00 + 1 / 21 * I20 + 1 / 70 * I40) 
        + J20 * I02
        )
end


function ξ_GNCxLD_Newtonian_LocalGP(s1, s2, y, cosmo::Cosmology; kwargs...)
    P1, P2 = Point(s1, cosmo), Point(s2, cosmo)
    return ξ_GNCxLD_Newtonian_LocalGP(P1, P2, y, cosmo; kwargs...)
end



##########################################################################################92

##########################################################################################92

##########################################################################################92



function ξ_LDxGNC_LocalGP_Newtonian(s1, s2, y, cosmo::Cosmology; 
    b1=nothing, b2=nothing, s_b1=nothing, s_b2=nothing,
    𝑓_evo1=nothing, 𝑓_evo2=nothing, s_lim=nothing, kwargs...)

    b1 = isnothing(b1) ? cosmo.params.b1 : b1
    b2 = isnothing(b2) ? cosmo.params.b2 : b2
    s_b1 = isnothing(s_b1) ? cosmo.params.s_b1 : s_b1
    s_b2 = isnothing(s_b2) ? cosmo.params.s_b2 : s_b2
    𝑓_evo1 = isnothing(𝑓_evo1) ? cosmo.params.𝑓_evo1 : 𝑓_evo1
    𝑓_evo2 = isnothing(𝑓_evo2) ? cosmo.params.𝑓_evo2 : 𝑓_evo2

    ξ_GNCxLD_Newtonian_LocalGP(s2, s1, y, cosmo;
        b1=b2, b2=b1, s_b1=s_b2, s_b2=s_b1,
        𝑓_evo1=𝑓_evo2, 𝑓_evo2=𝑓_evo1, s_lim=s_lim, kwargs...)
end



