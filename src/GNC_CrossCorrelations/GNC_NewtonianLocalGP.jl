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


@doc raw"""
     ξ_GNC_Newtonian_LocalGP(P1::Point, P2::Point, y, cosmo::Cosmology) :: Float64

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
function ξ_GNC_Newtonian_LocalGP(P1::Point, P2::Point, y, cosmo::Cosmology)
     s1, D1, f1 = P1.comdist, P1.D, P1.f
     s2, D2, f2, a2, ℋ2, ℛ2 = P2.comdist, P2.D, P2.f, P2.a, P2.ℋ, P2.ℛ_GNC
     b1 = cosmo.params.b
     𝑓_evo2 = cosmo.params.𝑓_evo
     Ω_M0 = cosmo.params.Ω_M0

     Δs = s(s1, s2, y)

     common = 2 * f2 * a2 * ℋ2^2 * (𝑓_evo2 - 3) + 3 * ℋ0^2 * Ω_M0 * (f2 + ℛ2 + 5 * s_b2 - 2)
     factor = f1 * ((3 * y^2 - 1) * s2^2 - 4 * y * s1 * s2 + 2 * s1^2)

     J20 = - 1 / 6 * (3 * b1 + f1) * (- 2 * y * s1 * s2 + s1^2 + s2^2)

     I00 = cosmo.tools.I00(Δs)
     I20 = cosmo.tools.I20(Δs)
     I40 = cosmo.tools.I40(Δs)
     I02 = cosmo.tools.I02(Δs)

     return D1 * D2 / a2 * common * (
          factor * (1 / 90 * I00 + 1 / 63 * I20 + 1 / 210 * I40) 
          + J20 * I02
          )
end


function ξ_GNC_Newtonian_LocalGP(s1, s2, y, cosmo::Cosmology)
     P1, P2 = Point(s1, cosmo), Point(s2, cosmo)
     return ξ_GNC_Newtonian_LocalGP(P1, P2, y, cosmo)
end




##########################################################################################92

##########################################################################################92

##########################################################################################92



function ξ_GNC_LocalGP_Newtonian(s1, s2, y, cosmo::Cosmology)
     ξ_GNC_Newtonian_LocalGP(s2, s1, y, cosmo)
end

