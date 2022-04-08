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


function ξ_GNC_Newtonian(P1::Point, P2::Point, y, cosmo::Cosmology)
     s1, D1, f1 = P1.comdist, P1.D, P1.f
     s2, D2, f2 = P2.comdist, P2.D, P1.f
     b1, b2 = cosmo.params.b, cosmo.params.b

     Δs = s(s1, s2, y)

     J00 = 1.0/15.0 * (15.0 * b1 * b2 + 5.0 * (b1 * f2 + b2 * f1) + (2.0 * y^2 + 1.0) * f1 * f2)
     J02 = -1.0/(21.0*Δs^2) * (
          s1^2 * (14.0 * b2 * f1 + 7.0 * b1 * f2 * (3.0 * y^2 - 1.0) + (11.0 * y^2 + 1.0) * f1 * f2)
          + s2^2 * (14.0 * b1 * f2 + 7.0 * b2 * f1 * (3.0 * y^2 - 1.0) + (11.0 * y^2 + 1.0) * f1 * f2)
          - 4.0 * y * s1 * s2 * (7.0 * b2 * f1 + 7.0 * b1 * f2 + (y^2 + 5.0) * f1 * f2)
     )
     J04 = f1 * f2/(35.0 * Δs^4) * (
          4.0 * (3.0 * y^2 - 1.0) * (s1^4 + s2^4) 
          + 3.0 * (3.0 + y^2)^2 * s1^2 * s2^2
          - 8.0 * y * s1 * s2 * (s1^2 + s2^2) * (3.0 + y^2)
     )

     I00 = cosmo.tools.I00(Δs)
     I20 = cosmo.tools.I20(Δs)
     I40 = cosmo.tools.I40(Δs)

     res = D1 * D2 * (J00 * I00 + J02 * I20 + J04 * I40)

     return res
end


function ξ_GNC_Newtonian(s1, s2, y, cosmo::Cosmology)
     P1, P2 = Point(s1, cosmo), Point(s2, cosmo)
     return ξ_GNC_Newtonian(P1, P2, y, cosmo)
end


"""
     ξ_GNC_Newtonian(P1::Point, P2::Point, y, cosmo::Cosmology) :: Float64

     ξ_GNC_Newtonian(s1, s2, y, cosmo::Cosmology) = 
          ξ_GNC_Newtonian(Point(s1, cosmo), Point(s2, cosmo), y, cosmo::Cosmology)

Return the local gravitational potential auto-correlation function concerning the perturbed
luminosity distance, defined as follows:
```math
\\xi^{\\phi\\phi} (s_1, s_2, \\cos{\\theta}) = 
     \\frac{9 \\mathcal{H}_0^4 \\Omega_{M0}^2 D(s_1) D(s_2)s^4}{4 a(s_1) a(s_2)}
     (1 + \\mathcal{R}_1 + \\mathcal{R}_2 + \\mathcal{R}_1\\mathcal{R}_2)
     \\tilde{I}^4_0(s)
```
where ``D_1 = D(s_1)``, ``D_2 = D(s_2)`` and so on, ``\\mathcal{H} = a H``, 
``y = \\cos{\\theta} = \\hat{\\mathbf{s}}_1 \\cdot \\hat{\\mathbf{s}}_2`` and:
```math
\\tilde{I}^4_0 (s) &= \\int_0^\\infty \\frac{\\mathrm{d}q}{2\\pi^2} 
          q^2 \\, P(q) \\, \\frac{j_0(q s) - 1}{(q s)^4}
```

## Inputs

- `P1::Point` and `P2::Point`: `Point` where the CF has to be calculated; they contain all the 
  data of interest needed for this calculus (comoving distance, growth factor and so on)
     
- `y`: the cosine of the angle between the two points `P1` and `P2`

- `cosmo::Cosmology`: cosmology to be used in this computation

See also: [`Point`](@ref), [`Cosmology`](@ref)
"""
ξ_GNC_Newtonian
