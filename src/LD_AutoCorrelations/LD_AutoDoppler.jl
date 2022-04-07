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


function ξ_LD_Doppler(P1::Point, P2::Point, y, cosmo::Cosmology)
     s1, D1, f1, ℋ1, ℛ1 = P1.comdist, P1.D, P1.f, P1.ℋ, P1.ℛ_LD
     s2, D2, f2, ℋ2, ℛ2 = P2.comdist, P2.D, P2.f, P2.ℋ, P2.ℛ_LD

     Δs = s(P1.comdist, P2.comdist, y)
     prefac = D1 * D2 * f1 * f2 * ℛ1 * ℛ2 * ℋ1 * ℋ2
     c1 = 3 * s1 * s2 - 2 * y * (s1^2 + s2^2) + s1 * s2 * y^2

     I00 = cosmo.tools.I00(Δs)
     I20 = cosmo.tools.I20(Δs)
     I40 = cosmo.tools.I40(Δs)
     I02 = cosmo.tools.I02(Δs)

     parenth = I00 / 45.0 + I20 / 31.5 + I40 / 105.0

     first = prefac * (c1 * parenth + I02 * y * Δs^2 / 3.0)

     return first
end


function ξ_LD_Doppler(s1, s2, y, cosmo::Cosmology)
     P1, P2 = Point(s1, cosmo), Point(s2, cosmo)
     return ξ_LD_Doppler(P1, P2, y, cosmo)
end



"""
     ξ_LD_Doppler(P1::Point, P2::Point, y, cosmo::Cosmology) :: Float64

     ξ_LD_Doppler(s1, s2, y, cosmo::Cosmology) = 
          ξ_LD_Doppler(Point(s1, cosmo), Point(s2, cosmo), y, cosmo)

Return the Doppler auto-correlation function concerning the perturbed
luminosity distance, defined as follows:
```math
\\xi^{v_{\\parallel}v_{\\parallel}} (s_1, s_2, \\cos{\\theta}) 
= D_1 D_2 f_1 f_2 \\mathcal{H}_1 \\mathcal{H}_2 \\mathcal{R}_1 \\mathcal{R}_2 
(J_{00}\\, I^0_0(s) + J_{02}\\,I^0_2(s) + J_{04}\\,I^0_4(s) + J_{20}\\,I^2_0(s))
```
where ``D_1 = D(s_1)``, ``D_2 = D(s_2)`` and so on, ``\\mathcal{H} = a H``, 
``y = \\cos{\\theta} = \\hat{\\mathbf{s}}_1 \\cdot \\hat{\\mathbf{s}}_2`` and 
the J coefficients are given by:
```math
\\begin{align*}
    J_{00} (s_1, s_2, y) & = \\frac{1}{45} (y^2 s_1 s_2 - 2y(s_1^2 + s_2^2) + 3s_1 s_2) \\\\
    J_{02} (s_1, s_2, y) & = \\frac{2}{63} (y^2 s_1 s_2 - 2y(s_1^2 + s_2^2) + 3s_1 s_2) \\\\
    J_{04} (s_1, s_2, y) & = \\frac{1}{105} (y^2 s_1 s_2 - 2y(s_1^2 + s_2^2) + 3s_1 s_2) \\\\
    J_{20} (s_1, s_2, y) & = \\frac{1}{3} y s^2
\\end{align*}
```

## Inputs

- `P1::Point` and `P2::Point`: `Point` where the CF has to be calculated; they contain all the 
  data of interest needed for this calculus (comoving distance, growth factor and so on)
     
- `y`: the cosine of the angle between the two points `P1` and `P2`

- `cosmo::Cosmology`: cosmology to be used in this computation


See also: [`Point`](@ref), [`Cosmology`](@ref)
"""
ξ_LD_Doppler
