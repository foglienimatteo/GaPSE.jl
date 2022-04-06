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


function my_Plm(l::Integer, m::Integer, x)
     if m ≥ 0
          return Plm(l, m, x)
     else
          return (-1)^m * factorial(l - m) / factorial(l + m) * Plm(l, -m, x)
     end
end


"""
     Aabcd(θ_max, a::Integer, b::Integer, c::Integer, d::Integer; kwargs...)

Calculates the following expression:

```math
\\begin{align*}
     A_{a b}^{c d}  
     &= \\int \\mathrm{d}\\theta \\sin \\left( \\theta \\right) W \\left( \\theta \\right) 
          P_a^c \\left( \\cos \\theta \\right)  P_b^d \\left( \\cos \\theta \\right) \\
     &= \\int_0^{\\theta_{\\mathrm{max}} \\mathrm{d}\\theta \\sin \\left( \\theta \\right) 
          P_a^c \\left( \\cos \\theta \\right)  P_b^d \\left( \\cos \\theta \\right) \\
     &= \\int_{\\cos \\theta_{\\mathrm{max}}}^1 \\mathrm{d}x 
          P_a^c \\left( x \\right)  P_b^d \\left( x \\right)
\\end{align*}
```

The integral is performed thanks to `quadgk`; all the `kwargs...` are passed to that function.
"""
function Aabcd(θ_max, a::Integer, b::Integer, c::Integer, d::Integer; kwargs...)
     quadgk(x -> my_Plm(a, c, x) * my_Plm(b, d, x), cos(θ_max), 1; kwargs...)[1]
end


"""
     Wllnn(θ_max, n1::Integer, n2::Integer, l3::Integer, l4::Integer; kwargs...)

Calculates the following function:

```math
     W_{n_1 n_2}^{\\ell_3 \\ell_4}  = 
          16 \\pi^4 \\sum_{m_3 m_4} 
          \\frac{\\left( n_1 -m_3 \\right)! }{ \\left(n_1 +m_3 \\right)!} 
          \\frac{\\left( \\ell_3 -m_3 \\right)! }{ \\left(\\ell_3 +m_3 \\right)!}
          \\frac{\\left( n_2 -m_4 \\right)! }{ \\left(n_2 +m_4 \\right)!} 
          \\frac{\\left( \\ell_4 -m_4 \\right)! }{ \\left(\\ell_4 +m_4 \\right)!} 
          A_{n_1 \\ell_3}^{m_3 m_3} A_{n_2 \\ell_3}^{m_4 m_3} 
          A_{\\ell_4 n_1 }^{m_4 m_3} A_{n_2 \\ell_4}^{m_4 m_4}
```
where 
```math
     A_{a b}^{c d} = 
          \\int_{\\cos \\theta_{\\mathrm{max}}}^1 \\mathrm{d}x 
          P_a^c \\left( x \\right)  P_b^d \\left( x \\right)
```
The ``A_{a b}^{c d}`` are computed through `Aabcd`.

See also: [`Aabcd`](@ref)
"""
function Wllnn(θ_max, n1::Integer, n2::Integer, l3::Integer, l4::Integer; rtol=1e-2, kwargs...)
     16 * π^4 * sum(
          [
          begin
               #println("n1 = $(n1) \t l3 = $(l3) \t m3 = $(m3)");
               #println("n1-m3 = $(n1-m3) \t l3-m3 = $(l3-m3) \t n2-m4 = $(n2-m4) \t l4-m4 = $(l4-m4)");
               #println("n1+m3 = $(n1+m3) \t l3+m3 = $(l3+m3) \t n2+m4 = $(n2+m4) \t l4+m4 = $(l4+m4)");
               factorial(n1 - m3) / factorial(n1 + m3) * factorial(l3 - m3) / factorial(l3 + m3) *
               factorial(n2 - m4) / factorial(n2 + m4) * factorial(l4 - m4) / factorial(l4 + m4)
               Aabcd(θ_max, n1, l3, m3, m3; rtol=rtol, kwargs...)
               Aabcd(θ_max, n2, l3, m4, m3; rtol=rtol, kwargs...) *
               Aabcd(θ_max, l4, n1, m4, m3; rtol=rtol, kwargs...) *
               Aabcd(θ_max, n2, l4, m4, m4; rtol=rtol, kwargs...)
          end
          for m3 in -l3:l3, m4 in -l4:l4
     ]
     )
end
