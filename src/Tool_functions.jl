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


F_map_data = readdlm(FILE_F_MAP, comments = true)
F_map_data_dict = Dict([name => F_map_data[2:end, i] for (i, name) in enumerate(NAMES_F_MAP)]...)

_xs = unique(F_map_data_dict["x"])
_μs = unique(F_map_data_dict["mu"])
_Fs = F_map_data_dict["F"]


# for my F map convenction
my_F_grid = GridInterpolations.RectangleGrid(_μs, _xs)
spline_F(x, μ) = GridInterpolations.interpolate(my_F_grid, _Fs, [μ, x])

# for the opposite convenction for F map
#other_F_grid = GridInterpolations.RectangleGrid(_xs, _μs)
#spline_F(x, μ) = GridInterpolations.interpolate(other_F_grid, _Fs, [x, μ])



@doc raw"""
     spline_F(x, μ) :: Float64

Return the 2-dim spline value of F in the given `(x,μ)` input.
The spline is obtained through the `interpolate` function of the 
[`GridInterpolations`](https://github.com/sisl/GridInterpolations.jl) Julia
package applied to the input `FILE_F_MAP`.
"""
spline_F

##########################################################################################92



ps = readdlm(FILE_PS, comments = true)
ps_dict = Dict([name => ps[:, i] for (i, name) in enumerate(NAMES_PS)]...)

PK = Spline1D(ps_dict["k (h/Mpc)"], ps_dict["P (Mpc/h)^3"])

N = 1024                                # number of points to use in the Fourier transform
k_max = ps_dict["k (h/Mpc)"][end]       # maximum k-value
k_min = ps_dict["k (h/Mpc)"][begin]     # minimum k-value
s0 = 1 / k_max;                         # minimum r-value (should be ~1/k_max)

I00 = Spline1D(xicalc(PK, 0, 0; N = N, kmin = k_min, kmax = k_max, r0 = s0)...)
I20 = Spline1D(xicalc(PK, 2, 0; N = N, kmin = k_min, kmax = k_max, r0 = s0)...)
I40 = Spline1D(xicalc(PK, 4, 0; N = N, kmin = k_min, kmax = k_max, r0 = s0)...)
I02 = Spline1D(xicalc(PK, 0, 2; N = N, kmin = k_min, kmax = k_max, r0 = s0)...)
I22 = Spline1D(xicalc(PK, 2, 2; N = N, kmin = k_min, kmax = k_max, r0 = s0)...)
I31 = Spline1D(xicalc(PK, 3, 1; N = N, kmin = k_min, kmax = k_max, r0 = s0)...)
I13 = Spline1D(xicalc(PK, 1, 3; N = N, kmin = k_min, kmax = k_max, r0 = s0)...)
I11 = Spline1D(xicalc(PK, 1, 1; N = N, kmin = k_min, kmax = k_max, r0 = s0)...)


@doc raw"""
     I00, I20, I40, I02, I22, I31, I13, I11 ::Float64

Return the value of the integral:

```math
I_\ell^n(s)=\int_0^\infty \frac{\mathrm{d} q}{2 \pi^2} q^2 \, P(q) 
    \, \frac{j_\ell(qs)}{(qs)^n}
```

where, for a generic Iab name, ``\ell`` is the FIRST number (`a`) and 
``n`` the second (`b`).

These function are obtained through a `Spline1D` (from the 
[Dierckx](https://github.com/kbarbary/Dierckx.jl) Julia package) of the Spherical
Bessel Transform function `xicalc` (from the 
[TwoFAST](https://github.com/hsgg/TwoFAST.jl) Julia package) applied to the 
input Power Spectrum `P(q)`.
"""
I00, I20, I40, I02, I22, I31, I13, I11



##########################################################################################92



σ_0 = quadgk(q -> PK(q) * q^2 / (2 * π^2), k_min, k_max)[1]
σ_1 = quadgk(q -> PK(q) * q / (2 * π^2), k_min, k_max)[1]
σ_2 = quadgk(q -> PK(q) / (2 * π^2), k_min, k_max)[1]
σ_3 = quadgk(q -> PK(q) / (2 * π^2 * q), k_min, k_max)[1]


@doc raw"""
     σ_0, σ_1, σ_2, σ_3

These are the results of the following integral:
```math
     \sigma_i = \int_0^\infty \frac{\mathrm{d} q}{2 \pi^2} q^{2-i} \, P(q) 
```
where  `P(q)` is the input Power Spectrum.
"""
σ_0, σ_1, σ_2, σ_3



##########################################################################################92


@doc raw"""
     ϕ(s; s_min = s_min, s_max = s_max) :: Float64

Radial part of the survey window function. Return `1.0` if is true that
``s_\mathrm{min} \le s \le s_\mathrm{max}`` and `0.0` otherwise.

In this software we made the assuption that the survey window function can be
separated into a radial and angular part, i.e.:

```math
     \phi(\mathbf{s}) = \phi(s) \, W(\hat{s})
```

See also: [`W`](@ref)
"""
ϕ(s; s_min = s_min, s_max = s_max) = s_min < s < s_max ? 1.0 : 0.0



@doc raw"""
     W(θ; θ_max = θ_MAX) :: Float64

Angular part of the survey window function. Return `1.0` if is true that
``0.0 \leq \theta \le \theta_\mathrm{max}`` and `0.0` otherwise. It is
implicitly assumed an azimutal simmetry of the survey.

In this software we made the assuption that the survey window function can be
separated into a radial and angular part, i.e.:

```math
     \phi(\mathbf{s}) = \phi(s) \, W(\hat{s})
```

See also: [`ϕ`](@ref)
"""
W(θ; θ_max = θ_MAX) = 0.0 ≤ θ < θ_max ? 1.0 : 0.0


@doc raw"""
     V_survey(s_min = s_min, s_max = s_max, θ_max = θ_MAX) :: Float64

Return the volume of a survey with azimutal simmetry, i.e.:

```math
\begin{split}
    V(s_\mathrm{max}, s_\mathrm{min}, \theta_\mathrm{max}) &= \; C_\mathrm{up} - C_\mathrm{down} + TC \\
    &C_\mathrm{up} = \frac{\pi}{3} s_\mathrm{max}^3 \, 
        (1 - \cos\theta_\mathrm{max})^2 \, (2 + \cos\theta_\mathrm{max}) \\
    &C_\mathrm{down} = \frac{\pi}{3} s_\mathrm{min}^3 \, 
        (1 - \cos\theta_\mathrm{max})^2 \, (2 + \cos\theta_\mathrm{max}) \\
    &TC = \frac{\pi}{3} (s_\mathrm{max}^2 + s_\mathrm{min}^2 + 
        s_\mathrm{max} \,s_\mathrm{min}) \,  (s_\mathrm{max} - s_\mathrm{min})\, 
        \cos\theta_\mathrm{max}\, \sin^2\theta_\mathrm{max}
\end{split}
```
"""
function V_survey(s_min = s_min, s_max = s_max, θ_max = θ_MAX)
     sin_θ, cos_θ = sin(θ_max), cos(θ_max)
     diff_up_down = (s_max^3 - s_min^3) * (1 - cos_θ)^2 * (2 + cos_θ)
     tr = (s_max^2 + s_min^2 + s_max * s_min) * (s_max - s_min) * cos_θ * sin_θ^2
     #r1, r2 = s_min * sin(θ_max), s_max * sin(θ_max)
     #d1, d2 = s_min * cos(θ_max), s_max * cos(θ_max)
     #calotta_up = π / 3 * (s_max - d2)^2 * (2 * s_max + d2)
     #calotta_down = π / 3 * (s_min - d1)^2 * (2 * s_min + d1)
     #tronco_cono = π / 3 * (r1^2 + r1 * r2 + r2^2) * (s_max - s_min) * cos(θ_max)
     return π / 3.0 * (diff_up_down + tr)
end


@doc raw"""
     A(s_min = s_min, s_max = s_max, θ_max = θ_MAX) :: Float64

Return the Power Spectrum multipole normalization coefficient `A`, i.e.:
```math
     A(s_\mathrm{max}, s_\mathrm{min}, \theta_\mathrm{max})= 2 \, \pi \, 
     V(s_\mathrm{max}, s_\mathrm{min}, \theta_\mathrm{max})
```
where ``V(s_\mathrm{max}, s_\mathrm{min}, \theta_\mathrm{max})`` is the 
survey volume.

Pay attention: this is NOT used for the normalization of [`PS`](@ref), see
instead [`A_prime`](@ref)

See also: [`V_survey`](@ref)
"""
function A(s_min = s_min, s_max = s_max, θ_max = θ_MAX)
     2.0 * π * V_survey(s_min, s_max, θ_max)
end


@doc raw"""
     A_prime :: Float64

It's the Power Spectrum multipole normalization coefficient ``A^{'}``, i.e.:
```math
     A^{'} = \frac{3 \, A}{ (s_\mathrm{max}^3 - s_\mathrm{min}^3)} = 4 \pi^2
```

See also: [`A`](@ref), [`V_survey`](@ref)
"""
const A_prime = 4.0 * π^2

