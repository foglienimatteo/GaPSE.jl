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
    WindowF(
        xs::Vector{Float64}
        μs::Vector{Float64}
        Fs::Matrix{Float64}
        )

Struct containing xs, μs and Fs values of the window function ``F(x, μ)``.
`xs` and `μs` are 1D vectors containing each value only once, while 
Fs values are contained in a matrix of size `(length(xs), length(μs))`, so:
- along a fixed column the changing value is `x`
- along a fixed row the changing value is `μ`
"""
struct WindowF
     xs::Vector{Float64}
     μs::Vector{Float64}
     Fs::Matrix{Float64}


     function WindowF(file::String)
          data = readdlm(file, comments = true)
          xs, μs, Fs = data[:, 1], data[:, 2], data[:, 3]
          @assert size(xs) == size(μs) == size(Fs) "xs, μs and Fs must have the same length!"

          new_xs = unique(xs)
          new_μs = unique(μs)
          new_Fs =
               if xs[2] == xs[1] && μs[2] ≠ μs[1]
                    transpose(reshape(Fs, (length(new_μs), length(new_xs))))
               elseif xs[2] ≠ xs[1] && μs[2] == μs[1]
                    reshape(Fs, (length(new_xs), length(new_μs)))
               else
                    throw(ErrorException("What kind of convenction for the file $file" *
                                         " are you using? I do not recognise it."))
               end
          new(new_xs, new_μs, new_Fs)
     end
end


@doc raw"""
     spline_F(x, μ, str::WindowF)) :: Float64

Return the 2-dim spline value of ``F`` in the given `(x,μ)`, where
``F`` is defined in the input `WindowF`.
The spline is obtained through the `interpolate` function of the 
[`GridInterpolations`](https://github.com/sisl/GridInterpolations.jl) Julia
package.
"""
function spline_F(x, μ, str::WindowF)
     grid = GridInterpolations.RectangleGrid(str.xs, str.μs)
     GridInterpolations.interpolate(grid, reshape(str.Fs, (:,1)), [x, μ])
end

#=
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
=#



##########################################################################################92




struct InputPS
     ks::Vector{Float64}
     pks::Vector{Float64}

     function InputPS(file::String)
          data = readdlm(file, comments = true)
          ks, pks = (data[:, 1], data[:, 2])
          @assert size(ks) == size(pks) "ks and pks must have the same length!"
          new(ks, pks)
     end

     function InputPS(ks::AbstractVector{T1}, pks::AbstractVector{T2}) where {T1,T2}
          @assert size(ks) == size(pks) "ks and pks must have the same length!"
          new(ks, pks)
     end
end



struct IPSTools
     I00::Dierckx.Spline1D
     I20::Dierckx.Spline1D
     I40::Dierckx.Spline1D
     I02::Dierckx.Spline1D
     I22::Dierckx.Spline1D
     I31::Dierckx.Spline1D
     I13::Dierckx.Spline1D
     I11::Dierckx.Spline1D

     σ_0::Float64
     σ_1::Float64
     σ_2::Float64
     σ_3::Float64

     function IPSTools(
          ips::InputPS;
          N = 1024,
          k_min::Union{Float64,Nothing} = nothing,
          k_max::Union{Float64,Nothing} = nothing,
          s_0::Union{Float64,Nothing} = nothing
     )

          PK = Spline1D(ips.ks, ips.pks)

          kmin = isnothing(k_min) ? min(ips.ks) : k_min
          kmax = isnothing(k_max) ? max(ips.ks) : k_max
          s0 = isnothing(s_0) ? 1.0 / kmax : s_0

          I00 = Spline1D(xicalc(PK, 0, 0; N = N, kmin = kmin, kmax = kmax, r0 = s0)...)
          I20 = Spline1D(xicalc(PK, 2, 0; N = N, kmin = kmin, kmax = kmax, r0 = s0)...)
          I40 = Spline1D(xicalc(PK, 4, 0; N = N, kmin = kmin, kmax = kmax, r0 = s0)...)
          I02 = Spline1D(xicalc(PK, 0, 2; N = N, kmin = kmin, kmax = kmax, r0 = s0)...)
          I22 = Spline1D(xicalc(PK, 2, 2; N = N, kmin = kmin, kmax = kmax, r0 = s0)...)
          I31 = Spline1D(xicalc(PK, 3, 1; N = N, kmin = kmin, kmax = kmax, r0 = s0)...)
          I13 = Spline1D(xicalc(PK, 1, 3; N = N, kmin = kmin, kmax = kmax, r0 = s0)...)
          I11 = Spline1D(xicalc(PK, 1, 1; N = N, kmin = kmin, kmax = kmax, r0 = s0)...)

          σ_0 = quadgk(q -> PK(q) * q^2 / (2 * π^2), kmin, kmax)[1]
          σ_1 = quadgk(q -> PK(q) * q / (2 * π^2), kmin, kmax)[1]
          σ_2 = quadgk(q -> PK(q) / (2 * π^2), kmin, kmax)[1]
          σ_3 = quadgk(q -> PK(q) / (2 * π^2 * q), kmin, kmax)[1]

          new(I00, I20, I40, I02, I22, I31, I13, I11, σ_0, σ_1, σ_2, σ_3)
     end

end

#=
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


@doc raw"""
     σ_0, σ_1, σ_2, σ_3

These are the results of the following integral:
```math
     \sigma_i = \int_0^\infty \frac{\mathrm{d} q}{2 \pi^2} q^{2-i} \, P(q) 
```
where  `P(q)` is the input Power Spectrum.
"""
σ_0, σ_1, σ_2, σ_3
=#


##########################################################################################92



@doc raw"""
     ϕ(s; s_min = s_MIN, s_max = s_MAX) :: Float64

Radial part of the survey window function. Return `1.0` if is true that
``s_\mathrm{min} \le s \le s_\mathrm{max}`` and `0.0` otherwise.

In this software we made the assuption that the survey window function can be
separated into a radial and angular part, i.e.:

```math
     \phi(\mathbf{s}) = \phi(s) \, W(\hat{s})
```

See also: [`W`](@ref)
"""
ϕ(s; s_min = s_MIN, s_max = s_MAX) = s_min < s < s_max ? 1.0 : 0.0



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
     V_survey(s_min = s_MIN, s_max = s_MAX, θ_max = θ_MAX) :: Float64

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
function V_survey(s_min = s_MIN, s_max = s_MAX, θ_max = θ_MAX)
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
     A(s_min = s_MIN, s_max = s_MAX, θ_max = θ_MAX) :: Float64

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
function A(s_min = s_MIN, s_max = s_MAX, θ_max = θ_MAX)
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

