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
     GridInterpolations.interpolate(grid, reshape(str.Fs, (:, 1)), [x, μ])
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



@doc raw"""
     InputPS(ks::Vector{Float64}, pks::Vector{Float64})

Store the Input Power Spectrum.

## Constructors

- `InputPS(file::String; expand::Bool = true)` : read the IPS from
  the given input `file`; it can contain comments (defined with a 
  starting # on each line), but the file structure mus be space-separated.

- `InputPS(ks::AbstractVector{T1}, pks::AbstractVector{T2}; expand::Bool = true) where {T1,T2}` :
  take `ks` as x-axis and `pks` as y-axis values of the IPS

## Optional arguments

- `expand::Bool = true` : if `true`, the IPS is expanded toward small (until 
  ``k \geq 10^{-8}  \; h_0 \, \mathrm{Mpc}^{-1}``) and high (until 
  ``k \seq 10^{3}  \; h_0 \, \mathrm{Mpc}^{-1}``) k-values; if the IPS already 
  arrives at this limits, nothing is changed.
  The expansion is made through `expanded_IPS`
  
See also: [`expanded_IPS`](@ref)
"""
struct InputPS
     ks::Vector{Float64}
     pks::Vector{Float64}

     function InputPS(file::String; expand::Bool = true)
          data = readdlm(file, comments = true)
          @assert size(data[:, 1]) == size(data[:, 2]) "ks and pks must have the same length!"

          ks, pks = expand ? begin
               println("I expand the input power spectrum at its extremes.")
               expanded_IPS(data[:, 1], data[:, 2])
          end : begin
               println("I take the input power spectrum as it is,without expanding.")
               (data[:, 1], data[:, 2])
          end

          new(ks, pks)
     end

     function InputPS(ks::AbstractVector{T1}, pks::AbstractVector{T2};
          expand::Bool = true) where {T1,T2}

          @assert size(ks) == size(pks) "ks and pks must have the same length!"
          new_ks, new_pks = expand ? begin
               println("I expand the input power spectrum at its extremes.")
               expanded_IPS(ks, pks; con = false)
          end : begin
               println("I take the input power spectrum as it is,without expanding.")
               (ks, pks)
          end

          new(new_ks, new_pks)
     end
end


##########################################################################################92




struct IntegralIPS
     si::Float64
     b::Float64
     a::Float64
     left::Float64
     spline::Dierckx.Spline1D

     function IntegralIPS(ips, l, n; N = 1024, kmin = 1e-4, kmax = 1e3, s0 = 1e-3,
          fit_min = 2.0, fit_max = 10.0, p0 = nothing, con = false)

          rs, xis = xicalc(ips, l, n; N = N, kmin = kmin, kmax = kmax, r0 = s0)

          p_0 = isnothing(p0) ? (con == true ? [-1.0, 1.0, 0.0] : [-1.0, 1.0]) : p0
          si, b, a = power_law_from_data(
               rs, xis, p_0, fit_min, fit_max; con = con)

          ind = findfirst(x->x>fit_min, rs)-1
          new_rs = vcat(rs[ind:end])
          new_Is = vcat(xis[ind:end])
          spline = Spline1D(new_rs, new_Is; bc="error")

          new(si, b, a, fit_min, spline)
     end

     function IntegralIPS(ips, func::Function; N = 1024, kmin = 1e-4, kmax = 1e3,
          fit_min = 0.1, fit_max = 1.0, p0 = nothing, con = false, kwargs...)

          ss = 10 .^ range(log10(0.9*fit_min), 4, length = 1024)
          Is = [func(ips, s, kmin, kmax; kwargs...) for s in ss]

          p_0 = isnothing(p0) ? (con == true ? [-2.0, -1.0, 0.0] : [-2.0, -1.0]) : p0
          si, b, a = power_law_from_data(
               ss, Is, p_0, fit_min, fit_max; con = con)
          #println("si, b, a = $si , $b , $a")

          spline = Spline1D(ss, Is; bc="error")

          new(si, b, a, fit_min, spline)
     end

end

function (Iln::IntegralIPS)(x)
     if x < Iln.left
          return power_law(x, Iln.si, Iln.b, Iln.a)
     else 
          return Iln.spline(x)
     end
end


@doc raw"""
     IPSTools(
          I00::Dierckx.Spline1D
          I20::Dierckx.Spline1D
          I40::Dierckx.Spline1D
          I02::Dierckx.Spline1D
          I22::Dierckx.Spline1D
          I31::Dierckx.Spline1D
          I13::Dierckx.Spline1D
          I11::Dierckx.Spline1D

          I04_tilde::Dierckx.Spline1D

          σ_0::Float64
          σ_1::Float64
          σ_2::Float64
          σ_3::Float64

          fit_min::Union{Float64,Nothing}
          fit_max::Union{Float64,Nothing}
          k_min::Float64
          k_max::Float64
          s_0::Float64
          )

Struct that contains all the function and values obtained from the 
Input Power Spectrum.

## Arguments

- `I00, I20, I40, I02, I22, I31, I13, I11 ::Dierckx.Spline1D`: spline that
  return the value of the integral:

  ```math
  I_\ell^n(s) = \int_0^\infty \frac{\mathrm{d} q}{2 \pi^2} q^2 \, P(q) 
    \, \frac{j_\ell(qs)}{(qs)^n}
  ```
  where, for a generic `Iab` name, ``\ell`` is the FIRST number (`a`) and 
  ``n`` the second (`b`).

  These function are obtained through a `Spline1D` (from the 
  [Dierckx](https://github.com/kbarbary/Dierckx.jl) Julia package) of the Spherical
  Bessel Transform function `xicalc` (from the 
  [TwoFAST](https://github.com/hsgg/TwoFAST.jl) Julia package) applied to the 
  input Power Spectrum `P(q)`. 

- `I04_tilde::Dierckx.Spline1D`: spline that
  return the value of the integral:

  ```math
  I_\ell^n(s) = \int_0^\infty \frac{\mathrm{d} q}{2 \pi^2} q^2 \, P(q) 
    \, \frac{j_\ell(qs)}{(qs)^n}
  ```

"""
struct IPSTools
     #=
     I00::Dierckx.Spline1D
     I20::Dierckx.Spline1D
     I40::Dierckx.Spline1D
     I02::Dierckx.Spline1D
     I22::Dierckx.Spline1D
     I31::Dierckx.Spline1D
     I13::Dierckx.Spline1D
     I11::Dierckx.Spline1D

     I04_tilde::Dierckx.Spline1D
     =#

     I00::IntegralIPS
     I20::IntegralIPS
     I40::IntegralIPS
     I02::IntegralIPS
     I22::IntegralIPS
     I31::IntegralIPS
     I13::IntegralIPS
     I11::IntegralIPS

     I04_tilde::IntegralIPS

     σ_0::Float64
     σ_1::Float64
     σ_2::Float64
     σ_3::Float64

     fit_min::Union{Float64,Nothing}
     fit_max::Union{Float64,Nothing}
     k_min::Float64
     k_max::Float64
     s_0::Float64

     function IPSTools(
          ips::InputPS;
          N = 1024,
          fit_min::Float64 = 0.05,
          fit_max::Float64 = 0.5,
          con::Bool = false,
          k_min::Float64 = 1e-6,
          k_max::Float64 = 10.0,
          lim::Float64 = 1e-8
     )

          PK = Spline1D(ips.ks, ips.pks; bc = "error")

          #kmin, kmax = min(ips.ks...), max(ips.ks...)
          kmin, kmax, s0 = 1e-5, 1e3, 1e-3

          p0 = con ? [-1.0, 1.0, 0.0] : [-1.0, 1.0]

          I00 = IntegralIPS(PK, 0, 0; N = N, kmin = kmin, kmax = kmax, s0 = s0,
                    fit_min = fit_min, fit_max = fit_max, p0 = p0, con = con)
          I20 = IntegralIPS(PK, 2, 0; N = N, kmin = kmin, kmax = kmax, s0 = s0,
                    fit_min = fit_min, fit_max = fit_max, p0 = p0, con = con)
          I40 = IntegralIPS(PK, 4, 0; N = N, kmin = kmin, kmax = kmax, s0 = s0,
                    fit_min = fit_min, fit_max = fit_max, p0 = p0, con = con)
          I02 = IntegralIPS(PK, 0, 2; N = N, kmin = kmin, kmax = kmax, s0 = s0,
                    fit_min = fit_min, fit_max = fit_max, p0 = p0, con = con)
          I22 = IntegralIPS(PK, 2, 2; N = N, kmin = kmin, kmax = kmax, s0 = s0,
                    fit_min = fit_min, fit_max = fit_max, p0 = p0, con = con)
          I31 = IntegralIPS(PK, 3, 1; N = N, kmin = kmin, kmax = kmax, s0 = s0,
                    fit_min = fit_min, fit_max = fit_max, p0 = p0, con = con)
          I13 = IntegralIPS(PK, 1, 3; N = N, kmin = kmin, kmax = kmax, s0 = s0,
                    fit_min = fit_min, fit_max = fit_max, p0 = p0, con = con)
          I11 = IntegralIPS(PK, 1, 1; N = N, kmin = kmin, kmax = kmax, s0 = s0,
                    fit_min = fit_min, fit_max = fit_max, p0 = p0, con = con)


          I04_tilde = IntegralIPS(PK, func_I04_tilde; N = N, kmin = kmin, kmax = kmax,
               fit_min = 0.1, fit_max = 1.0, p0 = nothing, con = con)


          #=
          I00 = Spline1D(expanded_Iln(PK, 0, 0; lim = lim, N = N, kmin = kmin, kmax = kmax, s0 = s0,
                    fit_min = fit_min, fit_max = fit_max, p0 = p0, con = con)...; bc = "error")
          I20 = Spline1D(expanded_Iln(PK, 2, 0; lim = lim, N = N, kmin = kmin, kmax = kmax, s0 = s0,
                    fit_min = fit_min, fit_max = fit_max, p0 = p0, con = con)...; bc = "error")
          I40 = Spline1D(expanded_Iln(PK, 4, 0; lim = lim, N = N, kmin = kmin, kmax = kmax, s0 = s0,
                    fit_min = fit_min, fit_max = fit_max, p0 = p0, con = con)...; bc = "error")
          I02 = Spline1D(expanded_Iln(PK, 0, 2; lim = lim, N = N, kmin = kmin, kmax = kmax, s0 = s0,
                    fit_min = fit_min, fit_max = fit_max, p0 = p0, con = con)...; bc = "error")
          I22 = Spline1D(expanded_Iln(PK, 2, 2; lim = lim, N = N, kmin = kmin, kmax = kmax, s0 = s0,
                    fit_min = fit_min, fit_max = fit_max, p0 = p0, con = con)...; bc = "error")
          I31 = Spline1D(expanded_Iln(PK, 3, 1; lim = lim, N = N, kmin = kmin, kmax = kmax, s0 = s0,
                    fit_min = fit_min, fit_max = fit_max, p0 = p0, con = con)...; bc = "error")
          I13 = Spline1D(expanded_Iln(PK, 1, 3; lim = lim, N = N, kmin = kmin, kmax = kmax, s0 = s0,
                    fit_min = fit_min, fit_max = fit_max, p0 = p0, con = con)...; bc = "error")
          I11 = Spline1D(expanded_Iln(PK, 1, 1; lim = lim, N = N, kmin = kmin, kmax = kmax, s0 = s0,
                    fit_min = fit_min, fit_max = fit_max, p0 = p0, con = con)...; bc = "error")

          #ss = 10 .^ range(log10(s0), log10(s0) - log10(kmin) - log10(kmax), length = N)
          ss = 10 .^ range(log10(lim), 6, length = 1024)
          I04_tildes = expanded_I04_tilde(PK, ss; kmin = kmin, kmax = kmax)
          #I04_tildes = [func_I04_tilde(PK, s, kmin, kmax) for s in ss]
          I04_tilde = Spline1D(ss, I04_tildes; bc = "error")
          =#

          σ_0 = quadgk(q -> PK(q) * q^2 / (2 * π^2), k_min, k_max)[1]
          σ_1 = quadgk(q -> PK(q) * q / (2 * π^2), k_min, k_max)[1]
          σ_2 = quadgk(q -> PK(q) / (2 * π^2), k_min, k_max)[1]
          σ_3 = quadgk(q -> PK(q) / (2 * π^2 * q), k_min, k_max)[1]

          new(I00, I20, I40, I02, I22, I31, I13, I11, I04_tilde, σ_0, σ_1, σ_2, σ_3,
               fit_min, fit_max, k_min, k_max, s0)
     end

     function IPSTools(ips::InputPS, iIs::String;
          k_min::Float64 = 1e-8,
          k_max::Float64 = 10.0
     )
          PK = Spline1D(ips.ks, ips.pks; bc = "error")

          tab_Is = readdlm(iIs, comments = true)
          ss = convert(Vector{Float64}, tab_Is[2:end, 1])

          #kmin, kmax = min(ips.ks...), max(ips.ks...)
          kmin, kmax, s0 = 1e-5, 1e3, 1e-3

          I00 = Spline1D(ss, convert(Vector{Float64}, tab_Is[2:end, 2]); bc = "error")
          I20 = Spline1D(ss, convert(Vector{Float64}, tab_Is[2:end, 3]); bc = "error")
          I40 = Spline1D(ss, convert(Vector{Float64}, tab_Is[2:end, 4]); bc = "error")
          I02 = Spline1D(ss, convert(Vector{Float64}, tab_Is[2:end, 5]) ./ ss .^ 2; bc = "error")
          I22 = Spline1D(ss, convert(Vector{Float64}, tab_Is[2:end, 6]) ./ ss .^ 2; bc = "error")
          I31 = Spline1D(ss, convert(Vector{Float64}, tab_Is[2:end, 7]) ./ ss; bc = "error")
          I11 = Spline1D(ss, convert(Vector{Float64}, tab_Is[2:end, 8]) ./ ss; bc = "error")
          I13 = Spline1D(xicalc(PK, 1, 3; N = 1024, kmin = kmin, kmax = kmax, r0 = s0)...; bc = "error")

          #ss = 10 .^ range(log10(s0), log10(s0) - log10(kmin) - log10(kmax), length = N)
          ss = 10 .^ range(log10(lim), 4, length = 1000)
          I04_tildes = expanded_I04_tilde(PK, ss; kmin = kmin, kmax = kmax)
          I04_tilde = Spline1D(ss, I04_tildes; bc = "error")

          σ_0 = quadgk(q -> PK(q) * q^2 / (2 * π^2), k_min, k_max)[1]
          σ_1 = quadgk(q -> PK(q) * q / (2 * π^2), k_min, k_max)[1]
          σ_2 = quadgk(q -> PK(q) / (2 * π^2), k_min, k_max)[1]
          σ_3 = quadgk(q -> PK(q) / (2 * π^2 * q), k_min, k_max)[1]

          new(I00, I20, I40, I02, I22, I31, I13, I11, I04_tilde, σ_0, σ_1, σ_2, σ_3,
               nothing, nothing, k_min, k_max, s0)
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
     ϕ(s, s_min, s_max) :: Float64

Radial part of the survey window function. Return `1.0` if is true that
``s_\mathrm{min} \le s \le s_\mathrm{max}`` and `0.0` otherwise.

In this software we made the assuption that the survey window function can be
separated into a radial and angular part, i.e.:

```math
     \phi(\mathbf{s}) = \phi(s) \, W(\hat{s})
```

See also: [`W`](@ref)
"""
ϕ(s, s_min, s_max) = s_min < s < s_max ? 1.0 : 0.0



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
W(θ, θ_max) = 0.0 ≤ θ < θ_max ? 1.0 : 0.0


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
function V_survey(s_min, s_max, θ_max)
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
     A(s_min, s_max, θ_max) :: Float64

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
function A(s_min, s_max, θ_max)
     2.0 * π * V_survey(s_min, s_max, θ_max)
end



@doc raw"""
     A_prime :: Float64

It's the Power Spectrum multipole normalization coefficient ``A^{'}``, i.e.:
```math
     A^{'} = \frac{3 \, A}{ (s_\mathrm{max}^3 - s_\mathrm{min}^3)} = 8 \pi^2
```

See also: [`A`](@ref), [`V_survey`](@ref)
"""
const A_prime = 8.0 * π^2

