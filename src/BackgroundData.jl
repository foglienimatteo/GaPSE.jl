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
     const f0 :: Float64

Linear growth rate at present time. Its value is equal to:
```math
     f_0 \simeq 0.5126998572951
```
"""
const f0 = 5.126998572951e-01


@doc raw"""
     const D0 :: Float64

Linear growth factor at present time. Its value is equal to:
```math
     D_0 = 1.0
```
"""
const D0 = 1.0


@doc raw"""
     const ℋ0 :: Float64

Comoving Hubble constant at present time. Its value is, in natural system
(where the speed of light c=1): 
``\mathcal{H}_0 \simeq 3.335641\times10^{-4} \; h_0^{-1}\mathrm{Mpc}``
"""
const ℋ0 = 3.3356409519815204e-4 # h_0/Mpc



##########################################################################################92



struct BackgroundData
     z::Vector{Float64}
     conftime::Vector{Float64}
     comdist::Vector{Float64}
     angdist::Vector{Float64}
     lumdist::Vector{Float64}
     D::Vector{Float64}
     f::Vector{Float64}
     ℋ::Vector{Float64}
     ℋ_p::Vector{Float64}

     function BackgroundData(file::String, z_max;
          names = NAMES_BACKGROUND, h = 0.7)
     
          I_redshift = findfirst(x -> x == "z", names)
          I_comdist = findfirst(x -> x == "comov. dist.", names)
     
          data = readdlm(file, comments = true)
     
          N_z_MAX = findfirst(z -> z <= z_max, data[:, I_redshift]) - 1
          com_dist_z_MAX = data[:, I_comdist][N_z_MAX]
          N_2_com_dist_z_MAX = findfirst(s -> s <= 3.0 * com_dist_z_MAX, data[:, I_comdist]) - 1
     
          data_dict = Dict([name => reverse(data[:, i][N_2_com_dist_z_MAX:end])
                            for (i, name) in enumerate(names)]...)
     
          com_H = data_dict["H [1/Mpc]"] ./ h ./ (1.0 .+ data_dict["z"])
          conf_time = data_dict["conf. time [Mpc]"] .* h
          spline_com_H = Spline1D(reverse(conf_time), reverse(com_H); bc = "nearest")
          com_H_p = [Dierckx.derivative(spline_com_H, t) for t in conf_time]
     
          new(
               data_dict["z"],
               conf_time,
               data_dict["comov. dist."] .* h,
               data_dict["ang.diam.dist."] .* h,
               data_dict["lum. dist."] .* h,
               data_dict["gr.fac. D"],
               data_dict["gr.fac. f"],
               com_H,
               com_H_p,
          )
     end
end



struct CosmoParams
     z_min::Float64
     z_max::Float64
     θ_max::Float64

     k_min::Float64
     k_max::Float64

     Ω_b::Float64
     Ω_cdm::Float64
     Ω_M0::Float64
     h_0::Float64

     N::Int
     fit_min::Union{Float64,Nothing}
     fit_max::Union{Float64,Nothing}
     con::Bool
     s_lim::Union{Float64, Nothing}


     function CosmoParams(z_min, z_max, θ_max;
          k_min = 1e-8, k_max = 10.0,
          Ω_b = 0.0489, Ω_cdm = 0.251020, h_0 = 0.70,
          N = 1024, fit_min = 2.0, fit_max = 10.0, con = true, s_lim = nothing)
     
          @assert 0.0 ≤ θ_max ≤ π / 2.0 " 0.0 ≤ θ_max ≤ π/2.0 must hold!"
          @assert 0.0 ≤ z_min < z_max " 0.0 ≤ z_min < z_max must hold!"
          @assert 0.0 ≤ k_min < k_max " 0.0 ≤ k_min < k_max must hold!"
          @assert Ω_b ≥ 0.0 " Ω_b ≥ 0.0 must hold!"
          @assert Ω_cdm ≥ 0.0 " Ω_cdm ≥ 0.0 must hold!"
          @assert 0.0 ≤ h_0 ≤ 1.0 " 0.0 ≤ h_0 ≤ 1.0 must hold!"
          @assert isnothing(s_lim) || (0.0 < s_lim < 50.0) "0.0 < s_lim < 50.0 must hold!"
     
          new(z_min, z_max, θ_max, k_min, k_max, Ω_b, Ω_cdm, Ω_cdm + Ω_b, h_0,
               N, fit_min, fit_max, con, s_lim)
     end
end


##########################################################################################92



@doc raw"""
     func_z_eff(s_min, s_max, z_of_s) :: Float64

Return the effective redshift ``z_\mathrm{eff}``, calcuated as follows:
```math
\begin{align*}
z_\mathrm{eff} := 
    \frac{
        \int \mathrm{d}^3\mathbf{s} \, \phi^2(\mathbf{s}) \, z(s)
     }{
         \int \mathrm{d}^3\mathbf{s}\, \phi^2(\mathbf{s}) 
      } &= \frac{
          \int_0^\infty \mathrm{d}s  \, s^2 \, \phi^2(s) \, z(s) \times
          \int_{4\pi}\mathrm{d}^2\hat{\mathbf{s}} \, W^2(\hat{\mathbf{s}})
      }{
          \int_0^\infty \mathrm{d}s \, s^2 \, \phi^2(s)\times
          \int_{4\pi}\mathrm{d}^2\hat{\mathbf{s}} \, W^2(\hat{\mathbf{s}})
      } \\[5pt]
      &= \frac{
          \int_0^\infty \mathrm{d}s  \, s^2 \, \phi^2(s) \, z(s)
      }{
          \int_0^\infty \mathrm{d}s \, s^2 \, \phi^2(s)
      } \\[4pt]
      &= \frac{3}{s_\mathrm{max}^3 - s_\mathrm{min}^3} \,
          \int_{s_\mathrm{min}}^{s_\mathrm{max}} \mathrm{d}s  \, s^2 \, z(s)
\end{align*}
```
where we have used our assuption on separability of the window function
```math
     \phi(\mathbf{s}) = \phi(s) \, W(\hat{s})
```
and their definitions.


See also: [`ϕ`](@ref), [`W`](@ref)
"""
function func_z_eff(s_min, s_max, z_of_s)
     3.0 / (s_max^3 - s_min^3) * quadgk(s -> s^2 * z_of_s(s), s_min, s_max)[1]
end


@doc raw"""
     s(s1, s2, y) :: Float64

Return the value ``s = \sqrt{s_1^2 + s_2^2 - 2 \, s_1 \, s_2 \, y}``

See also: [`μ`](@ref), [`s2`](@ref), [`y`](@ref)
"""
s(s1, s2, y) = √(s1^2 + s2^2 - 2 * s1 * s2 * y)


@doc raw"""
     μ(s1, s2, y) :: Float64

Return the value ``\mu=\hat{\mathbf{s}}_1\dot\hat{\mathbf{s}}``, defined as:
```math
\mu = \mu(s_1, s_2, y) = \frac{y \, s_2 - s_1}{s(s_1, s_2, y)} \;,
\quad s(s_1, s_2, y) = \sqrt{s_1^2 + s^2 - 2 \, s_1 \, s_2 \, y}
```
with ``y=\cos\theta=\hat{\mathbf{s}}_1\dot\hat{\mathbf{s}}`` and where ``s`` is 
obtained from the function `s`

See also: [`s`](@ref), [`s2`](@ref), [`y`](@ref)
"""
μ(s1, s2, y) = (y * s2 - s1) / s(s1, s2, y)



@doc raw"""
     s2(s1, s, μ) :: Float64

Return the value ``s_2 = \sqrt{s_1^2 + s^2 + 2 \, s_1 \, s \, \mu}``

See also: [`s`](@ref), [`μ`](@ref), [`y`](@ref)
"""
s2(s1, s, μ) = √(s1^2 + s^2 + 2 * s1 * s * μ)



@doc raw"""
     y(s1, s, μ) :: Float64

Return the value ``y=\cos\theta``, defined as:
```math
y = y(s_1, s, \mu) = \frac{\mu \, s + s_1}{s2(s_1, s, \mu)} \;,
\quad s_2 = \sqrt{s_1^2 + s^2 + 2 \, s_1 \, s \, \mu}
```
with ``\mu=\hat{\mathbf{s}}_1\dot\hat{\mathbf{s}}_2`` and 
where ``s_2`` is btained from the function `s2`

See also: [`s`](@ref), [`μ`](@ref), [`s2`](@ref)
"""
y(s1, s, μ) = (μ * s + s1) / s2(s, s1, μ)




