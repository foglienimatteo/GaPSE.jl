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


data = readdlm(FILE_BACKGROUND, comments = true)
N_z_MAX = findfirst(z -> z <= z_MAX, data[:, 1]) - 1
N_z_MIN = findfirst(z -> z <= z_MIN, data[:, 1]) + 1
#println("\nN_z_MAX = ", N_z_MAX, ", \t\t z[N_z_MAX] = ", data[:, 1][N_z_MAX])
#println("N_z_MIN = ", N_z_MIN, ", \t\t z[N_z_MIN] = ", data[:, 1][N_z_MIN])
data_dict = Dict([name => reverse(data[:, i][N_z_MAX:N_z_MIN]) for (i, name) in enumerate(NAMES_BACKGROUND)]...)


comdist_array = data_dict["comov. dist."] .* h_0
#angdist_array = angdist_array .* h_0
#lumdist_array = lumdist_array .* h_0
D = Spline1D(comdist_array, data_dict["gr.fac. D"])
f = Spline1D(comdist_array, data_dict["gr.fac. f"])
ℋ = Spline1D(comdist_array, data_dict["H [1/Mpc]"] ./ h_0 ./ (1.0 .+ data_dict["z"]))
ℋ_p(s) = Dierckx.derivative(ℋ, s)
#s_b = Spline1D(comdist_array, [0.0 for i in 1:length(data_dict["comov. dist."])])
s_b(s) = 0.0
s_of_z = Spline1D(data_dict["z"], comdist_array)
z_of_s = Spline1D(comdist_array, data_dict["z"])
f_evo = 0

ℛ(s) = 5 * s_b(s) + (2 - 5 * s_b(s)) / (ℋ(s) * s) + ℋ_p(s) / (ℋ(s)^2) - f_evo
function ℛ(s, ℋ, ℋ_p, s_b, f_evo = f_evo)
     5 * s_b + (2 - 5 * s_b) / (ℋ * s) + ℋ_p / (ℋ^2) - f_evo
end



@doc raw"""
     const z_eff :: Float64

The effective redshift ``z_\mathrm{eff}`` is calcuated as follows:
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
const z_eff = 3.0 / (s_max^3 - s_min^3) * quadgk(s -> s^2 * z_of_s(s), s_min, s_max)[1]


const f0 = data[:, column_NAMES_BACKGROUND["gr.fac. f"]][end]
const D0 = data[:, column_NAMES_BACKGROUND["gr.fac. D"]][end]
const ℋ0 = data[:, column_NAMES_BACKGROUND["H [1/Mpc]"]][end] / h_0
const ℋ0_p = 0.0
const s_b0 = 0.0
const s_min = s_of_z(z_MIN)
const s_max = s_of_z(z_MAX)
const s_eff = s_of_z(z_eff)



@doc raw"""
     const ℋ0 :: Float64

Comoving Hubble constant at present time. Its value is, in natural system
(where the speed of light c=1): 
``\mathcal{H}_0 \simeq 3.335641\times10^{-4} \; h_0^{-1}\mathrm{Mpc}``
"""
ℋ0

##########################################################################################92



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




