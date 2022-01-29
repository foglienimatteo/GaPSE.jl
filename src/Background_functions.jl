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

function z_eff(s_min = s_min, s_max = s_max, θ_max = θ_MAX)
     int_w2 = 2 * π * (1 - cos(θ_max))
     int_z_ϕ2 = quadgk(s -> s^2 * z_of_s(s), s_min, s_max)[1]
     int_ϕ2 = quadgk(s -> s^2, s_min, s_max)[1]

     return int_z_ϕ2 / int_ϕ2
end

const f0 = data[:, column_NAMES_BACKGROUND["gr.fac. f"]][end]
const D0 = data[:, column_NAMES_BACKGROUND["gr.fac. D"]][end]
const ℋ0 = data[:, column_NAMES_BACKGROUND["H [1/Mpc]"]][end] / h_0
const ℋ0_p = 0.0
const s_b0 = 0.0
const s_min = s_of_z(z_MIN)
const s_max = s_of_z(z_MAX)
const s_eff = s_of_z(z_eff())



##########################################################################################92


s(s1, s2, y) = √(s1^2 + s2^2 - 2 * s1 * s2 * y)
s2(s1, s, μ) = √(s1^2 + s^2 + 2 * s1 * s * μ)
y(s1, s, μ) = (μ * s + s1) / s2(s, s1, μ)




