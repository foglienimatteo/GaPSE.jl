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


#FILE_F_MAP = "data/F_map_stable_2.txt"
#NAMES_F_MAP = ["x", "mu", "F", "F_error"]
FILE_F_MAP = "/Users/matteofoglieni/AAA_TESI_MAGISTRALE/GaPSE-free-ipynb/PANTIRI_F_x_mu.txt"
NAMES_F_MAP = ["x", "mu", "F"]

F_map_data = readdlm(FILE_F_MAP, comments = true)
F_map_data_dict = Dict([name => F_map_data[2:end, i] for (i, name) in enumerate(NAMES_F_MAP)]...)

_xs = unique(F_map_data_dict["x"])
_μs = unique(F_map_data_dict["mu"])
#_Fs = reshape(F_map_data_dict["F"], (length(_μs), length(_xs)))' # FOR SciPy DOES NOT WORK
_Fs = F_map_data_dict["F"]


# for my F map with GridInterpolations
#my_F_grid = GridInterpolations.RectangleGrid(_μs, _xs)
#spline_F(x, μ) = GridInterpolations.interpolate(my_F_grid, _Fs, [μ, x])

# for mattia F map with GridInterpolations
mattia_F_grid = GridInterpolations.RectangleGrid(_xs, _μs)
spline_F(x, μ) = GridInterpolations.interpolate(mattia_F_grid, _Fs, [x, μ])

# for my F with RectBivariateSpline
#my_scipy_grid_Fs = reshape(_Fs, (length(_μs), length(_xs)))
#not_my_spline = SciPy.interpolate.RectBivariateSpline(_μs, _xs, my_scipy_grid_Fs)
#spline_F(x, μ) = not_my_spline(μ, x)[1]

# for mattia F with RectBivariateSpline
#mattia_scipy_grid_Fs = reshape(_Fs, (length(_xs), length(_μs)))
#not_mattia_spline = SciPy.interpolate.RectBivariateSpline(_xs, _μs, mattia_scipy_grid_Fs)
#spline_F(x, μ) = not_mattia_spline(x, μ)[1]


##########################################################################################92


data = readdlm(FILE_BACKGROUND, comments = true)
N_z_MAX = findfirst(z -> z <= z_MAX, data[:, 1]) - 1
N_z_MIN = findfirst(z -> z <= z_MIN, data[:, 1]) + 1
#println("\nN_z_MAX = ", N_z_MAX, ", \t\t z[N_z_MAX] = ", data[:, 1][N_z_MAX])
#println("N_z_MIN = ", N_z_MIN, ", \t\t z[N_z_MIN] = ", data[:, 1][N_z_MIN])
data_dict = Dict([name => reverse(data[:, i][N_z_MAX:N_z_MIN]) for (i, name) in enumerate(NAMES_BACKGROUND)]...)


comdist_array = data_dict["comov. dist."] .* h_0
#angdist_array = angdist_array*h_0
#lumdist_array = lumdist_array*h_0
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

const f0 = data[:, column_NAMES_BACKGROUND["gr.fac. f"]][end]
const D0 = data[:, column_NAMES_BACKGROUND["gr.fac. D"]][end]
const ℋ0 = data[:, column_NAMES_BACKGROUND["H [1/Mpc]"]][end] / h_0
const ℋ0_p = 0.0
const s_b0 = 0.0


##########################################################################################92


ps = readdlm(FILE_PS, comments = true)
ps_dict = Dict([name => ps[:, i] for (i, name) in enumerate(NAMES_PS)]...)

PK = Spline1D(ps_dict["k (h/Mpc)"], ps_dict["P (Mpc/h)^3"])

N = 1024                                # number of points to use in the Fourier transform
k_max = ps_dict["k (h/Mpc)"][end]        # maximum k-value
k_min = ps_dict["k (h/Mpc)"][begin]      # minimum k-value
s0 = 1 / k_max;                          # minimum r-value (should be ~1/k_max)

I00 = Spline1D(xicalc(PK, 0, 0; N = N, kmin = k_min, kmax = k_max, r0 = s0)...)
I20 = Spline1D(xicalc(PK, 2, 0; N = N, kmin = k_min, kmax = k_max, r0 = s0)...)
I40 = Spline1D(xicalc(PK, 4, 0; N = N, kmin = k_min, kmax = k_max, r0 = s0)...)
I02 = Spline1D(xicalc(PK, 0, 2; N = N, kmin = k_min, kmax = k_max, r0 = s0)...)
I22 = Spline1D(xicalc(PK, 2, 2; N = N, kmin = k_min, kmax = k_max, r0 = s0)...)
I31 = Spline1D(xicalc(PK, 3, 1; N = N, kmin = k_min, kmax = k_max, r0 = s0)...)
I13 = Spline1D(xicalc(PK, 1, 3; N = N, kmin = k_min, kmax = k_max, r0 = s0)...)
I11 = Spline1D(xicalc(PK, 1, 1; N = N, kmin = k_min, kmax = k_max, r0 = s0)...)

σ_0 = quadgk(q -> PK(q) * q^2 / (2 * π^2), k_min, k_max)[1]
σ_1 = quadgk(q -> PK(q) * q / (2 * π^2), k_min, k_max)[1]
σ_2 = quadgk(q -> PK(q) / (2 * π^2), k_min, k_max)[1]
σ_3 = quadgk(q -> PK(q) / (2 * π^2 * q), k_min, k_max)[1]


##########################################################################################92


s(s1, s2, y) = √(s1^2 + s2^2 - 2 * s1 * s2 * y)
s2(s1, s, μ) = √(s1^2 + s^2 + 2 * s1 * s * μ)
y(s1, s, μ) = (μ * s + s1) / s2(s, s1, μ)


const s_min = s_of_z(z_MIN)
const s_max = s_of_z(z_MAX)
ϕ(s; s_min = s_min, s_max = s_max) = s_min < s < s_max ? 1.0 : 0.0
W(θ; θ_max = θ_MAX) = 0 < θ < θ_max ? 1.0 : 0.0
function V(s_min = s_min, s_max = s_max, θ_max = θ_MAX)
     r1, r2 = s_min * sin(θ_max), s_max * sin(θ_max)
     d1, d2 = s_min * cos(θ_max), s_max * cos(θ_max)
     calotta_up = π / 3 * (s_max - d2)^2 * (2 * s_max + d2)
     calotta_down = π / 3 * (s_min - d1)^2 * (2 * s_min + d1)
     tronco_cono = π / 3 * (r1^2 + r1 * r2 + r2^2) * (s_max - s_min) * cos(θ_max)
     return calotta_up + tronco_cono - calotta_down
end

function A(s_min = s_min, s_max = s_max, θ_max = θ_MAX)
     2.0 * π * V(s_min, s_max, θ_max)
end

function z_eff(s_min = s_min, s_max = s_max, θ_max = θ_MAX)
     int_w2 = 2 * π * (1 - cos(θ_max))
     int_z_ϕ2 = quadgk(s -> s^2 * z_of_s(s), s_min, s_max)[1]
     int_ϕ2 = quadgk(s -> s^2, s_min, s_max)[1]

     return int_z_ϕ2 / int_ϕ2
end

const s_eff = s_of_z(z_eff())


##########################################################################################92


function parameters_used(io::IO)
     println(io, "# The following parameters were used for this computation: ")
     println(io, "# CLASS Power Spectrum input file : \"$(FILE_PS)\"")
     println(io, "# k_min = $k_min \t k_max = $k_max")
     println(io, "# F window function input file : \"$(FILE_F_MAP)\"")
     println(io, "# CLASS Background input file: \"$(FILE_BACKGROUND)\"")
     println(io, "# \t h_0 = $h_0 \t \t EVERYTHING IS MEASURED WITHOUT h_0!")
     println(io, "# \t comoving H_0 = $(@sprintf("%.6e", ℋ0)) h_0/Mpc")
     println(io, "# \t growth factor D_0 = $D0")
     println(io, "# \t growth rate f_0 = $(@sprintf("%.6f", f0))")
     println(io, "# \t z_min = $z_MIN \t\tcomoving s_min = " *
                 "$(@sprintf("%.5f", s_min)) Mpc/h_0")
     println(io, "# \t z_max = $z_MAX \t\tcomoving s_max = " *
                 "$(@sprintf("%.5f", s_max)) Mpc/h_0")
     println(io, "# \t z_eff = $(@sprintf("%.5f", z_eff())) \tcomoving s_eff = " *
                 "$(@sprintf("%.5f", s_eff)) Mpc/h_0")
     println(io, "# \t Ω_b = $Ω_b \t Ω_cdm = $Ω_cdm \t Ω_M0 = $Ω_M0")
     println(io, "# \t Volume of the survey V = $(@sprintf("%.6e", V()))")
     println(io, "# ")
end



##########################################################################################92


parameters_used(stdout)
