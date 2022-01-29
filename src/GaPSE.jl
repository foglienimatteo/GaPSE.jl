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



module GaPSE

using TwoFAST # Licence: MIT "Expat" (o GPL ?)
using Dierckx # Licence: BSD
using HCubature, QuadGK, LegendrePolynomials  # Licence: MIT "Expat"
using GridInterpolations  # Licence: MIT "Expat"
using ProgressMeter, Documenter  # Licence: MIT "Expat"
using Test, Printf, DelimitedFiles  # Licence: MIT "Expat"

const Ω_b = 0.0489
const Ω_cdm = 0.251020
const Ω_M0 = Ω_b + Ω_cdm

const h_0 = 0.70
const z_MIN = 0.05
const z_MAX = 0.2
const θ_MAX = π / 2

#FILE_F_MAP = "data/F_map_stable_2.txt"
#NAMES_F_MAP = ["x", "mu", "F", "F_error"]
FILE_F_MAP = "/Users/matteofoglieni/AAA_TESI_MAGISTRALE/GaPSE-free-ipynb/PANTIRI_F_x_mu.txt"
NAMES_F_MAP = ["x", "mu", "F"]
FILE_PS = "data/WideA_ZA_pk.dat"
NAMES_PS = ["k (h/Mpc)", "P (Mpc/h)^3"]
FILE_BACKGROUND = "data/WideA_ZA_background.dat"
NAMES_BACKGROUND = ["z", "proper time [Gyr]", "conf. time [Mpc]", "H [1/Mpc]",
     "comov. dist.", "ang.diam.dist.", "lum. dist.", "comov.snd.hrz.",
     "(.)rho_g", "(.)rho_b", "(.)rho_cdm", "(.)rho_lambda", "(.)rho_ur",
     "(.)rho_crit", "gr.fac. D", "gr.fac. f"]
column_NAMES_BACKGROUND = Dict([x => i for (i, x) in enumerate(NAMES_BACKGROUND)])

include("F_evaluation.jl")
include("Background_functions.jl")
include("Auto_doppler.jl")
include("Auto_lensing.jl")
include("Power_Spectrum.jl")

end # module
