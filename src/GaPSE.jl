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

FILE_F_MAP = "data/F_REFERENCE.txt"
NAMES_F_MAP = ["x", "mu", "F", "F_error"]

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
include("Tool_functions.jl")
include("Auto_doppler.jl")
include("Auto_lensing.jl")
include("Power_Spectrum.jl")

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
     println(io, "# \t Volume of the survey V_survey = $(@sprintf("%.6e", V_survey()))")
     println(io, "# \t σ_0 = $(@sprintf("%.3e", σ_0)) \t σ_1 = $(@sprintf("%.3e", σ_1)) \t ")
     println(io, "# \t σ_2 = $(@sprintf("%.3e", σ_2)) \t σ_3 = $(@sprintf("%.3e", σ_3)) \t ")
     println(io, "# ")
end

parameters_used(stdout)

end # module
