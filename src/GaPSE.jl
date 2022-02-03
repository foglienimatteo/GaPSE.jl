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

const θ_MAX = π / 2.0

const z_MIN = 0.05
const z_MAX = 0.2
const z_EFF = 0.15045636097417317

const s_MIN = 148.1920001343431
const s_MAX = 571.7022420911966
const s_EFF = 435.37470960794167

NAMES_F_MAP = ["x", "mu", "F", "F_error"]
NAMES_PS = ["k (h/Mpc)", "P (Mpc/h)^3"]
NAMES_BACKGROUND = ["z", "proper time [Gyr]", "conf. time [Mpc]", "H [1/Mpc]",
     "comov. dist.", "ang.diam.dist.", "lum. dist.", "comov.snd.hrz.",
     "(.)rho_g", "(.)rho_b", "(.)rho_cdm", "(.)rho_lambda", "(.)rho_ur",
     "(.)rho_crit", "gr.fac. D", "gr.fac. f"]

column_NAMES_BACKGROUND = Dict([x => i for (i, x) in enumerate(NAMES_BACKGROUND)])

include("F_evaluation.jl")
include("Background_functions.jl")
include("Tool_functions.jl")
include("Cosmology.jl")

include("Auto_doppler.jl")
include("Auto_lensing.jl")

IMPLEMENTED_GR_EFFECTS = ["auto_doppler", "auto_lensing"]
dict_gr_mu = Dict(
     "auto_doppler" => integrand_on_mu_doppler,
     "auto_lensing" => integrand_on_mu_lensing,
)

include("Power_Spectrum.jl")
include("Integral_on_mu.jl")

function parameters_used(io::IO, cosmo::Cosmology)
     println(io, "# The following parameters were used for this computation: ")
     println(io, "# CLASS Power Spectrum input file : \"$(cosmo.file_ips)\"")
     println(io, "# F window function input file : \"$(cosmo.file_windowF)\"")
     println(io, "# CLASS Background input file: \"$(cosmo.file_data)\"")
     println(io, "# \t z_min = $(cosmo.params.z_min) \t z_max = $(cosmo.params.z_max)")
     println(io, "# \t k_min = $(cosmo.params.k_min) \t k_max = $(cosmo.params.k_max)")
     println(io, "# \t h_0 = $(cosmo.params.h_0) \t Ω_b = $(cosmo.params.Ω_b) \t " *
                 "Ω_cdm = $(cosmo.params.Ω_cdm) \t Ω_M0 = $(cosmo.params.Ω_M0)")
     println(io, "# \t comoving s_min = $(cosmo.s_min) Mpc/h_0")
     println(io, "# \t comoving s_max = $(cosmo.s_max) Mpc/h_0")
     println(io, "# \t comoving s_eff = $(cosmo.s_eff) Mpc/h_0")
     println(io, "# \t comoving z_eff = $(cosmo.z_eff) ")
     println(io, "# \t Volume of the survey V_survey = $(cosmo.volume)")
     println(io, "# \t σ_0 = $(cosmo.tools.σ_0)")
     println(io, "# \t σ_1 = $(cosmo.tools.σ_1)")
     println(io, "# \t σ_2 = $(cosmo.tools.σ_2)")
     println(io, "# \t σ_3 = $(cosmo.tools.σ_0)")
     println(io, "# ")
end


end # module