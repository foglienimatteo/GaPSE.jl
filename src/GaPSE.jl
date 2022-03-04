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
using HCubature, QuadGK, LegendrePolynomials # Licence: MIT "Expat"
using SpecialFunctions, Trapz, LsqFit # Licence: MIT
using GridInterpolations  # Licence: MIT "Expat"
using ProgressMeter, Documenter  # Licence: MIT "Expat"
using Test, Printf, DelimitedFiles  # Licence: MIT "Expat"

#=
const θ_MAX = π / 2.0

const z_MIN = 0.05
const z_MAX = 0.2
const z_EFF = 0.15045636097417317

const s_MIN = 148.1920001343431
const s_MAX = 571.7022420911966
const s_EFF = 435.37470960794167
=#

NAMES_F_MAP = ["x", "mu", "F", "F_error"]
NAMES_PS = ["k (h/Mpc)", "P (Mpc/h)^3"]
NAMES_BACKGROUND = ["z", "proper time [Gyr]", "conf. time [Mpc]", "H [1/Mpc]",
     "comov. dist.", "ang.diam.dist.", "lum. dist.", "comov.snd.hrz.",
     "(.)rho_g", "(.)rho_b", "(.)rho_cdm", "(.)rho_lambda", "(.)rho_ur",
     "(.)rho_crit", "gr.fac. D", "gr.fac. f"]

include("MathUtils.jl")
include("F_evaluation.jl")
include("BackgroundData.jl")
include("IPSTools.jl")
include("Cosmology.jl")

include("AutoCorrelations/AutoDoppler.jl")
include("AutoCorrelations/AutoLensing.jl")
include("AutoCorrelations/AutoLocalGP.jl")
include("AutoCorrelations/AutoIntegratedGP.jl")

include("CrossCorrelations/LensingDoppler.jl")
include("CrossCorrelations/DopplerLocalGP.jl")
include("CrossCorrelations/DopplerIntegratedGP.jl")
include("CrossCorrelations/LensingLocalGP.jl")
include("CrossCorrelations/LensingIntegratedGP.jl")
include("CrossCorrelations/LocalGPIntegratedGP.jl")


IMPLEMENTED_GR_EFFECTS = [
     "auto_doppler", "auto_lensing",
     "auto_localgp", "auto_integratedgp", "lensing_doppler", "doppler_lensing",
     "doppler_localgp", "localgp_doppler",
     "doppler_integratedgp", "integratedgp_doppler",
     "lensing_localgp", "localgp_lensing",
     "lensing_integratedgp", "integratedgp_lensing",
     "localgp_integratedgp", "integratedgp_localgp",
]

IMPLEMENTED_ξs = [
     ξ_Doppler, ξ_Lensing, ξ_LocalGP, ξ_IntegratedGP, ξ_Lensing_Doppler, ξ_Doppler_Lensing,
     ξ_Doppler_LocalGP, ξ_LocalGP_Doppler,
     ξ_Doppler_IntegratedGP, ξ_IntegratedGP_Doppler,
     ξ_Lensing_LocalGP, ξ_LocalGP_Lensing,
     ξ_Lensing_IntegratedGP, ξ_IntegratedGP_Lensing,
     ξ_LocalGP_IntegratedGP, ξ_IntegratedGP_LocalGP
]

DICT_GR_ξs = Dict([k => v for (k, v) in zip(IMPLEMENTED_GR_EFFECTS, IMPLEMENTED_ξs)]...)
INDEX_GR_EFFECT = Dict([name => i for (i, name) in
                        enumerate(IMPLEMENTED_GR_EFFECTS)]...)
GR_EFFECT_INDEXED = Dict([i => name for (i, name) in
                          enumerate(IMPLEMENTED_GR_EFFECTS)]...)


include("PowerSpectrum.jl")
include("XiMultipoles.jl")
include("SumXiMultipoles.jl")


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
     println(io, "# \t σ_3 = $(cosmo.tools.σ_3)")
     println(io, "# ")
end


end # module
