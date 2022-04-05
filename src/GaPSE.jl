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


BRAND =
"""
###############
#    GaPSE    #
############### \n#"""

NAMES_F_MAP = ["x", "mu", "F", "F_error"]
NAMES_PS = ["k (h/Mpc)", "P (Mpc/h)^3"]
NAMES_BACKGROUND = ["z", "proper time [Gyr]", "conf. time [Mpc]", "H [1/Mpc]",
     "comov. dist.", "ang.diam.dist.", "lum. dist.", "comov.snd.hrz.",
     "(.)rho_g", "(.)rho_b", "(.)rho_cdm", "(.)rho_lambda", "(.)rho_ur",
     "(.)rho_crit", "gr.fac. D", "gr.fac. f"]

include("OtherUtils.jl")
include("MathUtils.jl")
include("WindowF.jl")
include("BackgroundData.jl")
include("CosmoParams.jl")
include("CosmoUtils.jl")
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

include("PP_Doppler.jl")

include("Dicts.jl")

include("PowerSpectrum.jl")
include("XiMultipoles.jl")
include("SumXiMultipoles.jl")



##########################################################################################92



function parameters_used(io::IO, cosmo::Cosmology)
     println(io, BRAND)
     println(io, "# The Cosmology considered had the following paremeters:\n# ")
     println(io, "# - Matter Power Spectrum input file: \"$(cosmo.file_ips)\"")
     println(io, "# - F window function input file: \"$(cosmo.file_windowF)\"")
     println(io, "# - Background data input file: \"$(cosmo.file_data)\"")
     println(io, "#")

     println(io, "# - Basic CosmoParams considered: ")
     println(io, "#\t z_min = $(cosmo.params.z_min) \t z_max = $(cosmo.params.z_max)")
     println(io, "#\t θ_max = $(cosmo.params.θ_max) [rad] \t h_0 = $(cosmo.params.h_0)")
     println(io, "#\t Ω_b = $(cosmo.params.Ω_b) \t " *
                 "Ω_cdm = $(cosmo.params.Ω_cdm) \t Ω_M0 = $(cosmo.params.Ω_M0)")
     println(io, "#")
    
     println(io, "# - CosmoParams about the Input Power Spectrum: ")
     my_println_dict(io, cosmo.params.IPS; pref ="#\t ", N = 2)
     println(io, "#")
    
     println(io, "# - CosmoParams about the Input Power Spectrum Tools: ")
     my_println_dict(io, cosmo.params.IPSTools; pref ="#\t ", N = 3)
     println(io, "#") 

     println(io, "# - Computed quantities: ")
     println(io, "# \t effective redshift z_eff = $(cosmo.z_eff) ")
     println(io, "# \t comoving s_min = $(cosmo.s_min) Mpc/h_0")
     println(io, "# \t comoving s_max = $(cosmo.s_max) Mpc/h_0")
     println(io, "# \t comoving s_eff = $(cosmo.s_eff) Mpc/h_0")
     println(io, "# \t Volume of the survey V_survey = $(cosmo.volume) Mpc^3/h_0^3")
     println(io, "# \t σ_0 = $(cosmo.tools.σ_0)")
     println(io, "# \t σ_1 = $(cosmo.tools.σ_1)")
     println(io, "# \t σ_2 = $(cosmo.tools.σ_2)")
     println(io, "# \t σ_3 = $(cosmo.tools.σ_3)")
     println(io, "# \t (where σ_i = \\int_{k_{min}}^{k_{max}}\\frac{dq}{2 π^2} q^{2-i} P(q))")
     println(io, "# ")
end

parameters_used(cosmo::Cosmology) = parameters_used(stdout, cosmo)


"""
     parameters_used(io::IO, cosmo::Cosmology)
     parameters_used(cosmo::Cosmology) = parameters_used(stdout, cosmo)

Writes in the `io` stream all the information concerning the input Cosmology `cosmo`.

See also: [`Cosmology`](@ref)
"""
parameters_used


end # module
