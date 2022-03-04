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


using GaPSE, Test
using Dierckx, DelimitedFiles, QuadGK, Suppressor


const FILE_F_MAP = "datatest/F_REFERENCE.txt"
const FILE_PS = "datatest/file_pk.txt"
#const FILE_ILN = "datatest/tab_xi.txt"
const FILE_ILN = "datatest/table_Iln.txt"
const FILE_BACKGROUND = "datatest/WideA_ZA_background.dat"

const Z_MIN = 0.05
const Z_MAX = 0.20
const S_MIN = 148.1920001343431
const S_MAX = 571.7022420911966
const S_EFF = 435.37470960794167
const Z_EFF = 0.15045636097417317
const VOLUME = 3.845366167993746e8
const HUBBLE_0 = 1e5 / 299792458.0

include("TEST_DATA.jl")

##########################################################################################92

@testset "test MathUtils" begin
     include("test_MathUtils.jl")
end

@testset "test IPSTools" begin
     include("test_IPSTools.jl")
end

@testset "test F_evaluation" begin
     include("test_F_evaluation.jl")
end

@testset "test BackgroundData" begin
     include("test_BackgroundData.jl")
end

@testset "test Cosmology" begin
     include("test_Cosmology.jl")
end


##########################################################################################92

const PARAMS = GaPSE.CosmoParams(Z_MIN, Z_MAX, π / 2.0;
     k_min = 1e-8, k_max = 10.0,
     Ω_b = 0.0489, Ω_cdm = 0.251020, h_0 = 0.70,
     N = 1024, fit_min = 0.05, fit_max = 0.5, con = true)

const COSMO = GaPSE.Cosmology(PARAMS, FILE_BACKGROUND, FILE_PS, FILE_F_MAP)


common_kwargs = Dict(
     :pr => false,
     :use_windows => false,
     :enhancer => 1e8, :N_μs => 30,
     :μ_atol => 0.0, :μ_rtol => 1e-2,
     :N_log => 100,
);

spec_effect = [
     "auto_lensing", "auto_integratedgp",
     "lensing_doppler", "doppler_lensing",
     "doppler_integratedgp", "integratedgp_doppler",
     "lensing_localgp", "localgp_lensing",
     "lensing_integratedgp", "integratedgp_lensing",
     "localgp_integratedgp", "integratedgp_localgp",
];

specific_kwargs = [effect ∈ spec_effect ? Dict(
     :en => 1e12, :N_χs => 30) : nothing for effect in GaPSE.IMPLEMENTED_GR_EFFECTS]

joint_kwargs = [isnothing(spec) ? common_kwargs : merge(common_kwargs, spec)
                for spec in specific_kwargs];

dict_L_dir = Dict(0 => "monopoles", 1 => "dipoles", 2 => "tripoles",
     3 => "quadrupoles", 4 => "pentapoles")




##############################


@testset "test XiMultipoles" begin
     include("test_XiMultipoles.jl")
end

@testset "test PowerSpectrum" begin
     include("test_PowerSpectrum.jl")
end

@testset "test SumXiMultipoles" begin
     include("test_SumXiMultipoles.jl")
end

##############################

@testset "test AutoDoppler" begin
     include("test_AutoCorrelations/test_AutoDoppler.jl")
end

@testset "test AutoIntegatedGP" begin
     include("test_AutoCorrelations/test_AutoIntegratedGP.jl")
end

@testset "test AutoLocalGP" begin
     include("test_AutoCorrelations/test_AutoLocalGP.jl")
end

@testset "test AutoLensing" begin
     include("test_AutoCorrelations/test_AutoLensing.jl")
end


##############################



@testset "test DopplerLensing" begin
     include("test_CrossCorrelations/test_DopplerLensing.jl")
end

@testset "test DopplerLocalGP" begin
     include("test_CrossCorrelations/test_DopplerLocalGP.jl")
end

@testset "test DopplerIntegratedGP" begin
     include("test_CrossCorrelations/test_DopplerIntegratedGP.jl")
end

@testset "test LensingIntegratedGP" begin
     include("test_CrossCorrelations/test_LensingIntegratedGP.jl")
end

@testset "test LensingLocalGP" begin
     include("test_CrossCorrelations/test_LensingLocalGP.jl")
end

@testset "test LocalGPIntegratedGP" begin
     include("test_CrossCorrelations/test_LocalGPIntegratedGP.jl")
end


##############################