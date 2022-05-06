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


const FILE_F_MAP = "datatest/F_REFERENCE_pi2.txt"
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


@testset "test OtherUtils" begin
     include("test_OtherUtils.jl")
end

@testset "test MathUtils" begin
     include("test_MathUtils.jl")
end

@testset "test CosmoUtils" begin
     include("test_CosmoUtils.jl")
end

@testset "test IPSTools" begin
     include("test_IPSTools.jl")
end

@testset "test BackgroundData" begin
     include("test_BackgroundData.jl")
end

@testset "test CosmoParams" begin
     include("test_CosmoParams.jl")
end

@testset "test Cosmology" begin
     include("test_Cosmology.jl")
end

@testset "test WindowF" begin
     include("test_WindowF.jl")
end


##########################################################################################92


const PARAMS = GaPSE.CosmoParams(Z_MIN, Z_MAX, π / 2.0;
     Ω_b=0.0489, Ω_cdm=0.251020, h_0=0.7, s_lim=1e-2,
     IPS_opts=Dict(
          :fit_left_min => 1e-6, :fit_left_max => 3e-6,
          :fit_right_min => 1e1, :fit_right_max => 2e1,
     ),
     IPSTools_opts=Dict(
          :N => 1024, :fit_min => 0.05, :fit_max => 0.5,
          :con => true, :k_min => 1e-8, :k_max => 10.0)
)
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
     :en => 1e12, :N_χs => 30) : nothing for effect in GaPSE.GR_EFFECTS_LD]

joint_kwargs = [isnothing(spec) ? common_kwargs : merge(common_kwargs, spec)
                for spec in specific_kwargs];

dict_L_dir = Dict(0 => "monopoles", 1 => "dipoles", 2 => "tripoles",
     3 => "quadrupoles", 4 => "pentapoles")


##############################

@testset "test PowerSpectrum" begin
     include("test_PowerSpectrum.jl")
end

@testset "test LD_XiMultipoles" begin
     include("test_LD_XiMultipoles.jl")
end

@testset "test LD_SumXiMultipoles" begin
     include("test_LD_SumXiMultipoles.jl")
end

##############################

@testset "test LD_AutoDoppler" begin
     include("test_LD_AutoCorrelations/test_LD_AutoDoppler.jl")
end

@testset "test LD_AutoIntegratedGP" begin
     include("test_LD_AutoCorrelations/test_LD_AutoIntegratedGP.jl")
end

@testset "test LD_AutoLocalGP" begin
     include("test_LD_AutoCorrelations/test_LD_AutoLocalGP.jl")
end

@testset "test LD_AutoLensing" begin
     include("test_LD_AutoCorrelations/test_LD_AutoLensing.jl")
end


##############################



@testset "test LD_DopplerLensing" begin
     include("test_LD_CrossCorrelations/test_LD_DopplerLensing.jl")
end

@testset "test LD_DopplerLocalGP" begin
     include("test_LD_CrossCorrelations/test_LD_DopplerLocalGP.jl")
end

@testset "test LD_DopplerIntegratedGP" begin
     include("test_LD_CrossCorrelations/test_LD_DopplerIntegratedGP.jl")
end

@testset "test LD_LensingIntegratedGP" begin
     include("test_LD_CrossCorrelations/test_LD_LensingIntegratedGP.jl")
end

@testset "test LD_LensingLocalGP" begin
     include("test_LD_CrossCorrelations/test_LD_LensingLocalGP.jl")
end

@testset "test LD_LocalGPIntegratedGP" begin
     include("test_LD_CrossCorrelations/test_LD_LocalGPIntegratedGP.jl")
end


##############################


@testset "test PPDoppler" begin
     include("test_PPDoppler.jl")
end

