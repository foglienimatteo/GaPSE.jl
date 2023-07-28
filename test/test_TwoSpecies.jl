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


println("\n\nNow I start the tests on a new Cosmology with two galaxy species, i.e.\n" *
        "  different values of b, s_b and f_evo between the first and the second effect in \n" *
        "  the TPCF evaluation; we will use again the functions print_map_sum_ξ_GNC_multipole, \n"*
        "  print_map_sum_ξ_GNCxLD_multipole, print_map_sum_ξ_LDxGNC_multipole.\n")

params_2species = GaPSE.CosmoParams(Z_MIN, Z_MAX, π / 2.0;
    Ω_b=0.0489, Ω_cdm=0.251020, h_0=0.7, s_lim=1e-2,
    b1=1.0, b2=3.0, s_b1=0.3, s_b2=2.0, 𝑓_evo1=30.0, 𝑓_evo2=40.0,
    IPS_opts=Dict(
        :fit_left_min => 1e-6, :fit_left_max => 3e-6,
        :fit_right_min => 1e1, :fit_right_max => 2e1,
    ),
    IPSTools_opts=Dict(
        :N => 1024, :fit_min => 0.05, :fit_max => 0.5,
        :con => true, :k_min => 1e-8, :k_max => 10.0
    )
    #=
     WFI_opts = Dict(:llim => 0.0, :rlim => Inf, :N => 1000,
          :trap => true, :rtol => 1e-2, :atol => 0.0,
          :ss_start => 0.0, :ss_step => 21.768735478453323,
          :ss_stop => 0.0)
    =#
)

const cosmo_2species = GaPSE.Cosmology(params_2species, FILE_BACKGROUND, FILE_PS, FILE_F_MAP, FILE_IF_MAP);
const ss_2species = 10 .^ range(0, log10(2.0 * cosmo_2species.s_max), length=100);

@testset "test print_map_sum_ξ_GNC_multipole L=0 no_window with obs" begin
    RTOL = 1.2e-2

    kwargs = Dict(
        :pr => false, :obs => :yes,
        :use_windows => false,
        :enhancer => 1e8,
        :N_trap => 30, :N_lob => 30,
        :atol_quad => 0.0, :rtol_quad => 1e-2,
        :N_χs => 100, :N_χs_2 => 60,
        :N_log => 3,
    )


    ##### quad ####
    ss_quad, res_sums_quad, res_xis_quad = GaPSE.readxyall(
        "datatest/TwoSpecies/xis_2species_GNC_L0_noF_withobs.txt", comments=true)

    name = "calc_xis_GNC_L0_quad_noF_withobs.txt"
    isfile(name) && rm(name)

    GaPSE.print_map_sum_ξ_GNC_multipole(cosmo_2species, name, ss_2species;
        single=true, L=0, alg=:quad, suit_sampling = true, kwargs...)
    calc_ss_quad, calc_sums_quad, calc_xis_quad = GaPSE.readxyall(name, comments=true)

    @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[1], calc_xis_quad[1])]) # auto_newton
    @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[2], calc_xis_quad[2])]) # auto_doppler
    @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[3], calc_xis_quad[3])]) # auto_lensing
    @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[4], calc_xis_quad[4])]) # auto_localgp
    @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[5], calc_xis_quad[5])]) # auto_integrated
    @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[6], calc_xis_quad[6])]) # newton_doppler
    @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[7], calc_xis_quad[7])]) # doppler_newton
    @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[8], calc_xis_quad[8])]) # newton_lensing
    @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[9], calc_xis_quad[9])]) # lensing_newton
    @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[10], calc_xis_quad[10])]) # newton_localgp
    @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[11], calc_xis_quad[11])]) # localgp_newton
    @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[12], calc_xis_quad[12])]) # newton_integratedgp
    @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[13], calc_xis_quad[13])]) # integratedgp_newton
    @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[14], calc_xis_quad[14])]) # lensing_doppler
    @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[15], calc_xis_quad[15])]) # doppler_lensing
    @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[16], calc_xis_quad[16])]) # doppler_localgp
    @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[17], calc_xis_quad[17])]) # localgp_doppler
    @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[18], calc_xis_quad[18])]) # doppler_integratedgp
    @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[19], calc_xis_quad[19])]) # integratedgp_doppler
    @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[20], calc_xis_quad[20])]) # lensing_localgp
    @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[21], calc_xis_quad[21])]) # localgp_lensing
    @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[22], calc_xis_quad[22])]) # lensing_integratedgp
    @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[23], calc_xis_quad[23])]) # integratedgp_lensing
    @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[24], calc_xis_quad[24])]) # localgp_integratedgp
    @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[25], calc_xis_quad[25])]) # integratedgp_localgp

    rm(name)
end

println("\n 1 / 3 ...")

@testset "test print_map_sum_ξ_GNCxLD_multipole L=0 no_window " begin
    RTOL = 1.2e-2

    kwargs = Dict(
        :pr => false,
        :use_windows => false,
        :enhancer => 1e8,
        :N_trap => 30, :N_lob => 30,
        :atol_quad => 0.0, :rtol_quad => 1e-2,
        :N_χs => 100, :N_χs_2 => 60,
        :N_log => 3,
    )


    ##### quad ####
    ss_quad, res_sums_quad, res_xis_quad = GaPSE.readxyall(
        "datatest/TwoSpecies/xis_2species_GNCxLD_L0_noF.txt", comments=true)

    name = "calc_xis_GNCxLD_L0_quad_noF.txt"
    isfile(name) && rm(name)

    GaPSE.print_map_sum_ξ_GNCxLD_multipole(cosmo_2species, name, ss_2species;
        single=true, L=0, alg=:quad, kwargs...)
    calc_ss_quad, calc_sums_quad, calc_xis_quad = GaPSE.readxyall(name, comments=true)

    @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[1], calc_xis_quad[1])]) # newton_doppler
    @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[2], calc_xis_quad[2])]) # newton_lensing
    @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[3], calc_xis_quad[3])]) # newton_localgp
    @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[4], calc_xis_quad[4])]) # newton_integratedgp
    @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[5], calc_xis_quad[5])]) # doppler_doppler
    @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[6], calc_xis_quad[6])]) # doppler_lensing
    @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[7], calc_xis_quad[7])]) # doppler_localgp
    @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[8], calc_xis_quad[8])]) # doppler_integratedgp
    @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[9], calc_xis_quad[9])]) # lensing_doppler
    @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[10], calc_xis_quad[10])]) # lensing_lensing
    @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[11], calc_xis_quad[11])]) # lensing_localgp
    @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[12], calc_xis_quad[12])]) # lensing_integratedgp
    @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[13], calc_xis_quad[13])]) # localgp_doppler
    @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[14], calc_xis_quad[14])]) # localgp_lensing
    @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[15], calc_xis_quad[15])]) # localgp_localgp
    @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[16], calc_xis_quad[16])]) # localgp_integratedgp
    @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[17], calc_xis_quad[17])]) # integratedgp_doppler
    @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[18], calc_xis_quad[18])]) # integratedgp_lensing
    @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[19], calc_xis_quad[19])]) # integratedgp_localgp
    @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[20], calc_xis_quad[20])]) # integratedgp_integratedgp

    rm(name)
end

println("\n 2 / 3 ...")

@testset "test print_map_sum_ξ_LDxGNC_multipole L=0 no_window " begin
    RTOL = 1.2e-2

    kwargs = Dict(
        :pr => false,
        :use_windows => false,
        :enhancer => 1e8,
        :N_trap => 30, :N_lob => 30,
        :atol_quad => 0.0, :rtol_quad => 1e-2,
        :N_χs => 100, :N_χs_2 => 60,
        :N_log => 3,
    )


    ##### quad ####
    ss_quad, res_sums_quad, res_xis_quad = GaPSE.readxyall(
        "datatest/TwoSpecies/xis_2species_LDxGNC_L0_noF.txt", comments=true)

    name = "calc_xis_LDxGNC_L0_quad_noF.txt"
    isfile(name) && rm(name)

    GaPSE.print_map_sum_ξ_LDxGNC_multipole(cosmo_2species, name, ss_2species;
        single=true, L=0, alg=:quad, kwargs...)
    calc_ss_quad, calc_sums_quad, calc_xis_quad = GaPSE.readxyall(name, comments=true)

    @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[1], calc_xis_quad[1])]) # doppler_newton
    @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[2], calc_xis_quad[2])]) # lensing_newton
    @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[3], calc_xis_quad[3])]) # localgp_newton
    @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[4], calc_xis_quad[4])]) # integratedgp_newton
    @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[5], calc_xis_quad[5])]) # doppler_doppler
    @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[6], calc_xis_quad[6])]) # lensing_doppler
    @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[7], calc_xis_quad[7])]) # localgp_doppler
    @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[8], calc_xis_quad[8])]) # integratedgp_doppler
    @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[9], calc_xis_quad[9])]) # doppler_lensing
    @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[10], calc_xis_quad[10])]) # lensing_lensing
    @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[11], calc_xis_quad[11])]) # localgp_lensing
    @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[12], calc_xis_quad[12])]) # integratedgp_lensing
    @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[13], calc_xis_quad[13])]) # doppler_localgp
    @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[14], calc_xis_quad[14])]) # lensing_localgp
    @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[15], calc_xis_quad[15])]) # localgp_localgp
    @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[16], calc_xis_quad[16])]) # integratedgp_localgp
    @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[17], calc_xis_quad[17])]) # doppler_integratedgp
    @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[18], calc_xis_quad[18])]) # lensing_integratedgp
    @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[19], calc_xis_quad[19])]) # localgp_integratedgp
    @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[20], calc_xis_quad[20])]) # integratedgp_integratedgp

    rm(name)
end

println("\n 3 / 3 ...")

println("\nEnded tests on 2 galaxies species!\n")
