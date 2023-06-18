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


println("\nNow I work on print_map_sum_ξ_GNC_multipole.")
println("This will be the longest test...")

@testset "test print_map_sum_ξ_GNC_multipole no_window no obs, first method" begin
    RTOL = 1.2e-2

    kwargs = Dict(
        :use_windows => false, :obs => :no,
        :enhancer => 1e8, :s1 => nothing,
        :N_trap => 30, :N_lob => 30,
        :atol_quad => 0.0, :rtol_quad => 1e-2,
        :N_χs => 40, :N_χs_2 => 20,
        :pr => false,
    )

    @testset "test_zeros" begin
        @test_throws AssertionError GaPSE.GaPSE.print_map_sum_ξ_GNC_multipole(COSMO, "NONEXISTINGDIR/file.txt", [1.0, 2.0, 3.0])
    end

    @testset "L = 0 no_window no obs" begin
        L = 0

        ##### quad ####
        ss_quad, res_sums_quad, res_xis_quad = GaPSE.readxyall(
            "datatest/GNC_SumXiMultipoles/xis_GNC_L$L" * "_quad_noF_noobs_specific_ss.txt", comments=true)

        name = "calc_xis_GNC_L$L" * "_quad_noF_noobs_specific_ss.txt"
        isfile(name) && rm(name)

        GaPSE.print_map_sum_ξ_GNC_multipole(COSMO, name, ss_quad;
            single=true, L=L, alg=:quad, kwargs...)
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

        #### trap ####
        ss_trap, res_sums_trap, res_xis_trap = GaPSE.readxyall(
            "datatest/GNC_SumXiMultipoles/xis_GNC_L$L" * "_trap_noF_noobs_specific_ss.txt", comments=true)

        name = "calc_xis_GNC_L$L" * "_trap_noF_noobs_specific_ss.txt"
        isfile(name) && rm(name)

        GaPSE.print_map_sum_ξ_GNC_multipole(COSMO, name, ss_trap;
            single=true, L=L, alg=:trap, kwargs...)
        calc_ss_trap, calc_sums_trap, calc_xis_trap = GaPSE.readxyall(name, comments=true)

        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(ss_trap, calc_ss_trap)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_sums_trap, calc_sums_trap)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[1], calc_xis_trap[1])]) # auto_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[2], calc_xis_trap[2])]) # auto_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[3], calc_xis_trap[3])]) # auto_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[4], calc_xis_trap[4])]) # auto_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[5], calc_xis_trap[5])]) # auto_integrated
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[6], calc_xis_trap[6])]) # newton_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[7], calc_xis_trap[7])]) # doppler_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[8], calc_xis_trap[8])]) # newton_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[9], calc_xis_trap[9])]) # lensing_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[10], calc_xis_trap[10])]) # newton_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[11], calc_xis_trap[11])]) # localgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[12], calc_xis_trap[12])]) # newton_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[13], calc_xis_trap[13])]) # integratedgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[14], calc_xis_trap[14])]) # lensing_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[15], calc_xis_trap[15])]) # doppler_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[16], calc_xis_trap[16])]) # doppler_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[17], calc_xis_trap[17])]) # localgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[18], calc_xis_trap[18])]) # doppler_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[19], calc_xis_trap[19])]) # integratedgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[20], calc_xis_trap[20])]) # lensing_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[21], calc_xis_trap[21])]) # localgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[22], calc_xis_trap[22])]) # lensing_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[23], calc_xis_trap[23])]) # integratedgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[24], calc_xis_trap[24])]) # localgp_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[25], calc_xis_trap[25])]) # integratedgp_localgp

        rm(name)

        #### lob ####
        ss_lob, res_sums_lob, res_xis_lob = GaPSE.readxyall(
            "datatest/GNC_SumXiMultipoles/xis_GNC_L$L" * "_lob_noF_noobs_specific_ss.txt", comments=true)

        name = "calc_xis_GNC_L$L" * "_lob_noF_noobs_specific_ss.txt"
        isfile(name) && rm(name)

        GaPSE.print_map_sum_ξ_GNC_multipole(COSMO, name, ss_lob;
            single=true, L=L, alg=:lobatto, kwargs...)
        calc_ss_lob, calc_sums_lob, calc_xis_lob = GaPSE.readxyall(name, comments=true)

        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(ss_lob, calc_ss_lob)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_sums_lob, calc_sums_lob)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[1], calc_xis_lob[1])]) # auto_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[2], calc_xis_lob[2])]) # auto_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[3], calc_xis_lob[3])]) # auto_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[4], calc_xis_lob[4])]) # auto_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[5], calc_xis_lob[5])]) # auto_integrated
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[6], calc_xis_lob[6])]) # newton_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[7], calc_xis_lob[7])]) # doppler_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[8], calc_xis_lob[8])]) # newton_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[9], calc_xis_lob[9])]) # lensing_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[10], calc_xis_lob[10])]) # newton_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[11], calc_xis_lob[11])]) # localgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[12], calc_xis_lob[12])]) # newton_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[13], calc_xis_lob[13])]) # integratedgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[14], calc_xis_lob[14])]) # lensing_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[15], calc_xis_lob[15])]) # doppler_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[16], calc_xis_lob[16])]) # doppler_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[17], calc_xis_lob[17])]) # localgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[18], calc_xis_lob[18])]) # doppler_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[19], calc_xis_lob[19])]) # integratedgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[20], calc_xis_lob[20])]) # lensing_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[21], calc_xis_lob[21])]) # localgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[22], calc_xis_lob[22])]) # lensing_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[23], calc_xis_lob[23])]) # integratedgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[24], calc_xis_lob[24])]) # localgp_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[25], calc_xis_lob[25])]) # integratedgp_localgp

        rm(name)

    end

    @testset "L = 1 no_window no obs" begin
        L = 1


        #### trap ####
        ss_trap, res_sums_trap, res_xis_trap = GaPSE.readxyall(
            "datatest/GNC_SumXiMultipoles/xis_GNC_L$L" * "_trap_noF_noobs_specific_ss.txt", comments=true)

        name = "calc_xis_GNC_L$L" * "_trap_noF_noobs_specific_ss.txt"
        isfile(name) && rm(name)

        GaPSE.print_map_sum_ξ_GNC_multipole(COSMO, name, ss_trap;
            single=true, L=L, alg=:trap, kwargs...)
        calc_ss_trap, calc_sums_trap, calc_xis_trap = GaPSE.readxyall(name, comments=true)

        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(ss_trap, calc_ss_trap)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_sums_trap, calc_sums_trap)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[1], calc_xis_trap[1])]) # auto_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[2], calc_xis_trap[2])]) # auto_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[3], calc_xis_trap[3])]) # auto_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[4], calc_xis_trap[4])]) # auto_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[5], calc_xis_trap[5])]) # auto_integrated
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[6], calc_xis_trap[6])]) # newton_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[7], calc_xis_trap[7])]) # doppler_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[8], calc_xis_trap[8])]) # newton_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[9], calc_xis_trap[9])]) # lensing_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[10], calc_xis_trap[10])]) # newton_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[11], calc_xis_trap[11])]) # localgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[12], calc_xis_trap[12])]) # newton_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[13], calc_xis_trap[13])]) # integratedgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[14], calc_xis_trap[14])]) # lensing_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[15], calc_xis_trap[15])]) # doppler_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[16], calc_xis_trap[16])]) # doppler_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[17], calc_xis_trap[17])]) # localgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[18], calc_xis_trap[18])]) # doppler_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[19], calc_xis_trap[19])]) # integratedgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[20], calc_xis_trap[20])]) # lensing_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[21], calc_xis_trap[21])]) # localgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[22], calc_xis_trap[22])]) # lensing_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[23], calc_xis_trap[23])]) # integratedgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[24], calc_xis_trap[24])]) # localgp_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[25], calc_xis_trap[25])]) # integratedgp_localgp

        rm(name)

        #### lob ####
        ss_lob, res_sums_lob, res_xis_lob = GaPSE.readxyall(
            "datatest/GNC_SumXiMultipoles/xis_GNC_L$L" * "_lob_noF_noobs_specific_ss.txt", comments=true)

        name = "calc_xis_GNC_L$L" * "_lob_noF_noobs_specific_ss.txt"
        isfile(name) && rm(name)

        GaPSE.print_map_sum_ξ_GNC_multipole(COSMO, name, ss_lob;
            single=true, L=L, alg=:lobatto, kwargs...)
        calc_ss_lob, calc_sums_lob, calc_xis_lob = GaPSE.readxyall(name, comments=true)

        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(ss_lob, calc_ss_lob)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_sums_lob, calc_sums_lob)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[1], calc_xis_lob[1])]) # auto_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[2], calc_xis_lob[2])]) # auto_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[3], calc_xis_lob[3])]) # auto_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[4], calc_xis_lob[4])]) # auto_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[5], calc_xis_lob[5])]) # auto_integrated
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[6], calc_xis_lob[6])]) # newton_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[7], calc_xis_lob[7])]) # doppler_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[8], calc_xis_lob[8])]) # newton_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[9], calc_xis_lob[9])]) # lensing_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[10], calc_xis_lob[10])]) # newton_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[11], calc_xis_lob[11])]) # localgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[12], calc_xis_lob[12])]) # newton_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[13], calc_xis_lob[13])]) # integratedgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[14], calc_xis_lob[14])]) # lensing_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[15], calc_xis_lob[15])]) # doppler_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[16], calc_xis_lob[16])]) # doppler_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[17], calc_xis_lob[17])]) # localgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[18], calc_xis_lob[18])]) # doppler_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[19], calc_xis_lob[19])]) # integratedgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[20], calc_xis_lob[20])]) # lensing_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[21], calc_xis_lob[21])]) # localgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[22], calc_xis_lob[22])]) # lensing_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[23], calc_xis_lob[23])]) # integratedgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[24], calc_xis_lob[24])]) # localgp_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[25], calc_xis_lob[25])]) # integratedgp_localgp

        rm(name)

    end

    @testset "L = 2 no_window no obs" begin
        L = 2

        ##### quad ####
        ss_quad, res_sums_quad, res_xis_quad = GaPSE.readxyall(
            "datatest/GNC_SumXiMultipoles/xis_GNC_L$L" * "_quad_noF_noobs_specific_ss.txt", comments=true)

        name = "calc_xis_GNC_L$L" * "_quad_noF_noobs_specific_ss.txt"
        isfile(name) && rm(name)

        GaPSE.print_map_sum_ξ_GNC_multipole(COSMO, name, ss_quad;
            single=true, L=L, alg=:quad, kwargs...)
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

        #### trap ####
        ss_trap, res_sums_trap, res_xis_trap = GaPSE.readxyall(
            "datatest/GNC_SumXiMultipoles/xis_GNC_L$L" * "_trap_noF_noobs_specific_ss.txt", comments=true)

        name = "calc_xis_GNC_L$L" * "_trap_noF_noobs_specific_ss.txt"
        isfile(name) && rm(name)

        GaPSE.print_map_sum_ξ_GNC_multipole(COSMO, name, ss_trap;
            single=true, L=L, alg=:trap, kwargs...)
        calc_ss_trap, calc_sums_trap, calc_xis_trap = GaPSE.readxyall(name, comments=true)

        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(ss_trap, calc_ss_trap)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_sums_trap, calc_sums_trap)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[1], calc_xis_trap[1])]) # auto_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[2], calc_xis_trap[2])]) # auto_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[3], calc_xis_trap[3])]) # auto_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[4], calc_xis_trap[4])]) # auto_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[5], calc_xis_trap[5])]) # auto_integrated
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[6], calc_xis_trap[6])]) # newton_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[7], calc_xis_trap[7])]) # doppler_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[8], calc_xis_trap[8])]) # newton_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[9], calc_xis_trap[9])]) # lensing_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[10], calc_xis_trap[10])]) # newton_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[11], calc_xis_trap[11])]) # localgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[12], calc_xis_trap[12])]) # newton_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[13], calc_xis_trap[13])]) # integratedgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[14], calc_xis_trap[14])]) # lensing_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[15], calc_xis_trap[15])]) # doppler_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[16], calc_xis_trap[16])]) # doppler_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[17], calc_xis_trap[17])]) # localgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[18], calc_xis_trap[18])]) # doppler_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[19], calc_xis_trap[19])]) # integratedgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[20], calc_xis_trap[20])]) # lensing_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[21], calc_xis_trap[21])]) # localgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[22], calc_xis_trap[22])]) # lensing_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[23], calc_xis_trap[23])]) # integratedgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[24], calc_xis_trap[24])]) # localgp_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[25], calc_xis_trap[25])]) # integratedgp_localgp

        rm(name)

        #### lob ####
        ss_lob, res_sums_lob, res_xis_lob = GaPSE.readxyall(
            "datatest/GNC_SumXiMultipoles/xis_GNC_L$L" * "_lob_noF_noobs_specific_ss.txt", comments=true)

        name = "calc_xis_GNC_L$L" * "_lob_noF_noobs_specific_ss.txt"
        isfile(name) && rm(name)

        GaPSE.print_map_sum_ξ_GNC_multipole(COSMO, name, ss_lob;
            single=true, L=L, alg=:lobatto, kwargs...)
        calc_ss_lob, calc_sums_lob, calc_xis_lob = GaPSE.readxyall(name, comments=true)

        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(ss_lob, calc_ss_lob)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_sums_lob, calc_sums_lob)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[1], calc_xis_lob[1])]) # auto_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[2], calc_xis_lob[2])]) # auto_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[3], calc_xis_lob[3])]) # auto_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[4], calc_xis_lob[4])]) # auto_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[5], calc_xis_lob[5])]) # auto_integrated
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[6], calc_xis_lob[6])]) # newton_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[7], calc_xis_lob[7])]) # doppler_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[8], calc_xis_lob[8])]) # newton_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[9], calc_xis_lob[9])]) # lensing_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[10], calc_xis_lob[10])]) # newton_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[11], calc_xis_lob[11])]) # localgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[12], calc_xis_lob[12])]) # newton_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[13], calc_xis_lob[13])]) # integratedgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[14], calc_xis_lob[14])]) # lensing_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[15], calc_xis_lob[15])]) # doppler_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[16], calc_xis_lob[16])]) # doppler_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[17], calc_xis_lob[17])]) # localgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[18], calc_xis_lob[18])]) # doppler_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[19], calc_xis_lob[19])]) # integratedgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[20], calc_xis_lob[20])]) # lensing_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[21], calc_xis_lob[21])]) # localgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[22], calc_xis_lob[22])]) # lensing_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[23], calc_xis_lob[23])]) # integratedgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[24], calc_xis_lob[24])]) # localgp_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[25], calc_xis_lob[25])]) # integratedgp_localgp

        rm(name)

    end

    @testset "L = 3 no_window no obs" begin
        L = 3


        #### trap ####
        ss_trap, res_sums_trap, res_xis_trap = GaPSE.readxyall(
            "datatest/GNC_SumXiMultipoles/xis_GNC_L$L" * "_trap_noF_noobs_specific_ss.txt", comments=true)

        name = "calc_xis_GNC_L$L" * "_trap_noF_noobs_specific_ss.txt"
        isfile(name) && rm(name)

        GaPSE.print_map_sum_ξ_GNC_multipole(COSMO, name, ss_trap;
            single=true, L=L, alg=:trap, kwargs...)
        calc_ss_trap, calc_sums_trap, calc_xis_trap = GaPSE.readxyall(name, comments=true)

        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(ss_trap, calc_ss_trap)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_sums_trap, calc_sums_trap)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[1], calc_xis_trap[1])]) # auto_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[2], calc_xis_trap[2])]) # auto_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[3], calc_xis_trap[3])]) # auto_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[4], calc_xis_trap[4])]) # auto_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[5], calc_xis_trap[5])]) # auto_integrated
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[6], calc_xis_trap[6])]) # newton_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[7], calc_xis_trap[7])]) # doppler_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[8], calc_xis_trap[8])]) # newton_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[9], calc_xis_trap[9])]) # lensing_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[10], calc_xis_trap[10])]) # newton_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[11], calc_xis_trap[11])]) # localgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[12], calc_xis_trap[12])]) # newton_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[13], calc_xis_trap[13])]) # integratedgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[14], calc_xis_trap[14])]) # lensing_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[15], calc_xis_trap[15])]) # doppler_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[16], calc_xis_trap[16])]) # doppler_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[17], calc_xis_trap[17])]) # localgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[18], calc_xis_trap[18])]) # doppler_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[19], calc_xis_trap[19])]) # integratedgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[20], calc_xis_trap[20])]) # lensing_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[21], calc_xis_trap[21])]) # localgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[22], calc_xis_trap[22])]) # lensing_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[23], calc_xis_trap[23])]) # integratedgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[24], calc_xis_trap[24])]) # localgp_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[25], calc_xis_trap[25])]) # integratedgp_localgp

        rm(name)

        #### lob ####
        ss_lob, res_sums_lob, res_xis_lob = GaPSE.readxyall(
            "datatest/GNC_SumXiMultipoles/xis_GNC_L$L" * "_lob_noF_noobs_specific_ss.txt", comments=true)

        name = "calc_xis_GNC_L$L" * "_lob_noF_noobs_specific_ss.txt"
        isfile(name) && rm(name)

        GaPSE.print_map_sum_ξ_GNC_multipole(COSMO, name, ss_lob;
            single=true, L=L, alg=:lobatto, kwargs...)
        calc_ss_lob, calc_sums_lob, calc_xis_lob = GaPSE.readxyall(name, comments=true)

        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(ss_lob, calc_ss_lob)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_sums_lob, calc_sums_lob)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[1], calc_xis_lob[1])]) # auto_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[2], calc_xis_lob[2])]) # auto_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[3], calc_xis_lob[3])]) # auto_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[4], calc_xis_lob[4])]) # auto_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[5], calc_xis_lob[5])]) # auto_integrated
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[6], calc_xis_lob[6])]) # newton_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[7], calc_xis_lob[7])]) # doppler_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[8], calc_xis_lob[8])]) # newton_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[9], calc_xis_lob[9])]) # lensing_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[10], calc_xis_lob[10])]) # newton_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[11], calc_xis_lob[11])]) # localgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[12], calc_xis_lob[12])]) # newton_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[13], calc_xis_lob[13])]) # integratedgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[14], calc_xis_lob[14])]) # lensing_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[15], calc_xis_lob[15])]) # doppler_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[16], calc_xis_lob[16])]) # doppler_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[17], calc_xis_lob[17])]) # localgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[18], calc_xis_lob[18])]) # doppler_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[19], calc_xis_lob[19])]) # integratedgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[20], calc_xis_lob[20])]) # lensing_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[21], calc_xis_lob[21])]) # localgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[22], calc_xis_lob[22])]) # lensing_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[23], calc_xis_lob[23])]) # integratedgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[24], calc_xis_lob[24])]) # localgp_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[25], calc_xis_lob[25])]) # integratedgp_localgp

        rm(name)

    end

    @testset "L = 4 no_window no obs" begin
        L = 4


        #### trap ####
        ss_trap, res_sums_trap, res_xis_trap = GaPSE.readxyall(
            "datatest/GNC_SumXiMultipoles/xis_GNC_L$L" * "_trap_noF_noobs_specific_ss.txt", comments=true)

        name = "calc_xis_GNC_L$L" * "_trap_noF_noobs_specific_ss.txt"
        isfile(name) && rm(name)

        GaPSE.print_map_sum_ξ_GNC_multipole(COSMO, name, ss_trap;
            single=true, L=L, alg=:trap, kwargs...)
        calc_ss_trap, calc_sums_trap, calc_xis_trap = GaPSE.readxyall(name, comments=true)

        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(ss_trap, calc_ss_trap)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_sums_trap, calc_sums_trap)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[1], calc_xis_trap[1])]) # auto_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[2], calc_xis_trap[2])]) # auto_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[3], calc_xis_trap[3])]) # auto_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[4], calc_xis_trap[4])]) # auto_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[5], calc_xis_trap[5])]) # auto_integrated
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[6], calc_xis_trap[6])]) # newton_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[7], calc_xis_trap[7])]) # doppler_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[8], calc_xis_trap[8])]) # newton_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[9], calc_xis_trap[9])]) # lensing_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[10], calc_xis_trap[10])]) # newton_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[11], calc_xis_trap[11])]) # localgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[12], calc_xis_trap[12])]) # newton_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[13], calc_xis_trap[13])]) # integratedgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[14], calc_xis_trap[14])]) # lensing_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[15], calc_xis_trap[15])]) # doppler_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[16], calc_xis_trap[16])]) # doppler_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[17], calc_xis_trap[17])]) # localgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[18], calc_xis_trap[18])]) # doppler_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[19], calc_xis_trap[19])]) # integratedgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[20], calc_xis_trap[20])]) # lensing_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[21], calc_xis_trap[21])]) # localgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[22], calc_xis_trap[22])]) # lensing_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[23], calc_xis_trap[23])]) # integratedgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[24], calc_xis_trap[24])]) # localgp_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[25], calc_xis_trap[25])]) # integratedgp_localgp

        rm(name)

        #### lob ####
        ss_lob, res_sums_lob, res_xis_lob = GaPSE.readxyall(
            "datatest/GNC_SumXiMultipoles/xis_GNC_L$L" * "_lob_noF_noobs_specific_ss.txt", comments=true)

        name = "calc_xis_GNC_L$L" * "_lob_noF_noobs_specific_ss.txt"
        isfile(name) && rm(name)

        GaPSE.print_map_sum_ξ_GNC_multipole(COSMO, name, ss_lob;
            single=true, L=L, alg=:lobatto, kwargs...)
        calc_ss_lob, calc_sums_lob, calc_xis_lob = GaPSE.readxyall(name, comments=true)

        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(ss_lob, calc_ss_lob)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_sums_lob, calc_sums_lob)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[1], calc_xis_lob[1])]) # auto_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[2], calc_xis_lob[2])]) # auto_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[3], calc_xis_lob[3])]) # auto_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[4], calc_xis_lob[4])]) # auto_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[5], calc_xis_lob[5])]) # auto_integrated
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[6], calc_xis_lob[6])]) # newton_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[7], calc_xis_lob[7])]) # doppler_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[8], calc_xis_lob[8])]) # newton_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[9], calc_xis_lob[9])]) # lensing_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[10], calc_xis_lob[10])]) # newton_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[11], calc_xis_lob[11])]) # localgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[12], calc_xis_lob[12])]) # newton_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[13], calc_xis_lob[13])]) # integratedgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[14], calc_xis_lob[14])]) # lensing_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[15], calc_xis_lob[15])]) # doppler_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[16], calc_xis_lob[16])]) # doppler_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[17], calc_xis_lob[17])]) # localgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[18], calc_xis_lob[18])]) # doppler_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[19], calc_xis_lob[19])]) # integratedgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[20], calc_xis_lob[20])]) # lensing_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[21], calc_xis_lob[21])]) # localgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[22], calc_xis_lob[22])]) # lensing_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[23], calc_xis_lob[23])]) # integratedgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[24], calc_xis_lob[24])]) # localgp_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[25], calc_xis_lob[25])]) # integratedgp_localgp

        rm(name)

    end
end

println("\n 1 / 8 ...")

@testset "test print_map_sum_ξ_GNC_multipole no_window no obs vel, first method" begin
    RTOL = 1.2e-2

    kwargs = Dict(
        :use_windows => false, :obs => :noobsvel,
        :enhancer => 1e8, :s1 => nothing,
        :N_trap => 30, :N_lob => 30,
        :atol_quad => 0.0, :rtol_quad => 1e-2,
        :N_χs => 40, :N_χs_2 => 20,
        :pr => false,
    )

    @testset "test_zeros" begin
        @test_throws AssertionError GaPSE.GaPSE.print_map_sum_ξ_GNC_multipole(COSMO, "NONEXISTINGDIR/file.txt", [1.0, 2.0, 3.0])
    end

    @testset "L = 0 no_window no obs vel" begin
        L = 0

        ##### quad ####
        ss_quad, res_sums_quad, res_xis_quad = GaPSE.readxyall(
            "datatest/GNC_SumXiMultipoles/xis_GNC_L$L" * "_quad_noF_noobsvel_specific_ss.txt", comments=true)

        name = "calc_xis_GNC_L$L" * "_quad_noF_noobsvel_specific_ss.txt"
        isfile(name) && rm(name)

        GaPSE.print_map_sum_ξ_GNC_multipole(COSMO, name, ss_quad;
            single=true, L=L, alg=:quad, kwargs...)
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

        #### trap ####
        ss_trap, res_sums_trap, res_xis_trap = GaPSE.readxyall(
            "datatest/GNC_SumXiMultipoles/xis_GNC_L$L" * "_trap_noF_noobsvel_specific_ss.txt", comments=true)

        name = "calc_xis_GNC_L$L" * "_trap_noF_noobsvel_specific_ss.txt"
        isfile(name) && rm(name)

        GaPSE.print_map_sum_ξ_GNC_multipole(COSMO, name, ss_trap;
            single=true, L=L, alg=:trap, kwargs...)
        calc_ss_trap, calc_sums_trap, calc_xis_trap = GaPSE.readxyall(name, comments=true)

        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(ss_trap, calc_ss_trap)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_sums_trap, calc_sums_trap)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[1], calc_xis_trap[1])]) # auto_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[2], calc_xis_trap[2])]) # auto_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[3], calc_xis_trap[3])]) # auto_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[4], calc_xis_trap[4])]) # auto_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[5], calc_xis_trap[5])]) # auto_integrated
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[6], calc_xis_trap[6])]) # newton_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[7], calc_xis_trap[7])]) # doppler_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[8], calc_xis_trap[8])]) # newton_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[9], calc_xis_trap[9])]) # lensing_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[10], calc_xis_trap[10])]) # newton_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[11], calc_xis_trap[11])]) # localgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[12], calc_xis_trap[12])]) # newton_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[13], calc_xis_trap[13])]) # integratedgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[14], calc_xis_trap[14])]) # lensing_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[15], calc_xis_trap[15])]) # doppler_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[16], calc_xis_trap[16])]) # doppler_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[17], calc_xis_trap[17])]) # localgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[18], calc_xis_trap[18])]) # doppler_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[19], calc_xis_trap[19])]) # integratedgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[20], calc_xis_trap[20])]) # lensing_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[21], calc_xis_trap[21])]) # localgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[22], calc_xis_trap[22])]) # lensing_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[23], calc_xis_trap[23])]) # integratedgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[24], calc_xis_trap[24])]) # localgp_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[25], calc_xis_trap[25])]) # integratedgp_localgp

        rm(name)

        #### lob ####
        ss_lob, res_sums_lob, res_xis_lob = GaPSE.readxyall(
            "datatest/GNC_SumXiMultipoles/xis_GNC_L$L" * "_lob_noF_noobsvel_specific_ss.txt", comments=true)

        name = "calc_xis_GNC_L$L" * "_lob_noF_noobsvel_specific_ss.txt"
        isfile(name) && rm(name)

        GaPSE.print_map_sum_ξ_GNC_multipole(COSMO, name, ss_lob;
            single=true, L=L, alg=:lobatto, kwargs...)
        calc_ss_lob, calc_sums_lob, calc_xis_lob = GaPSE.readxyall(name, comments=true)

        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(ss_lob, calc_ss_lob)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_sums_lob, calc_sums_lob)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[1], calc_xis_lob[1])]) # auto_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[2], calc_xis_lob[2])]) # auto_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[3], calc_xis_lob[3])]) # auto_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[4], calc_xis_lob[4])]) # auto_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[5], calc_xis_lob[5])]) # auto_integrated
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[6], calc_xis_lob[6])]) # newton_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[7], calc_xis_lob[7])]) # doppler_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[8], calc_xis_lob[8])]) # newton_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[9], calc_xis_lob[9])]) # lensing_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[10], calc_xis_lob[10])]) # newton_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[11], calc_xis_lob[11])]) # localgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[12], calc_xis_lob[12])]) # newton_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[13], calc_xis_lob[13])]) # integratedgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[14], calc_xis_lob[14])]) # lensing_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[15], calc_xis_lob[15])]) # doppler_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[16], calc_xis_lob[16])]) # doppler_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[17], calc_xis_lob[17])]) # localgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[18], calc_xis_lob[18])]) # doppler_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[19], calc_xis_lob[19])]) # integratedgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[20], calc_xis_lob[20])]) # lensing_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[21], calc_xis_lob[21])]) # localgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[22], calc_xis_lob[22])]) # lensing_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[23], calc_xis_lob[23])]) # integratedgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[24], calc_xis_lob[24])]) # localgp_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[25], calc_xis_lob[25])]) # integratedgp_localgp

        rm(name)

    end

    @testset "L = 1 no_window no obs vel" begin
        L = 1

        #### trap ####
        ss_trap, res_sums_trap, res_xis_trap = GaPSE.readxyall(
            "datatest/GNC_SumXiMultipoles/xis_GNC_L$L" * "_trap_noF_noobsvel_specific_ss.txt", comments=true)

        name = "calc_xis_GNC_L$L" * "_trap_noF_noobsvel_specific_ss.txt"
        isfile(name) && rm(name)

        GaPSE.print_map_sum_ξ_GNC_multipole(COSMO, name, ss_trap;
            single=true, L=L, alg=:trap, kwargs...)
        calc_ss_trap, calc_sums_trap, calc_xis_trap = GaPSE.readxyall(name, comments=true)

        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(ss_trap, calc_ss_trap)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_sums_trap, calc_sums_trap)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[1], calc_xis_trap[1])]) # auto_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[2], calc_xis_trap[2])]) # auto_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[3], calc_xis_trap[3])]) # auto_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[4], calc_xis_trap[4])]) # auto_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[5], calc_xis_trap[5])]) # auto_integrated
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[6], calc_xis_trap[6])]) # newton_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[7], calc_xis_trap[7])]) # doppler_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[8], calc_xis_trap[8])]) # newton_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[9], calc_xis_trap[9])]) # lensing_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[10], calc_xis_trap[10])]) # newton_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[11], calc_xis_trap[11])]) # localgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[12], calc_xis_trap[12])]) # newton_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[13], calc_xis_trap[13])]) # integratedgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[14], calc_xis_trap[14])]) # lensing_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[15], calc_xis_trap[15])]) # doppler_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[16], calc_xis_trap[16])]) # doppler_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[17], calc_xis_trap[17])]) # localgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[18], calc_xis_trap[18])]) # doppler_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[19], calc_xis_trap[19])]) # integratedgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[20], calc_xis_trap[20])]) # lensing_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[21], calc_xis_trap[21])]) # localgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[22], calc_xis_trap[22])]) # lensing_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[23], calc_xis_trap[23])]) # integratedgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[24], calc_xis_trap[24])]) # localgp_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[25], calc_xis_trap[25])]) # integratedgp_localgp

        rm(name)

        #### lob ####
        ss_lob, res_sums_lob, res_xis_lob = GaPSE.readxyall(
            "datatest/GNC_SumXiMultipoles/xis_GNC_L$L" * "_lob_noF_noobsvel_specific_ss.txt", comments=true)

        name = "calc_xis_GNC_L$L" * "_lob_noF_noobsvel_specific_ss.txt"
        isfile(name) && rm(name)

        GaPSE.print_map_sum_ξ_GNC_multipole(COSMO, name, ss_lob;
            single=true, L=L, alg=:lobatto, kwargs...)
        calc_ss_lob, calc_sums_lob, calc_xis_lob = GaPSE.readxyall(name, comments=true)

        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(ss_lob, calc_ss_lob)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_sums_lob, calc_sums_lob)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[1], calc_xis_lob[1])]) # auto_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[2], calc_xis_lob[2])]) # auto_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[3], calc_xis_lob[3])]) # auto_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[4], calc_xis_lob[4])]) # auto_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[5], calc_xis_lob[5])]) # auto_integrated
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[6], calc_xis_lob[6])]) # newton_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[7], calc_xis_lob[7])]) # doppler_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[8], calc_xis_lob[8])]) # newton_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[9], calc_xis_lob[9])]) # lensing_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[10], calc_xis_lob[10])]) # newton_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[11], calc_xis_lob[11])]) # localgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[12], calc_xis_lob[12])]) # newton_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[13], calc_xis_lob[13])]) # integratedgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[14], calc_xis_lob[14])]) # lensing_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[15], calc_xis_lob[15])]) # doppler_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[16], calc_xis_lob[16])]) # doppler_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[17], calc_xis_lob[17])]) # localgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[18], calc_xis_lob[18])]) # doppler_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[19], calc_xis_lob[19])]) # integratedgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[20], calc_xis_lob[20])]) # lensing_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[21], calc_xis_lob[21])]) # localgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[22], calc_xis_lob[22])]) # lensing_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[23], calc_xis_lob[23])]) # integratedgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[24], calc_xis_lob[24])]) # localgp_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[25], calc_xis_lob[25])]) # integratedgp_localgp

        rm(name)

    end

    @testset "L = 2 no_window no obs vel" begin
        L = 2

        ##### quad ####
        ss_quad, res_sums_quad, res_xis_quad = GaPSE.readxyall(
            "datatest/GNC_SumXiMultipoles/xis_GNC_L$L" * "_quad_noF_noobsvel_specific_ss.txt", comments=true)

        name = "calc_xis_GNC_L$L" * "_quad_noF_noobsvel_specific_ss.txt"
        isfile(name) && rm(name)

        GaPSE.print_map_sum_ξ_GNC_multipole(COSMO, name, ss_quad;
            single=true, L=L, alg=:quad, kwargs...)
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

        #### trap ####
        ss_trap, res_sums_trap, res_xis_trap = GaPSE.readxyall(
            "datatest/GNC_SumXiMultipoles/xis_GNC_L$L" * "_trap_noF_noobsvel_specific_ss.txt", comments=true)

        name = "calc_xis_GNC_L$L" * "_trap_noF_noobsvel_specific_ss.txt"
        isfile(name) && rm(name)

        GaPSE.print_map_sum_ξ_GNC_multipole(COSMO, name, ss_trap;
            single=true, L=L, alg=:trap, kwargs...)
        calc_ss_trap, calc_sums_trap, calc_xis_trap = GaPSE.readxyall(name, comments=true)

        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(ss_trap, calc_ss_trap)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_sums_trap, calc_sums_trap)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[1], calc_xis_trap[1])]) # auto_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[2], calc_xis_trap[2])]) # auto_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[3], calc_xis_trap[3])]) # auto_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[4], calc_xis_trap[4])]) # auto_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[5], calc_xis_trap[5])]) # auto_integrated
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[6], calc_xis_trap[6])]) # newton_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[7], calc_xis_trap[7])]) # doppler_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[8], calc_xis_trap[8])]) # newton_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[9], calc_xis_trap[9])]) # lensing_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[10], calc_xis_trap[10])]) # newton_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[11], calc_xis_trap[11])]) # localgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[12], calc_xis_trap[12])]) # newton_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[13], calc_xis_trap[13])]) # integratedgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[14], calc_xis_trap[14])]) # lensing_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[15], calc_xis_trap[15])]) # doppler_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[16], calc_xis_trap[16])]) # doppler_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[17], calc_xis_trap[17])]) # localgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[18], calc_xis_trap[18])]) # doppler_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[19], calc_xis_trap[19])]) # integratedgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[20], calc_xis_trap[20])]) # lensing_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[21], calc_xis_trap[21])]) # localgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[22], calc_xis_trap[22])]) # lensing_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[23], calc_xis_trap[23])]) # integratedgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[24], calc_xis_trap[24])]) # localgp_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[25], calc_xis_trap[25])]) # integratedgp_localgp

        rm(name)

        #### lob ####
        ss_lob, res_sums_lob, res_xis_lob = GaPSE.readxyall(
            "datatest/GNC_SumXiMultipoles/xis_GNC_L$L" * "_lob_noF_noobsvel_specific_ss.txt", comments=true)

        name = "calc_xis_GNC_L$L" * "_lob_noF_noobsvel_specific_ss.txt"
        isfile(name) && rm(name)

        GaPSE.print_map_sum_ξ_GNC_multipole(COSMO, name, ss_lob;
            single=true, L=L, alg=:lobatto, kwargs...)
        calc_ss_lob, calc_sums_lob, calc_xis_lob = GaPSE.readxyall(name, comments=true)

        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(ss_lob, calc_ss_lob)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_sums_lob, calc_sums_lob)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[1], calc_xis_lob[1])]) # auto_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[2], calc_xis_lob[2])]) # auto_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[3], calc_xis_lob[3])]) # auto_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[4], calc_xis_lob[4])]) # auto_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[5], calc_xis_lob[5])]) # auto_integrated
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[6], calc_xis_lob[6])]) # newton_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[7], calc_xis_lob[7])]) # doppler_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[8], calc_xis_lob[8])]) # newton_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[9], calc_xis_lob[9])]) # lensing_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[10], calc_xis_lob[10])]) # newton_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[11], calc_xis_lob[11])]) # localgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[12], calc_xis_lob[12])]) # newton_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[13], calc_xis_lob[13])]) # integratedgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[14], calc_xis_lob[14])]) # lensing_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[15], calc_xis_lob[15])]) # doppler_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[16], calc_xis_lob[16])]) # doppler_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[17], calc_xis_lob[17])]) # localgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[18], calc_xis_lob[18])]) # doppler_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[19], calc_xis_lob[19])]) # integratedgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[20], calc_xis_lob[20])]) # lensing_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[21], calc_xis_lob[21])]) # localgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[22], calc_xis_lob[22])]) # lensing_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[23], calc_xis_lob[23])]) # integratedgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[24], calc_xis_lob[24])]) # localgp_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[25], calc_xis_lob[25])]) # integratedgp_localgp

        rm(name)

    end

    @testset "L = 3 no_window no obs vel" begin
        L = 3

        #### trap ####
        ss_trap, res_sums_trap, res_xis_trap = GaPSE.readxyall(
            "datatest/GNC_SumXiMultipoles/xis_GNC_L$L" * "_trap_noF_noobsvel_specific_ss.txt", comments=true)

        name = "calc_xis_GNC_L$L" * "_trap_noF_noobsvel_specific_ss.txt"
        isfile(name) && rm(name)

        GaPSE.print_map_sum_ξ_GNC_multipole(COSMO, name, ss_trap;
            single=true, L=L, alg=:trap, kwargs...)
        calc_ss_trap, calc_sums_trap, calc_xis_trap = GaPSE.readxyall(name, comments=true)

        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(ss_trap, calc_ss_trap)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_sums_trap, calc_sums_trap)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[1], calc_xis_trap[1])]) # auto_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[2], calc_xis_trap[2])]) # auto_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[3], calc_xis_trap[3])]) # auto_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[4], calc_xis_trap[4])]) # auto_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[5], calc_xis_trap[5])]) # auto_integrated
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[6], calc_xis_trap[6])]) # newton_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[7], calc_xis_trap[7])]) # doppler_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[8], calc_xis_trap[8])]) # newton_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[9], calc_xis_trap[9])]) # lensing_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[10], calc_xis_trap[10])]) # newton_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[11], calc_xis_trap[11])]) # localgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[12], calc_xis_trap[12])]) # newton_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[13], calc_xis_trap[13])]) # integratedgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[14], calc_xis_trap[14])]) # lensing_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[15], calc_xis_trap[15])]) # doppler_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[16], calc_xis_trap[16])]) # doppler_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[17], calc_xis_trap[17])]) # localgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[18], calc_xis_trap[18])]) # doppler_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[19], calc_xis_trap[19])]) # integratedgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[20], calc_xis_trap[20])]) # lensing_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[21], calc_xis_trap[21])]) # localgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[22], calc_xis_trap[22])]) # lensing_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[23], calc_xis_trap[23])]) # integratedgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[24], calc_xis_trap[24])]) # localgp_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[25], calc_xis_trap[25])]) # integratedgp_localgp

        rm(name)

        #### lob ####
        ss_lob, res_sums_lob, res_xis_lob = GaPSE.readxyall(
            "datatest/GNC_SumXiMultipoles/xis_GNC_L$L" * "_lob_noF_noobsvel_specific_ss.txt", comments=true)

        name = "calc_xis_GNC_L$L" * "_lob_noF_noobsvel_specific_ss.txt"
        isfile(name) && rm(name)

        GaPSE.print_map_sum_ξ_GNC_multipole(COSMO, name, ss_lob;
            single=true, L=L, alg=:lobatto, kwargs...)
        calc_ss_lob, calc_sums_lob, calc_xis_lob = GaPSE.readxyall(name, comments=true)

        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(ss_lob, calc_ss_lob)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_sums_lob, calc_sums_lob)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[1], calc_xis_lob[1])]) # auto_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[2], calc_xis_lob[2])]) # auto_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[3], calc_xis_lob[3])]) # auto_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[4], calc_xis_lob[4])]) # auto_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[5], calc_xis_lob[5])]) # auto_integrated
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[6], calc_xis_lob[6])]) # newton_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[7], calc_xis_lob[7])]) # doppler_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[8], calc_xis_lob[8])]) # newton_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[9], calc_xis_lob[9])]) # lensing_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[10], calc_xis_lob[10])]) # newton_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[11], calc_xis_lob[11])]) # localgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[12], calc_xis_lob[12])]) # newton_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[13], calc_xis_lob[13])]) # integratedgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[14], calc_xis_lob[14])]) # lensing_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[15], calc_xis_lob[15])]) # doppler_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[16], calc_xis_lob[16])]) # doppler_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[17], calc_xis_lob[17])]) # localgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[18], calc_xis_lob[18])]) # doppler_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[19], calc_xis_lob[19])]) # integratedgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[20], calc_xis_lob[20])]) # lensing_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[21], calc_xis_lob[21])]) # localgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[22], calc_xis_lob[22])]) # lensing_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[23], calc_xis_lob[23])]) # integratedgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[24], calc_xis_lob[24])]) # localgp_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[25], calc_xis_lob[25])]) # integratedgp_localgp

        rm(name)

    end

    @testset "L = 4 no_window no obs vel" begin
        L = 4

        #### trap ####
        ss_trap, res_sums_trap, res_xis_trap = GaPSE.readxyall(
            "datatest/GNC_SumXiMultipoles/xis_GNC_L$L" * "_trap_noF_noobsvel_specific_ss.txt", comments=true)

        name = "calc_xis_GNC_L$L" * "_trap_noF_noobsvel_specific_ss.txt"
        isfile(name) && rm(name)

        GaPSE.print_map_sum_ξ_GNC_multipole(COSMO, name, ss_trap;
            single=true, L=L, alg=:trap, kwargs...)
        calc_ss_trap, calc_sums_trap, calc_xis_trap = GaPSE.readxyall(name, comments=true)

        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(ss_trap, calc_ss_trap)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_sums_trap, calc_sums_trap)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[1], calc_xis_trap[1])]) # auto_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[2], calc_xis_trap[2])]) # auto_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[3], calc_xis_trap[3])]) # auto_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[4], calc_xis_trap[4])]) # auto_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[5], calc_xis_trap[5])]) # auto_integrated
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[6], calc_xis_trap[6])]) # newton_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[7], calc_xis_trap[7])]) # doppler_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[8], calc_xis_trap[8])]) # newton_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[9], calc_xis_trap[9])]) # lensing_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[10], calc_xis_trap[10])]) # newton_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[11], calc_xis_trap[11])]) # localgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[12], calc_xis_trap[12])]) # newton_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[13], calc_xis_trap[13])]) # integratedgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[14], calc_xis_trap[14])]) # lensing_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[15], calc_xis_trap[15])]) # doppler_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[16], calc_xis_trap[16])]) # doppler_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[17], calc_xis_trap[17])]) # localgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[18], calc_xis_trap[18])]) # doppler_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[19], calc_xis_trap[19])]) # integratedgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[20], calc_xis_trap[20])]) # lensing_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[21], calc_xis_trap[21])]) # localgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[22], calc_xis_trap[22])]) # lensing_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[23], calc_xis_trap[23])]) # integratedgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[24], calc_xis_trap[24])]) # localgp_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[25], calc_xis_trap[25])]) # integratedgp_localgp

        rm(name)

        #### lob ####
        ss_lob, res_sums_lob, res_xis_lob = GaPSE.readxyall(
            "datatest/GNC_SumXiMultipoles/xis_GNC_L$L" * "_lob_noF_noobsvel_specific_ss.txt", comments=true)

        name = "calc_xis_GNC_L$L" * "_lob_noF_noobsvel_specific_ss.txt"
        isfile(name) && rm(name)

        GaPSE.print_map_sum_ξ_GNC_multipole(COSMO, name, ss_lob;
            single=true, L=L, alg=:lobatto, kwargs...)
        calc_ss_lob, calc_sums_lob, calc_xis_lob = GaPSE.readxyall(name, comments=true)

        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(ss_lob, calc_ss_lob)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_sums_lob, calc_sums_lob)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[1], calc_xis_lob[1])]) # auto_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[2], calc_xis_lob[2])]) # auto_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[3], calc_xis_lob[3])]) # auto_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[4], calc_xis_lob[4])]) # auto_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[5], calc_xis_lob[5])]) # auto_integrated
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[6], calc_xis_lob[6])]) # newton_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[7], calc_xis_lob[7])]) # doppler_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[8], calc_xis_lob[8])]) # newton_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[9], calc_xis_lob[9])]) # lensing_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[10], calc_xis_lob[10])]) # newton_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[11], calc_xis_lob[11])]) # localgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[12], calc_xis_lob[12])]) # newton_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[13], calc_xis_lob[13])]) # integratedgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[14], calc_xis_lob[14])]) # lensing_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[15], calc_xis_lob[15])]) # doppler_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[16], calc_xis_lob[16])]) # doppler_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[17], calc_xis_lob[17])]) # localgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[18], calc_xis_lob[18])]) # doppler_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[19], calc_xis_lob[19])]) # integratedgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[20], calc_xis_lob[20])]) # lensing_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[21], calc_xis_lob[21])]) # localgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[22], calc_xis_lob[22])]) # lensing_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[23], calc_xis_lob[23])]) # integratedgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[24], calc_xis_lob[24])]) # localgp_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[25], calc_xis_lob[25])]) # integratedgp_localgp

        rm(name)

    end
end

println("\n 2 / 8 ...")

@testset "test print_map_sum_ξ_GNC_multipole no_window with obs, first method" begin
    RTOL = 1.2e-2

    kwargs = Dict(
        :use_windows => false, :obs => :yes,
        :enhancer => 1e8, :s1 => nothing,
        :N_trap => 30, :N_lob => 30,
        :atol_quad => 0.0, :rtol_quad => 1e-2,
        :N_χs => 40, :N_χs_2 => 20,
        :pr => false,
    )

    @testset "test_zeros" begin
        @test_throws AssertionError GaPSE.GaPSE.print_map_sum_ξ_GNC_multipole(COSMO, "NONEXISTINGDIR/file.txt", [1.0, 2.0, 3.0])
    end

    @testset "L = 0 no_window with obs" begin
        L = 0

        ##### quad ####
        ss_quad, res_sums_quad, res_xis_quad = GaPSE.readxyall(
            "datatest/GNC_SumXiMultipoles/xis_GNC_L$L" * "_quad_noF_withobs_specific_ss.txt", comments=true)

        name = "calc_xis_GNC_L$L" * "_quad_noF_withobs_specific_ss.txt"
        isfile(name) && rm(name)

        GaPSE.print_map_sum_ξ_GNC_multipole(COSMO, name, ss_quad;
            single=true, L=L, alg=:quad, kwargs...)
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

        #### trap ####
        ss_trap, res_sums_trap, res_xis_trap = GaPSE.readxyall(
            "datatest/GNC_SumXiMultipoles/xis_GNC_L$L" * "_trap_noF_withobs_specific_ss.txt", comments=true)

        name = "calc_xis_GNC_L$L" * "_trap_noF_withobs_specific_ss.txt"
        isfile(name) && rm(name)

        GaPSE.print_map_sum_ξ_GNC_multipole(COSMO, name, ss_trap;
            single=true, L=L, alg=:trap, kwargs...)
        calc_ss_trap, calc_sums_trap, calc_xis_trap = GaPSE.readxyall(name, comments=true)

        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(ss_trap, calc_ss_trap)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_sums_trap, calc_sums_trap)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[1], calc_xis_trap[1])]) # auto_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[2], calc_xis_trap[2])]) # auto_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[3], calc_xis_trap[3])]) # auto_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[4], calc_xis_trap[4])]) # auto_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[5], calc_xis_trap[5])]) # auto_integrated
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[6], calc_xis_trap[6])]) # newton_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[7], calc_xis_trap[7])]) # doppler_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[8], calc_xis_trap[8])]) # newton_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[9], calc_xis_trap[9])]) # lensing_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[10], calc_xis_trap[10])]) # newton_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[11], calc_xis_trap[11])]) # localgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[12], calc_xis_trap[12])]) # newton_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[13], calc_xis_trap[13])]) # integratedgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[14], calc_xis_trap[14])]) # lensing_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[15], calc_xis_trap[15])]) # doppler_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[16], calc_xis_trap[16])]) # doppler_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[17], calc_xis_trap[17])]) # localgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[18], calc_xis_trap[18])]) # doppler_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[19], calc_xis_trap[19])]) # integratedgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[20], calc_xis_trap[20])]) # lensing_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[21], calc_xis_trap[21])]) # localgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[22], calc_xis_trap[22])]) # lensing_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[23], calc_xis_trap[23])]) # integratedgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[24], calc_xis_trap[24])]) # localgp_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[25], calc_xis_trap[25])]) # integratedgp_localgp

        rm(name)

        #### lob ####
        ss_lob, res_sums_lob, res_xis_lob = GaPSE.readxyall(
            "datatest/GNC_SumXiMultipoles/xis_GNC_L$L" * "_lob_noF_withobs_specific_ss.txt", comments=true)

        name = "calc_xis_GNC_L$L" * "_lob_noF_withobs_specific_ss.txt"
        isfile(name) && rm(name)

        GaPSE.print_map_sum_ξ_GNC_multipole(COSMO, name, ss_lob;
            single=true, L=L, alg=:lobatto, kwargs...)
        calc_ss_lob, calc_sums_lob, calc_xis_lob = GaPSE.readxyall(name, comments=true)

        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(ss_lob, calc_ss_lob)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_sums_lob, calc_sums_lob)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[1], calc_xis_lob[1])]) # auto_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[2], calc_xis_lob[2])]) # auto_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[3], calc_xis_lob[3])]) # auto_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[4], calc_xis_lob[4])]) # auto_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[5], calc_xis_lob[5])]) # auto_integrated
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[6], calc_xis_lob[6])]) # newton_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[7], calc_xis_lob[7])]) # doppler_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[8], calc_xis_lob[8])]) # newton_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[9], calc_xis_lob[9])]) # lensing_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[10], calc_xis_lob[10])]) # newton_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[11], calc_xis_lob[11])]) # localgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[12], calc_xis_lob[12])]) # newton_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[13], calc_xis_lob[13])]) # integratedgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[14], calc_xis_lob[14])]) # lensing_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[15], calc_xis_lob[15])]) # doppler_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[16], calc_xis_lob[16])]) # doppler_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[17], calc_xis_lob[17])]) # localgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[18], calc_xis_lob[18])]) # doppler_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[19], calc_xis_lob[19])]) # integratedgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[20], calc_xis_lob[20])]) # lensing_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[21], calc_xis_lob[21])]) # localgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[22], calc_xis_lob[22])]) # lensing_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[23], calc_xis_lob[23])]) # integratedgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[24], calc_xis_lob[24])]) # localgp_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[25], calc_xis_lob[25])]) # integratedgp_localgp

        rm(name)

    end

    @testset "L = 1 no_window with obs" begin
        L = 1

        #### trap ####
        ss_trap, res_sums_trap, res_xis_trap = GaPSE.readxyall(
            "datatest/GNC_SumXiMultipoles/xis_GNC_L$L" * "_trap_noF_withobs_specific_ss.txt", comments=true)

        name = "calc_xis_GNC_L$L" * "_trap_noF_withobs_specific_ss.txt"
        isfile(name) && rm(name)

        GaPSE.print_map_sum_ξ_GNC_multipole(COSMO, name, ss_trap;
            single=true, L=L, alg=:trap, kwargs...)
        calc_ss_trap, calc_sums_trap, calc_xis_trap = GaPSE.readxyall(name, comments=true)

        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(ss_trap, calc_ss_trap)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_sums_trap, calc_sums_trap)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[1], calc_xis_trap[1])]) # auto_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[2], calc_xis_trap[2])]) # auto_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[3], calc_xis_trap[3])]) # auto_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[4], calc_xis_trap[4])]) # auto_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[5], calc_xis_trap[5])]) # auto_integrated
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[6], calc_xis_trap[6])]) # newton_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[7], calc_xis_trap[7])]) # doppler_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[8], calc_xis_trap[8])]) # newton_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[9], calc_xis_trap[9])]) # lensing_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[10], calc_xis_trap[10])]) # newton_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[11], calc_xis_trap[11])]) # localgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[12], calc_xis_trap[12])]) # newton_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[13], calc_xis_trap[13])]) # integratedgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[14], calc_xis_trap[14])]) # lensing_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[15], calc_xis_trap[15])]) # doppler_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[16], calc_xis_trap[16])]) # doppler_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[17], calc_xis_trap[17])]) # localgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[18], calc_xis_trap[18])]) # doppler_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[19], calc_xis_trap[19])]) # integratedgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[20], calc_xis_trap[20])]) # lensing_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[21], calc_xis_trap[21])]) # localgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[22], calc_xis_trap[22])]) # lensing_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[23], calc_xis_trap[23])]) # integratedgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[24], calc_xis_trap[24])]) # localgp_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[25], calc_xis_trap[25])]) # integratedgp_localgp

        rm(name)

        #### lob ####
        ss_lob, res_sums_lob, res_xis_lob = GaPSE.readxyall(
            "datatest/GNC_SumXiMultipoles/xis_GNC_L$L" * "_lob_noF_withobs_specific_ss.txt", comments=true)

        name = "calc_xis_GNC_L$L" * "_lob_noF_withobs_specific_ss.txt"
        isfile(name) && rm(name)

        GaPSE.print_map_sum_ξ_GNC_multipole(COSMO, name, ss_lob;
            single=true, L=L, alg=:lobatto, kwargs...)
        calc_ss_lob, calc_sums_lob, calc_xis_lob = GaPSE.readxyall(name, comments=true)

        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(ss_lob, calc_ss_lob)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_sums_lob, calc_sums_lob)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[1], calc_xis_lob[1])]) # auto_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[2], calc_xis_lob[2])]) # auto_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[3], calc_xis_lob[3])]) # auto_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[4], calc_xis_lob[4])]) # auto_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[5], calc_xis_lob[5])]) # auto_integrated
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[6], calc_xis_lob[6])]) # newton_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[7], calc_xis_lob[7])]) # doppler_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[8], calc_xis_lob[8])]) # newton_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[9], calc_xis_lob[9])]) # lensing_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[10], calc_xis_lob[10])]) # newton_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[11], calc_xis_lob[11])]) # localgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[12], calc_xis_lob[12])]) # newton_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[13], calc_xis_lob[13])]) # integratedgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[14], calc_xis_lob[14])]) # lensing_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[15], calc_xis_lob[15])]) # doppler_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[16], calc_xis_lob[16])]) # doppler_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[17], calc_xis_lob[17])]) # localgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[18], calc_xis_lob[18])]) # doppler_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[19], calc_xis_lob[19])]) # integratedgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[20], calc_xis_lob[20])]) # lensing_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[21], calc_xis_lob[21])]) # localgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[22], calc_xis_lob[22])]) # lensing_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[23], calc_xis_lob[23])]) # integratedgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[24], calc_xis_lob[24])]) # localgp_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[25], calc_xis_lob[25])]) # integratedgp_localgp

        rm(name)

    end

    @testset "L = 2 no_window with obs" begin
        L = 2

        ##### quad ####
        ss_quad, res_sums_quad, res_xis_quad = GaPSE.readxyall(
            "datatest/GNC_SumXiMultipoles/xis_GNC_L$L" * "_quad_noF_withobs_specific_ss.txt", comments=true)

        name = "calc_xis_GNC_L$L" * "_quad_noF_withobs_specific_ss.txt"
        isfile(name) && rm(name)

        GaPSE.print_map_sum_ξ_GNC_multipole(COSMO, name, ss_quad;
            single=true, L=L, alg=:quad, kwargs...)
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

        #### trap ####
        ss_trap, res_sums_trap, res_xis_trap = GaPSE.readxyall(
            "datatest/GNC_SumXiMultipoles/xis_GNC_L$L" * "_trap_noF_withobs_specific_ss.txt", comments=true)

        name = "calc_xis_GNC_L$L" * "_trap_noF_withobs_specific_ss.txt"
        isfile(name) && rm(name)

        GaPSE.print_map_sum_ξ_GNC_multipole(COSMO, name, ss_trap;
            single=true, L=L, alg=:trap, kwargs...)
        calc_ss_trap, calc_sums_trap, calc_xis_trap = GaPSE.readxyall(name, comments=true)

        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(ss_trap, calc_ss_trap)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_sums_trap, calc_sums_trap)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[1], calc_xis_trap[1])]) # auto_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[2], calc_xis_trap[2])]) # auto_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[3], calc_xis_trap[3])]) # auto_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[4], calc_xis_trap[4])]) # auto_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[5], calc_xis_trap[5])]) # auto_integrated
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[6], calc_xis_trap[6])]) # newton_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[7], calc_xis_trap[7])]) # doppler_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[8], calc_xis_trap[8])]) # newton_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[9], calc_xis_trap[9])]) # lensing_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[10], calc_xis_trap[10])]) # newton_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[11], calc_xis_trap[11])]) # localgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[12], calc_xis_trap[12])]) # newton_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[13], calc_xis_trap[13])]) # integratedgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[14], calc_xis_trap[14])]) # lensing_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[15], calc_xis_trap[15])]) # doppler_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[16], calc_xis_trap[16])]) # doppler_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[17], calc_xis_trap[17])]) # localgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[18], calc_xis_trap[18])]) # doppler_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[19], calc_xis_trap[19])]) # integratedgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[20], calc_xis_trap[20])]) # lensing_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[21], calc_xis_trap[21])]) # localgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[22], calc_xis_trap[22])]) # lensing_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[23], calc_xis_trap[23])]) # integratedgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[24], calc_xis_trap[24])]) # localgp_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[25], calc_xis_trap[25])]) # integratedgp_localgp

        rm(name)

        #### lob ####
        ss_lob, res_sums_lob, res_xis_lob = GaPSE.readxyall(
            "datatest/GNC_SumXiMultipoles/xis_GNC_L$L" * "_lob_noF_withobs_specific_ss.txt", comments=true)

        name = "calc_xis_GNC_L$L" * "_lob_noF_withobs_specific_ss.txt"
        isfile(name) && rm(name)

        GaPSE.print_map_sum_ξ_GNC_multipole(COSMO, name, ss_lob;
            single=true, L=L, alg=:lobatto, kwargs...)
        calc_ss_lob, calc_sums_lob, calc_xis_lob = GaPSE.readxyall(name, comments=true)

        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(ss_lob, calc_ss_lob)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_sums_lob, calc_sums_lob)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[1], calc_xis_lob[1])]) # auto_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[2], calc_xis_lob[2])]) # auto_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[3], calc_xis_lob[3])]) # auto_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[4], calc_xis_lob[4])]) # auto_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[5], calc_xis_lob[5])]) # auto_integrated
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[6], calc_xis_lob[6])]) # newton_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[7], calc_xis_lob[7])]) # doppler_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[8], calc_xis_lob[8])]) # newton_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[9], calc_xis_lob[9])]) # lensing_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[10], calc_xis_lob[10])]) # newton_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[11], calc_xis_lob[11])]) # localgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[12], calc_xis_lob[12])]) # newton_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[13], calc_xis_lob[13])]) # integratedgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[14], calc_xis_lob[14])]) # lensing_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[15], calc_xis_lob[15])]) # doppler_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[16], calc_xis_lob[16])]) # doppler_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[17], calc_xis_lob[17])]) # localgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[18], calc_xis_lob[18])]) # doppler_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[19], calc_xis_lob[19])]) # integratedgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[20], calc_xis_lob[20])]) # lensing_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[21], calc_xis_lob[21])]) # localgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[22], calc_xis_lob[22])]) # lensing_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[23], calc_xis_lob[23])]) # integratedgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[24], calc_xis_lob[24])]) # localgp_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[25], calc_xis_lob[25])]) # integratedgp_localgp

        rm(name)

    end

    @testset "L = 3 no_window with obs" begin
        L = 3

        #### trap ####
        ss_trap, res_sums_trap, res_xis_trap = GaPSE.readxyall(
            "datatest/GNC_SumXiMultipoles/xis_GNC_L$L" * "_trap_noF_withobs_specific_ss.txt", comments=true)

        name = "calc_xis_GNC_L$L" * "_trap_noF_withobs_specific_ss.txt"
        isfile(name) && rm(name)

        GaPSE.print_map_sum_ξ_GNC_multipole(COSMO, name, ss_trap;
            single=true, L=L, alg=:trap, kwargs...)
        calc_ss_trap, calc_sums_trap, calc_xis_trap = GaPSE.readxyall(name, comments=true)

        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(ss_trap, calc_ss_trap)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_sums_trap, calc_sums_trap)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[1], calc_xis_trap[1])]) # auto_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[2], calc_xis_trap[2])]) # auto_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[3], calc_xis_trap[3])]) # auto_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[4], calc_xis_trap[4])]) # auto_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[5], calc_xis_trap[5])]) # auto_integrated
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[6], calc_xis_trap[6])]) # newton_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[7], calc_xis_trap[7])]) # doppler_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[8], calc_xis_trap[8])]) # newton_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[9], calc_xis_trap[9])]) # lensing_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[10], calc_xis_trap[10])]) # newton_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[11], calc_xis_trap[11])]) # localgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[12], calc_xis_trap[12])]) # newton_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[13], calc_xis_trap[13])]) # integratedgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[14], calc_xis_trap[14])]) # lensing_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[15], calc_xis_trap[15])]) # doppler_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[16], calc_xis_trap[16])]) # doppler_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[17], calc_xis_trap[17])]) # localgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[18], calc_xis_trap[18])]) # doppler_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[19], calc_xis_trap[19])]) # integratedgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[20], calc_xis_trap[20])]) # lensing_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[21], calc_xis_trap[21])]) # localgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[22], calc_xis_trap[22])]) # lensing_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[23], calc_xis_trap[23])]) # integratedgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[24], calc_xis_trap[24])]) # localgp_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[25], calc_xis_trap[25])]) # integratedgp_localgp

        rm(name)

        #### lob ####
        ss_lob, res_sums_lob, res_xis_lob = GaPSE.readxyall(
            "datatest/GNC_SumXiMultipoles/xis_GNC_L$L" * "_lob_noF_withobs_specific_ss.txt", comments=true)

        name = "calc_xis_GNC_L$L" * "_lob_noF_withobs_specific_ss.txt"
        isfile(name) && rm(name)

        GaPSE.print_map_sum_ξ_GNC_multipole(COSMO, name, ss_lob;
            single=true, L=L, alg=:lobatto, kwargs...)
        calc_ss_lob, calc_sums_lob, calc_xis_lob = GaPSE.readxyall(name, comments=true)

        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(ss_lob, calc_ss_lob)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_sums_lob, calc_sums_lob)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[1], calc_xis_lob[1])]) # auto_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[2], calc_xis_lob[2])]) # auto_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[3], calc_xis_lob[3])]) # auto_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[4], calc_xis_lob[4])]) # auto_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[5], calc_xis_lob[5])]) # auto_integrated
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[6], calc_xis_lob[6])]) # newton_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[7], calc_xis_lob[7])]) # doppler_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[8], calc_xis_lob[8])]) # newton_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[9], calc_xis_lob[9])]) # lensing_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[10], calc_xis_lob[10])]) # newton_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[11], calc_xis_lob[11])]) # localgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[12], calc_xis_lob[12])]) # newton_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[13], calc_xis_lob[13])]) # integratedgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[14], calc_xis_lob[14])]) # lensing_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[15], calc_xis_lob[15])]) # doppler_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[16], calc_xis_lob[16])]) # doppler_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[17], calc_xis_lob[17])]) # localgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[18], calc_xis_lob[18])]) # doppler_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[19], calc_xis_lob[19])]) # integratedgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[20], calc_xis_lob[20])]) # lensing_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[21], calc_xis_lob[21])]) # localgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[22], calc_xis_lob[22])]) # lensing_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[23], calc_xis_lob[23])]) # integratedgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[24], calc_xis_lob[24])]) # localgp_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[25], calc_xis_lob[25])]) # integratedgp_localgp

        rm(name)

    end

    @testset "L = 4 no_window with obs" begin
        L = 4

        #### trap ####
        ss_trap, res_sums_trap, res_xis_trap = GaPSE.readxyall(
            "datatest/GNC_SumXiMultipoles/xis_GNC_L$L" * "_trap_noF_withobs_specific_ss.txt", comments=true)

        name = "calc_xis_GNC_L$L" * "_trap_noF_withobs_specific_ss.txt"
        isfile(name) && rm(name)

        GaPSE.print_map_sum_ξ_GNC_multipole(COSMO, name, ss_trap;
            single=true, L=L, alg=:trap, kwargs...)
        calc_ss_trap, calc_sums_trap, calc_xis_trap = GaPSE.readxyall(name, comments=true)

        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(ss_trap, calc_ss_trap)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_sums_trap, calc_sums_trap)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[1], calc_xis_trap[1])]) # auto_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[2], calc_xis_trap[2])]) # auto_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[3], calc_xis_trap[3])]) # auto_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[4], calc_xis_trap[4])]) # auto_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[5], calc_xis_trap[5])]) # auto_integrated
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[6], calc_xis_trap[6])]) # newton_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[7], calc_xis_trap[7])]) # doppler_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[8], calc_xis_trap[8])]) # newton_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[9], calc_xis_trap[9])]) # lensing_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[10], calc_xis_trap[10])]) # newton_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[11], calc_xis_trap[11])]) # localgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[12], calc_xis_trap[12])]) # newton_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[13], calc_xis_trap[13])]) # integratedgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[14], calc_xis_trap[14])]) # lensing_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[15], calc_xis_trap[15])]) # doppler_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[16], calc_xis_trap[16])]) # doppler_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[17], calc_xis_trap[17])]) # localgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[18], calc_xis_trap[18])]) # doppler_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[19], calc_xis_trap[19])]) # integratedgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[20], calc_xis_trap[20])]) # lensing_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[21], calc_xis_trap[21])]) # localgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[22], calc_xis_trap[22])]) # lensing_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[23], calc_xis_trap[23])]) # integratedgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[24], calc_xis_trap[24])]) # localgp_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[25], calc_xis_trap[25])]) # integratedgp_localgp

        rm(name)

        #### lob ####
        ss_lob, res_sums_lob, res_xis_lob = GaPSE.readxyall(
            "datatest/GNC_SumXiMultipoles/xis_GNC_L$L" * "_lob_noF_withobs_specific_ss.txt", comments=true)

        name = "calc_xis_GNC_L$L" * "_lob_noF_withobs_specific_ss.txt"
        isfile(name) && rm(name)

        GaPSE.print_map_sum_ξ_GNC_multipole(COSMO, name, ss_lob;
            single=true, L=L, alg=:lobatto, kwargs...)
        calc_ss_lob, calc_sums_lob, calc_xis_lob = GaPSE.readxyall(name, comments=true)

        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(ss_lob, calc_ss_lob)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_sums_lob, calc_sums_lob)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[1], calc_xis_lob[1])]) # auto_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[2], calc_xis_lob[2])]) # auto_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[3], calc_xis_lob[3])]) # auto_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[4], calc_xis_lob[4])]) # auto_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[5], calc_xis_lob[5])]) # auto_integrated
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[6], calc_xis_lob[6])]) # newton_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[7], calc_xis_lob[7])]) # doppler_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[8], calc_xis_lob[8])]) # newton_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[9], calc_xis_lob[9])]) # lensing_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[10], calc_xis_lob[10])]) # newton_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[11], calc_xis_lob[11])]) # localgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[12], calc_xis_lob[12])]) # newton_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[13], calc_xis_lob[13])]) # integratedgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[14], calc_xis_lob[14])]) # lensing_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[15], calc_xis_lob[15])]) # doppler_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[16], calc_xis_lob[16])]) # doppler_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[17], calc_xis_lob[17])]) # localgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[18], calc_xis_lob[18])]) # doppler_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[19], calc_xis_lob[19])]) # integratedgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[20], calc_xis_lob[20])]) # lensing_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[21], calc_xis_lob[21])]) # localgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[22], calc_xis_lob[22])]) # lensing_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[23], calc_xis_lob[23])]) # integratedgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[24], calc_xis_lob[24])]) # localgp_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[25], calc_xis_lob[25])]) # integratedgp_localgp

        rm(name)

    end
end

println("\n 3 / 8 ...")

@testset "test print_map_sum_ξ_GNC_multipole with_window no obs, first method" begin
    RTOL = 1.2e-2

    kwargs = Dict(
        :use_windows => true, :obs => :no,
        :enhancer => 1e8, :s1 => nothing,
        :N_trap => 30, :N_lob => 30,
        :atol_quad => 0.0, :rtol_quad => 1e-2,
        :N_χs => 40, :N_χs_2 => 20,
        :pr => false,
    )

    @testset "test_zeros" begin
        @test_throws AssertionError GaPSE.GaPSE.print_map_sum_ξ_GNC_multipole(COSMO, "NONEXISTINGDIR/file.txt", [1.0, 2.0, 3.0])
    end

    @testset "L = 0 with_window no obs" begin
        L = 0

        ##### quad ####
        ss_quad, res_sums_quad, res_xis_quad = GaPSE.readxyall(
            "datatest/GNC_SumXiMultipoles/xis_GNC_L$L" * "_quad_withF_noobs_specific_ss.txt", comments=true)

        name = "calc_xis_GNC_L$L" * "_quad_withF_noobs_specific_ss.txt"
        isfile(name) && rm(name)

        GaPSE.print_map_sum_ξ_GNC_multipole(COSMO, name, ss_quad;
            single=true, L=L, alg=:quad, kwargs...)
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

        #### trap ####
        ss_trap, res_sums_trap, res_xis_trap = GaPSE.readxyall(
            "datatest/GNC_SumXiMultipoles/xis_GNC_L$L" * "_trap_withF_noobs_specific_ss.txt", comments=true)

        name = "calc_xis_GNC_L$L" * "_trap_withF_noobs_specific_ss.txt"
        isfile(name) && rm(name)

        GaPSE.print_map_sum_ξ_GNC_multipole(COSMO, name, ss_trap;
            single=true, L=L, alg=:trap, kwargs...)
        calc_ss_trap, calc_sums_trap, calc_xis_trap = GaPSE.readxyall(name, comments=true)

        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(ss_trap, calc_ss_trap)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_sums_trap, calc_sums_trap)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[1], calc_xis_trap[1])]) # auto_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[2], calc_xis_trap[2])]) # auto_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[3], calc_xis_trap[3])]) # auto_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[4], calc_xis_trap[4])]) # auto_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[5], calc_xis_trap[5])]) # auto_integrated
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[6], calc_xis_trap[6])]) # newton_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[7], calc_xis_trap[7])]) # doppler_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[8], calc_xis_trap[8])]) # newton_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[9], calc_xis_trap[9])]) # lensing_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[10], calc_xis_trap[10])]) # newton_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[11], calc_xis_trap[11])]) # localgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[12], calc_xis_trap[12])]) # newton_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[13], calc_xis_trap[13])]) # integratedgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[14], calc_xis_trap[14])]) # lensing_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[15], calc_xis_trap[15])]) # doppler_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[16], calc_xis_trap[16])]) # doppler_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[17], calc_xis_trap[17])]) # localgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[18], calc_xis_trap[18])]) # doppler_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[19], calc_xis_trap[19])]) # integratedgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[20], calc_xis_trap[20])]) # lensing_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[21], calc_xis_trap[21])]) # localgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[22], calc_xis_trap[22])]) # lensing_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[23], calc_xis_trap[23])]) # integratedgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[24], calc_xis_trap[24])]) # localgp_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[25], calc_xis_trap[25])]) # integratedgp_localgp

        rm(name)

        #### lob ####
        ss_lob, res_sums_lob, res_xis_lob = GaPSE.readxyall(
            "datatest/GNC_SumXiMultipoles/xis_GNC_L$L" * "_lob_withF_noobs_specific_ss.txt", comments=true)

        name = "calc_xis_GNC_L$L" * "_lob_withF_noobs_specific_ss.txt"
        isfile(name) && rm(name)

        GaPSE.print_map_sum_ξ_GNC_multipole(COSMO, name, ss_lob;
            single=true, L=L, alg=:lobatto, kwargs...)
        calc_ss_lob, calc_sums_lob, calc_xis_lob = GaPSE.readxyall(name, comments=true)

        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(ss_lob, calc_ss_lob)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_sums_lob, calc_sums_lob)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[1], calc_xis_lob[1])]) # auto_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[2], calc_xis_lob[2])]) # auto_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[3], calc_xis_lob[3])]) # auto_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[4], calc_xis_lob[4])]) # auto_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[5], calc_xis_lob[5])]) # auto_integrated
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[6], calc_xis_lob[6])]) # newton_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[7], calc_xis_lob[7])]) # doppler_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[8], calc_xis_lob[8])]) # newton_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[9], calc_xis_lob[9])]) # lensing_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[10], calc_xis_lob[10])]) # newton_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[11], calc_xis_lob[11])]) # localgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[12], calc_xis_lob[12])]) # newton_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[13], calc_xis_lob[13])]) # integratedgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[14], calc_xis_lob[14])]) # lensing_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[15], calc_xis_lob[15])]) # doppler_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[16], calc_xis_lob[16])]) # doppler_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[17], calc_xis_lob[17])]) # localgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[18], calc_xis_lob[18])]) # doppler_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[19], calc_xis_lob[19])]) # integratedgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[20], calc_xis_lob[20])]) # lensing_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[21], calc_xis_lob[21])]) # localgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[22], calc_xis_lob[22])]) # lensing_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[23], calc_xis_lob[23])]) # integratedgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[24], calc_xis_lob[24])]) # localgp_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[25], calc_xis_lob[25])]) # integratedgp_localgp

        rm(name)

    end

    @testset "L = 1 with_window no obs" begin
        L = 1

        #### trap ####
        ss_trap, res_sums_trap, res_xis_trap = GaPSE.readxyall(
            "datatest/GNC_SumXiMultipoles/xis_GNC_L$L" * "_trap_withF_noobs_specific_ss.txt", comments=true)

        name = "calc_xis_GNC_L$L" * "_trap_withF_noobs_specific_ss.txt"
        isfile(name) && rm(name)

        GaPSE.print_map_sum_ξ_GNC_multipole(COSMO, name, ss_trap;
            single=true, L=L, alg=:trap, kwargs...)
        calc_ss_trap, calc_sums_trap, calc_xis_trap = GaPSE.readxyall(name, comments=true)

        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(ss_trap, calc_ss_trap)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_sums_trap, calc_sums_trap)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[1], calc_xis_trap[1])]) # auto_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[2], calc_xis_trap[2])]) # auto_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[3], calc_xis_trap[3])]) # auto_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[4], calc_xis_trap[4])]) # auto_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[5], calc_xis_trap[5])]) # auto_integrated
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[6], calc_xis_trap[6])]) # newton_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[7], calc_xis_trap[7])]) # doppler_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[8], calc_xis_trap[8])]) # newton_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[9], calc_xis_trap[9])]) # lensing_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[10], calc_xis_trap[10])]) # newton_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[11], calc_xis_trap[11])]) # localgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[12], calc_xis_trap[12])]) # newton_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[13], calc_xis_trap[13])]) # integratedgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[14], calc_xis_trap[14])]) # lensing_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[15], calc_xis_trap[15])]) # doppler_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[16], calc_xis_trap[16])]) # doppler_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[17], calc_xis_trap[17])]) # localgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[18], calc_xis_trap[18])]) # doppler_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[19], calc_xis_trap[19])]) # integratedgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[20], calc_xis_trap[20])]) # lensing_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[21], calc_xis_trap[21])]) # localgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[22], calc_xis_trap[22])]) # lensing_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[23], calc_xis_trap[23])]) # integratedgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[24], calc_xis_trap[24])]) # localgp_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[25], calc_xis_trap[25])]) # integratedgp_localgp

        rm(name)

        #### lob ####
        ss_lob, res_sums_lob, res_xis_lob = GaPSE.readxyall(
            "datatest/GNC_SumXiMultipoles/xis_GNC_L$L" * "_lob_withF_noobs_specific_ss.txt", comments=true)

        name = "calc_xis_GNC_L$L" * "_lob_withF_noobs_specific_ss.txt"
        isfile(name) && rm(name)

        GaPSE.print_map_sum_ξ_GNC_multipole(COSMO, name, ss_lob;
            single=true, L=L, alg=:lobatto, kwargs...)
        calc_ss_lob, calc_sums_lob, calc_xis_lob = GaPSE.readxyall(name, comments=true)

        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(ss_lob, calc_ss_lob)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_sums_lob, calc_sums_lob)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[1], calc_xis_lob[1])]) # auto_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[2], calc_xis_lob[2])]) # auto_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[3], calc_xis_lob[3])]) # auto_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[4], calc_xis_lob[4])]) # auto_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[5], calc_xis_lob[5])]) # auto_integrated
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[6], calc_xis_lob[6])]) # newton_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[7], calc_xis_lob[7])]) # doppler_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[8], calc_xis_lob[8])]) # newton_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[9], calc_xis_lob[9])]) # lensing_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[10], calc_xis_lob[10])]) # newton_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[11], calc_xis_lob[11])]) # localgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[12], calc_xis_lob[12])]) # newton_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[13], calc_xis_lob[13])]) # integratedgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[14], calc_xis_lob[14])]) # lensing_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[15], calc_xis_lob[15])]) # doppler_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[16], calc_xis_lob[16])]) # doppler_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[17], calc_xis_lob[17])]) # localgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[18], calc_xis_lob[18])]) # doppler_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[19], calc_xis_lob[19])]) # integratedgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[20], calc_xis_lob[20])]) # lensing_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[21], calc_xis_lob[21])]) # localgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[22], calc_xis_lob[22])]) # lensing_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[23], calc_xis_lob[23])]) # integratedgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[24], calc_xis_lob[24])]) # localgp_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[25], calc_xis_lob[25])]) # integratedgp_localgp

        rm(name)

    end

    @testset "L = 2 with_window no obs" begin
        L = 2

        ##### quad ####
        ss_quad, res_sums_quad, res_xis_quad = GaPSE.readxyall(
            "datatest/GNC_SumXiMultipoles/xis_GNC_L$L" * "_quad_withF_noobs_specific_ss.txt", comments=true)

        name = "calc_xis_GNC_L$L" * "_quad_withF_noobs_specific_ss.txt"
        isfile(name) && rm(name)

        GaPSE.print_map_sum_ξ_GNC_multipole(COSMO, name, ss_quad;
            single=true, L=L, alg=:quad, kwargs...)
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

        #### trap ####
        ss_trap, res_sums_trap, res_xis_trap = GaPSE.readxyall(
            "datatest/GNC_SumXiMultipoles/xis_GNC_L$L" * "_trap_withF_noobs_specific_ss.txt", comments=true)

        name = "calc_xis_GNC_L$L" * "_trap_withF_noobs_specific_ss.txt"
        isfile(name) && rm(name)

        GaPSE.print_map_sum_ξ_GNC_multipole(COSMO, name, ss_trap;
            single=true, L=L, alg=:trap, kwargs...)
        calc_ss_trap, calc_sums_trap, calc_xis_trap = GaPSE.readxyall(name, comments=true)

        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(ss_trap, calc_ss_trap)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_sums_trap, calc_sums_trap)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[1], calc_xis_trap[1])]) # auto_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[2], calc_xis_trap[2])]) # auto_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[3], calc_xis_trap[3])]) # auto_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[4], calc_xis_trap[4])]) # auto_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[5], calc_xis_trap[5])]) # auto_integrated
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[6], calc_xis_trap[6])]) # newton_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[7], calc_xis_trap[7])]) # doppler_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[8], calc_xis_trap[8])]) # newton_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[9], calc_xis_trap[9])]) # lensing_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[10], calc_xis_trap[10])]) # newton_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[11], calc_xis_trap[11])]) # localgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[12], calc_xis_trap[12])]) # newton_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[13], calc_xis_trap[13])]) # integratedgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[14], calc_xis_trap[14])]) # lensing_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[15], calc_xis_trap[15])]) # doppler_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[16], calc_xis_trap[16])]) # doppler_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[17], calc_xis_trap[17])]) # localgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[18], calc_xis_trap[18])]) # doppler_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[19], calc_xis_trap[19])]) # integratedgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[20], calc_xis_trap[20])]) # lensing_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[21], calc_xis_trap[21])]) # localgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[22], calc_xis_trap[22])]) # lensing_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[23], calc_xis_trap[23])]) # integratedgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[24], calc_xis_trap[24])]) # localgp_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[25], calc_xis_trap[25])]) # integratedgp_localgp

        rm(name)

        #### lob ####
        ss_lob, res_sums_lob, res_xis_lob = GaPSE.readxyall(
            "datatest/GNC_SumXiMultipoles/xis_GNC_L$L" * "_lob_withF_noobs_specific_ss.txt", comments=true)

        name = "calc_xis_GNC_L$L" * "_lob_withF_noobs_specific_ss.txt"
        isfile(name) && rm(name)

        GaPSE.print_map_sum_ξ_GNC_multipole(COSMO, name, ss_lob;
            single=true, L=L, alg=:lobatto, kwargs...)
        calc_ss_lob, calc_sums_lob, calc_xis_lob = GaPSE.readxyall(name, comments=true)

        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(ss_lob, calc_ss_lob)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_sums_lob, calc_sums_lob)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[1], calc_xis_lob[1])]) # auto_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[2], calc_xis_lob[2])]) # auto_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[3], calc_xis_lob[3])]) # auto_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[4], calc_xis_lob[4])]) # auto_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[5], calc_xis_lob[5])]) # auto_integrated
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[6], calc_xis_lob[6])]) # newton_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[7], calc_xis_lob[7])]) # doppler_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[8], calc_xis_lob[8])]) # newton_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[9], calc_xis_lob[9])]) # lensing_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[10], calc_xis_lob[10])]) # newton_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[11], calc_xis_lob[11])]) # localgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[12], calc_xis_lob[12])]) # newton_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[13], calc_xis_lob[13])]) # integratedgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[14], calc_xis_lob[14])]) # lensing_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[15], calc_xis_lob[15])]) # doppler_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[16], calc_xis_lob[16])]) # doppler_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[17], calc_xis_lob[17])]) # localgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[18], calc_xis_lob[18])]) # doppler_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[19], calc_xis_lob[19])]) # integratedgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[20], calc_xis_lob[20])]) # lensing_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[21], calc_xis_lob[21])]) # localgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[22], calc_xis_lob[22])]) # lensing_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[23], calc_xis_lob[23])]) # integratedgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[24], calc_xis_lob[24])]) # localgp_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[25], calc_xis_lob[25])]) # integratedgp_localgp

        rm(name)

    end

    @testset "L = 3 with_window no obs" begin
        L = 3

        #### trap ####
        ss_trap, res_sums_trap, res_xis_trap = GaPSE.readxyall(
            "datatest/GNC_SumXiMultipoles/xis_GNC_L$L" * "_trap_withF_noobs_specific_ss.txt", comments=true)

        name = "calc_xis_GNC_L$L" * "_trap_withF_noobs_specific_ss.txt"
        isfile(name) && rm(name)

        GaPSE.print_map_sum_ξ_GNC_multipole(COSMO, name, ss_trap;
            single=true, L=L, alg=:trap, kwargs...)
        calc_ss_trap, calc_sums_trap, calc_xis_trap = GaPSE.readxyall(name, comments=true)

        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(ss_trap, calc_ss_trap)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_sums_trap, calc_sums_trap)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[1], calc_xis_trap[1])]) # auto_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[2], calc_xis_trap[2])]) # auto_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[3], calc_xis_trap[3])]) # auto_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[4], calc_xis_trap[4])]) # auto_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[5], calc_xis_trap[5])]) # auto_integrated
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[6], calc_xis_trap[6])]) # newton_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[7], calc_xis_trap[7])]) # doppler_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[8], calc_xis_trap[8])]) # newton_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[9], calc_xis_trap[9])]) # lensing_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[10], calc_xis_trap[10])]) # newton_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[11], calc_xis_trap[11])]) # localgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[12], calc_xis_trap[12])]) # newton_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[13], calc_xis_trap[13])]) # integratedgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[14], calc_xis_trap[14])]) # lensing_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[15], calc_xis_trap[15])]) # doppler_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[16], calc_xis_trap[16])]) # doppler_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[17], calc_xis_trap[17])]) # localgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[18], calc_xis_trap[18])]) # doppler_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[19], calc_xis_trap[19])]) # integratedgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[20], calc_xis_trap[20])]) # lensing_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[21], calc_xis_trap[21])]) # localgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[22], calc_xis_trap[22])]) # lensing_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[23], calc_xis_trap[23])]) # integratedgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[24], calc_xis_trap[24])]) # localgp_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[25], calc_xis_trap[25])]) # integratedgp_localgp

        rm(name)

        #### lob ####
        ss_lob, res_sums_lob, res_xis_lob = GaPSE.readxyall(
            "datatest/GNC_SumXiMultipoles/xis_GNC_L$L" * "_lob_withF_noobs_specific_ss.txt", comments=true)

        name = "calc_xis_GNC_L$L" * "_lob_withF_noobs_specific_ss.txt"
        isfile(name) && rm(name)

        GaPSE.print_map_sum_ξ_GNC_multipole(COSMO, name, ss_lob;
            single=true, L=L, alg=:lobatto, kwargs...)
        calc_ss_lob, calc_sums_lob, calc_xis_lob = GaPSE.readxyall(name, comments=true)

        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(ss_lob, calc_ss_lob)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_sums_lob, calc_sums_lob)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[1], calc_xis_lob[1])]) # auto_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[2], calc_xis_lob[2])]) # auto_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[3], calc_xis_lob[3])]) # auto_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[4], calc_xis_lob[4])]) # auto_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[5], calc_xis_lob[5])]) # auto_integrated
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[6], calc_xis_lob[6])]) # newton_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[7], calc_xis_lob[7])]) # doppler_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[8], calc_xis_lob[8])]) # newton_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[9], calc_xis_lob[9])]) # lensing_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[10], calc_xis_lob[10])]) # newton_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[11], calc_xis_lob[11])]) # localgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[12], calc_xis_lob[12])]) # newton_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[13], calc_xis_lob[13])]) # integratedgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[14], calc_xis_lob[14])]) # lensing_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[15], calc_xis_lob[15])]) # doppler_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[16], calc_xis_lob[16])]) # doppler_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[17], calc_xis_lob[17])]) # localgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[18], calc_xis_lob[18])]) # doppler_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[19], calc_xis_lob[19])]) # integratedgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[20], calc_xis_lob[20])]) # lensing_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[21], calc_xis_lob[21])]) # localgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[22], calc_xis_lob[22])]) # lensing_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[23], calc_xis_lob[23])]) # integratedgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[24], calc_xis_lob[24])]) # localgp_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[25], calc_xis_lob[25])]) # integratedgp_localgp

        rm(name)

    end

    @testset "L = 4 with_window no obs" begin
        L = 4

        #### trap ####
        ss_trap, res_sums_trap, res_xis_trap = GaPSE.readxyall(
            "datatest/GNC_SumXiMultipoles/xis_GNC_L$L" * "_trap_withF_noobs_specific_ss.txt", comments=true)

        name = "calc_xis_GNC_L$L" * "_trap_withF_noobs_specific_ss.txt"
        isfile(name) && rm(name)

        GaPSE.print_map_sum_ξ_GNC_multipole(COSMO, name, ss_trap;
            single=true, L=L, alg=:trap, kwargs...)
        calc_ss_trap, calc_sums_trap, calc_xis_trap = GaPSE.readxyall(name, comments=true)

        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(ss_trap, calc_ss_trap)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_sums_trap, calc_sums_trap)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[1], calc_xis_trap[1])]) # auto_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[2], calc_xis_trap[2])]) # auto_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[3], calc_xis_trap[3])]) # auto_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[4], calc_xis_trap[4])]) # auto_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[5], calc_xis_trap[5])]) # auto_integrated
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[6], calc_xis_trap[6])]) # newton_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[7], calc_xis_trap[7])]) # doppler_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[8], calc_xis_trap[8])]) # newton_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[9], calc_xis_trap[9])]) # lensing_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[10], calc_xis_trap[10])]) # newton_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[11], calc_xis_trap[11])]) # localgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[12], calc_xis_trap[12])]) # newton_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[13], calc_xis_trap[13])]) # integratedgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[14], calc_xis_trap[14])]) # lensing_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[15], calc_xis_trap[15])]) # doppler_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[16], calc_xis_trap[16])]) # doppler_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[17], calc_xis_trap[17])]) # localgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[18], calc_xis_trap[18])]) # doppler_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[19], calc_xis_trap[19])]) # integratedgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[20], calc_xis_trap[20])]) # lensing_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[21], calc_xis_trap[21])]) # localgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[22], calc_xis_trap[22])]) # lensing_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[23], calc_xis_trap[23])]) # integratedgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[24], calc_xis_trap[24])]) # localgp_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[25], calc_xis_trap[25])]) # integratedgp_localgp

        rm(name)

        #### lob ####
        ss_lob, res_sums_lob, res_xis_lob = GaPSE.readxyall(
            "datatest/GNC_SumXiMultipoles/xis_GNC_L$L" * "_lob_withF_noobs_specific_ss.txt", comments=true)

        name = "calc_xis_GNC_L$L" * "_lob_withF_noobs_specific_ss.txt"
        isfile(name) && rm(name)

        GaPSE.print_map_sum_ξ_GNC_multipole(COSMO, name, ss_lob;
            single=true, L=L, alg=:lobatto, kwargs...)
        calc_ss_lob, calc_sums_lob, calc_xis_lob = GaPSE.readxyall(name, comments=true)

        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(ss_lob, calc_ss_lob)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_sums_lob, calc_sums_lob)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[1], calc_xis_lob[1])]) # auto_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[2], calc_xis_lob[2])]) # auto_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[3], calc_xis_lob[3])]) # auto_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[4], calc_xis_lob[4])]) # auto_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[5], calc_xis_lob[5])]) # auto_integrated
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[6], calc_xis_lob[6])]) # newton_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[7], calc_xis_lob[7])]) # doppler_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[8], calc_xis_lob[8])]) # newton_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[9], calc_xis_lob[9])]) # lensing_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[10], calc_xis_lob[10])]) # newton_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[11], calc_xis_lob[11])]) # localgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[12], calc_xis_lob[12])]) # newton_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[13], calc_xis_lob[13])]) # integratedgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[14], calc_xis_lob[14])]) # lensing_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[15], calc_xis_lob[15])]) # doppler_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[16], calc_xis_lob[16])]) # doppler_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[17], calc_xis_lob[17])]) # localgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[18], calc_xis_lob[18])]) # doppler_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[19], calc_xis_lob[19])]) # integratedgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[20], calc_xis_lob[20])]) # lensing_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[21], calc_xis_lob[21])]) # localgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[22], calc_xis_lob[22])]) # lensing_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[23], calc_xis_lob[23])]) # integratedgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[24], calc_xis_lob[24])]) # localgp_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[25], calc_xis_lob[25])]) # integratedgp_localgp

        rm(name)

    end
end

println("\nI am at half of this path...")

##########################################################################################92


@testset "test print_map_sum_ξ_GNC_multipole no_window no obs, second method" begin
    RTOL = 1.2e-2

    kwargs = Dict(
        :use_windows => false, :obs => :no,
        :enhancer => 1e8, :s1 => nothing,
        :N_trap => 30, :N_lob => 30,
        :atol_quad => 0.0, :rtol_quad => 1e-2,
        :N_χs => 40, :N_χs_2 => 20,
        :pr => false,
    )

    @testset "test_zeros" begin
        @test_throws AssertionError GaPSE.GaPSE.print_map_sum_ξ_GNC_multipole(COSMO, "NONEXISTINGDIR/file.txt", [1.0, 2.0, 3.0])
    end

    @testset "L = 0 no_window no obs" begin
        L = 0

        ##### quad ####
        ss_quad, res_sums_quad, res_xis_quad = GaPSE.readxyall(
            "datatest/GNC_SumXiMultipoles/xis_GNC_L$L" * "_quad_noF_noobs_specific_ss.txt", comments=true)

        name = "calc_xis_GNC_L$L" * "_quad_noF_noobs_specific_ss.txt"
        isfile(name) && rm(name)

        GaPSE.print_map_sum_ξ_GNC_multipole(COSMO, name, ss_quad;
            single=false, L=L, alg=:quad, kwargs...)

        calc_ss_quad, calc_sums_quad = GaPSE.readxy(name, comments=true)
        calc_xis_quad = [
            begin
                a_name = "all_GNC_standalones_CFs/xi_GNC_" * effect * "_L$L" * ".txt"
                _, xis = GaPSE.readxy(a_name, comments=true)
                xis
            end for effect in GaPSE.GR_EFFECTS_GNC
        ]

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
        for effect in GaPSE.GR_EFFECTS_GNC
            rm("all_GNC_standalones_CFs/xi_GNC_" * effect * "_L$L" * ".txt")
        end
        rm("all_GNC_standalones_CFs")

        #### trap ####
        ss_trap, res_sums_trap, res_xis_trap = GaPSE.readxyall(
            "datatest/GNC_SumXiMultipoles/xis_GNC_L$L" * "_trap_noF_noobs_specific_ss.txt", comments=true)

        name = "calc_xis_GNC_L$L" * "_trap_noF_noobs_specific_ss.txt"
        isfile(name) && rm(name)

        GaPSE.print_map_sum_ξ_GNC_multipole(COSMO, name, ss_trap;
            single=false, L=L, alg=:trap, kwargs...)
        calc_ss_trap, calc_sums_trap = GaPSE.readxy(name, comments=true)
        calc_xis_trap = [
            begin
                a_name = "all_GNC_standalones_CFs/xi_GNC_" * effect * "_L$L" * ".txt"
                _, xis = GaPSE.readxy(a_name, comments=true)
                xis
            end for effect in GaPSE.GR_EFFECTS_GNC
        ]

        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(ss_trap, calc_ss_trap)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_sums_trap, calc_sums_trap)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[1], calc_xis_trap[1])]) # auto_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[2], calc_xis_trap[2])]) # auto_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[3], calc_xis_trap[3])]) # auto_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[4], calc_xis_trap[4])]) # auto_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[5], calc_xis_trap[5])]) # auto_integrated
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[6], calc_xis_trap[6])]) # newton_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[7], calc_xis_trap[7])]) # doppler_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[8], calc_xis_trap[8])]) # newton_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[9], calc_xis_trap[9])]) # lensing_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[10], calc_xis_trap[10])]) # newton_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[11], calc_xis_trap[11])]) # localgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[12], calc_xis_trap[12])]) # newton_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[13], calc_xis_trap[13])]) # integratedgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[14], calc_xis_trap[14])]) # lensing_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[15], calc_xis_trap[15])]) # doppler_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[16], calc_xis_trap[16])]) # doppler_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[17], calc_xis_trap[17])]) # localgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[18], calc_xis_trap[18])]) # doppler_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[19], calc_xis_trap[19])]) # integratedgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[20], calc_xis_trap[20])]) # lensing_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[21], calc_xis_trap[21])]) # localgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[22], calc_xis_trap[22])]) # lensing_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[23], calc_xis_trap[23])]) # integratedgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[24], calc_xis_trap[24])]) # localgp_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[25], calc_xis_trap[25])]) # integratedgp_localgp

        rm(name)
        for effect in GaPSE.GR_EFFECTS_GNC
            rm("all_GNC_standalones_CFs/xi_GNC_" * effect * "_L$L" * ".txt")
        end
        rm("all_GNC_standalones_CFs")

        #### lob ####
        ss_lob, res_sums_lob, res_xis_lob = GaPSE.readxyall(
            "datatest/GNC_SumXiMultipoles/xis_GNC_L$L" * "_lob_noF_noobs_specific_ss.txt", comments=true)

        name = "calc_xis_GNC_L$L" * "_lob_noF_noobs_specific_ss.txt"
        isfile(name) && rm(name)

        GaPSE.print_map_sum_ξ_GNC_multipole(COSMO, name, ss_lob;
            single=false, L=L, alg=:lobatto, kwargs...)
        calc_ss_lob, calc_sums_lob, calc_xis_lob = GaPSE.readxyall(name, comments=true)
        calc_ss_lob, calc_sums_lob = GaPSE.readxy(name, comments=true)
        calc_xis_lob = [
            begin
                a_name = "all_GNC_standalones_CFs/xi_GNC_" * effect * "_L$L" * ".txt"
                _, xis = GaPSE.readxy(a_name, comments=true)
                xis
            end for effect in GaPSE.GR_EFFECTS_GNC
        ]

        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(ss_lob, calc_ss_lob)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_sums_lob, calc_sums_lob)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[1], calc_xis_lob[1])]) # auto_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[2], calc_xis_lob[2])]) # auto_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[3], calc_xis_lob[3])]) # auto_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[4], calc_xis_lob[4])]) # auto_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[5], calc_xis_lob[5])]) # auto_integrated
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[6], calc_xis_lob[6])]) # newton_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[7], calc_xis_lob[7])]) # doppler_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[8], calc_xis_lob[8])]) # newton_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[9], calc_xis_lob[9])]) # lensing_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[10], calc_xis_lob[10])]) # newton_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[11], calc_xis_lob[11])]) # localgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[12], calc_xis_lob[12])]) # newton_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[13], calc_xis_lob[13])]) # integratedgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[14], calc_xis_lob[14])]) # lensing_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[15], calc_xis_lob[15])]) # doppler_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[16], calc_xis_lob[16])]) # doppler_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[17], calc_xis_lob[17])]) # localgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[18], calc_xis_lob[18])]) # doppler_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[19], calc_xis_lob[19])]) # integratedgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[20], calc_xis_lob[20])]) # lensing_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[21], calc_xis_lob[21])]) # localgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[22], calc_xis_lob[22])]) # lensing_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[23], calc_xis_lob[23])]) # integratedgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[24], calc_xis_lob[24])]) # localgp_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[25], calc_xis_lob[25])]) # integratedgp_localgp

        rm(name)
        for effect in GaPSE.GR_EFFECTS_GNC
            rm("all_GNC_standalones_CFs/xi_GNC_" * effect * "_L$L" * ".txt")
        end
        rm("all_GNC_standalones_CFs")

    end

    @testset "L = 1 no_window no obs" begin
        L = 1

        #### trap ####
        ss_trap, res_sums_trap, res_xis_trap = GaPSE.readxyall(
            "datatest/GNC_SumXiMultipoles/xis_GNC_L$L" * "_trap_noF_noobs_specific_ss.txt", comments=true)

        name = "calc_xis_GNC_L$L" * "_trap_noF_noobs_specific_ss.txt"
        isfile(name) && rm(name)

        GaPSE.print_map_sum_ξ_GNC_multipole(COSMO, name, ss_trap;
            single=false, L=L, alg=:trap, kwargs...)
        calc_ss_trap, calc_sums_trap = GaPSE.readxy(name, comments=true)
        calc_xis_trap = [
            begin
                a_name = "all_GNC_standalones_CFs/xi_GNC_" * effect * "_L$L" * ".txt"
                _, xis = GaPSE.readxy(a_name, comments=true)
                xis
            end for effect in GaPSE.GR_EFFECTS_GNC
        ]

        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(ss_trap, calc_ss_trap)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_sums_trap, calc_sums_trap)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[1], calc_xis_trap[1])]) # auto_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[2], calc_xis_trap[2])]) # auto_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[3], calc_xis_trap[3])]) # auto_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[4], calc_xis_trap[4])]) # auto_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[5], calc_xis_trap[5])]) # auto_integrated
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[6], calc_xis_trap[6])]) # newton_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[7], calc_xis_trap[7])]) # doppler_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[8], calc_xis_trap[8])]) # newton_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[9], calc_xis_trap[9])]) # lensing_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[10], calc_xis_trap[10])]) # newton_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[11], calc_xis_trap[11])]) # localgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[12], calc_xis_trap[12])]) # newton_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[13], calc_xis_trap[13])]) # integratedgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[14], calc_xis_trap[14])]) # lensing_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[15], calc_xis_trap[15])]) # doppler_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[16], calc_xis_trap[16])]) # doppler_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[17], calc_xis_trap[17])]) # localgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[18], calc_xis_trap[18])]) # doppler_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[19], calc_xis_trap[19])]) # integratedgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[20], calc_xis_trap[20])]) # lensing_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[21], calc_xis_trap[21])]) # localgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[22], calc_xis_trap[22])]) # lensing_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[23], calc_xis_trap[23])]) # integratedgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[24], calc_xis_trap[24])]) # localgp_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[25], calc_xis_trap[25])]) # integratedgp_localgp

        rm(name)
        for effect in GaPSE.GR_EFFECTS_GNC
            rm("all_GNC_standalones_CFs/xi_GNC_" * effect * "_L$L" * ".txt")
        end
        rm("all_GNC_standalones_CFs")

        #### lob ####
        ss_lob, res_sums_lob, res_xis_lob = GaPSE.readxyall(
            "datatest/GNC_SumXiMultipoles/xis_GNC_L$L" * "_lob_noF_noobs_specific_ss.txt", comments=true)

        name = "calc_xis_GNC_L$L" * "_lob_noF_noobs_specific_ss.txt"
        isfile(name) && rm(name)

        GaPSE.print_map_sum_ξ_GNC_multipole(COSMO, name, ss_lob;
            single=false, L=L, alg=:lobatto, kwargs...)
        calc_ss_lob, calc_sums_lob = GaPSE.readxy(name, comments=true)
        calc_xis_lob = [
            begin
                a_name = "all_GNC_standalones_CFs/xi_GNC_" * effect * "_L$L" * ".txt"
                _, xis = GaPSE.readxy(a_name, comments=true)
                xis
            end for effect in GaPSE.GR_EFFECTS_GNC
        ]

        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(ss_lob, calc_ss_lob)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_sums_lob, calc_sums_lob)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[1], calc_xis_lob[1])]) # auto_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[2], calc_xis_lob[2])]) # auto_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[3], calc_xis_lob[3])]) # auto_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[4], calc_xis_lob[4])]) # auto_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[5], calc_xis_lob[5])]) # auto_integrated
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[6], calc_xis_lob[6])]) # newton_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[7], calc_xis_lob[7])]) # doppler_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[8], calc_xis_lob[8])]) # newton_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[9], calc_xis_lob[9])]) # lensing_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[10], calc_xis_lob[10])]) # newton_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[11], calc_xis_lob[11])]) # localgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[12], calc_xis_lob[12])]) # newton_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[13], calc_xis_lob[13])]) # integratedgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[14], calc_xis_lob[14])]) # lensing_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[15], calc_xis_lob[15])]) # doppler_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[16], calc_xis_lob[16])]) # doppler_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[17], calc_xis_lob[17])]) # localgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[18], calc_xis_lob[18])]) # doppler_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[19], calc_xis_lob[19])]) # integratedgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[20], calc_xis_lob[20])]) # lensing_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[21], calc_xis_lob[21])]) # localgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[22], calc_xis_lob[22])]) # lensing_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[23], calc_xis_lob[23])]) # integratedgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[24], calc_xis_lob[24])]) # localgp_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[25], calc_xis_lob[25])]) # integratedgp_localgp

        rm(name)
        for effect in GaPSE.GR_EFFECTS_GNC
            rm("all_GNC_standalones_CFs/xi_GNC_" * effect * "_L$L" * ".txt")
        end
        rm("all_GNC_standalones_CFs")

    end

    @testset "L = 2 no_window no obs" begin
        L = 2

        ##### quad ####
        ss_quad, res_sums_quad, res_xis_quad = GaPSE.readxyall(
            "datatest/GNC_SumXiMultipoles/xis_GNC_L$L" * "_quad_noF_noobs_specific_ss.txt", comments=true)

        name = "calc_xis_GNC_L$L" * "_quad_noF_noobs_specific_ss.txt"
        isfile(name) && rm(name)

        GaPSE.print_map_sum_ξ_GNC_multipole(COSMO, name, ss_quad;
            single=false, L=L, alg=:quad, kwargs...)
        calc_ss_quad, calc_sums_quad = GaPSE.readxy(name, comments=true)
        calc_xis_quad = [
            begin
                a_name = "all_GNC_standalones_CFs/xi_GNC_" * effect * "_L$L" * ".txt"
                _, xis = GaPSE.readxy(a_name, comments=true)
                xis
            end for effect in GaPSE.GR_EFFECTS_GNC
        ]

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
        for effect in GaPSE.GR_EFFECTS_GNC
            rm("all_GNC_standalones_CFs/xi_GNC_" * effect * "_L$L" * ".txt")
        end
        rm("all_GNC_standalones_CFs")

        #### trap ####
        ss_trap, res_sums_trap, res_xis_trap = GaPSE.readxyall(
            "datatest/GNC_SumXiMultipoles/xis_GNC_L$L" * "_trap_noF_noobs_specific_ss.txt", comments=true)

        name = "calc_xis_GNC_L$L" * "_trap_noF_noobs_specific_ss.txt"
        isfile(name) && rm(name)

        GaPSE.print_map_sum_ξ_GNC_multipole(COSMO, name, ss_trap;
            single=false, L=L, alg=:trap, kwargs...)
        calc_ss_trap, calc_sums_trap = GaPSE.readxy(name, comments=true)
        calc_xis_trap = [
            begin
                a_name = "all_GNC_standalones_CFs/xi_GNC_" * effect * "_L$L" * ".txt"
                _, xis = GaPSE.readxy(a_name, comments=true)
                xis
            end for effect in GaPSE.GR_EFFECTS_GNC
        ]

        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(ss_trap, calc_ss_trap)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_sums_trap, calc_sums_trap)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[1], calc_xis_trap[1])]) # auto_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[2], calc_xis_trap[2])]) # auto_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[3], calc_xis_trap[3])]) # auto_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[4], calc_xis_trap[4])]) # auto_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[5], calc_xis_trap[5])]) # auto_integrated
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[6], calc_xis_trap[6])]) # newton_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[7], calc_xis_trap[7])]) # doppler_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[8], calc_xis_trap[8])]) # newton_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[9], calc_xis_trap[9])]) # lensing_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[10], calc_xis_trap[10])]) # newton_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[11], calc_xis_trap[11])]) # localgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[12], calc_xis_trap[12])]) # newton_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[13], calc_xis_trap[13])]) # integratedgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[14], calc_xis_trap[14])]) # lensing_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[15], calc_xis_trap[15])]) # doppler_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[16], calc_xis_trap[16])]) # doppler_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[17], calc_xis_trap[17])]) # localgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[18], calc_xis_trap[18])]) # doppler_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[19], calc_xis_trap[19])]) # integratedgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[20], calc_xis_trap[20])]) # lensing_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[21], calc_xis_trap[21])]) # localgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[22], calc_xis_trap[22])]) # lensing_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[23], calc_xis_trap[23])]) # integratedgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[24], calc_xis_trap[24])]) # localgp_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[25], calc_xis_trap[25])]) # integratedgp_localgp

        rm(name)
        for effect in GaPSE.GR_EFFECTS_GNC
            rm("all_GNC_standalones_CFs/xi_GNC_" * effect * "_L$L" * ".txt")
        end
        rm("all_GNC_standalones_CFs")

        #### lob ####
        ss_lob, res_sums_lob, res_xis_lob = GaPSE.readxyall(
            "datatest/GNC_SumXiMultipoles/xis_GNC_L$L" * "_lob_noF_noobs_specific_ss.txt", comments=true)

        name = "calc_xis_GNC_L$L" * "_lob_noF_noobs_specific_ss.txt"
        isfile(name) && rm(name)

        GaPSE.print_map_sum_ξ_GNC_multipole(COSMO, name, ss_lob;
            single=false, L=L, alg=:lobatto, kwargs...)
        calc_ss_lob, calc_sums_lob = GaPSE.readxy(name, comments=true)
        calc_xis_lob = [
            begin
                a_name = "all_GNC_standalones_CFs/xi_GNC_" * effect * "_L$L" * ".txt"
                _, xis = GaPSE.readxy(a_name, comments=true)
                xis
            end for effect in GaPSE.GR_EFFECTS_GNC
        ]

        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(ss_lob, calc_ss_lob)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_sums_lob, calc_sums_lob)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[1], calc_xis_lob[1])]) # auto_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[2], calc_xis_lob[2])]) # auto_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[3], calc_xis_lob[3])]) # auto_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[4], calc_xis_lob[4])]) # auto_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[5], calc_xis_lob[5])]) # auto_integrated
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[6], calc_xis_lob[6])]) # newton_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[7], calc_xis_lob[7])]) # doppler_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[8], calc_xis_lob[8])]) # newton_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[9], calc_xis_lob[9])]) # lensing_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[10], calc_xis_lob[10])]) # newton_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[11], calc_xis_lob[11])]) # localgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[12], calc_xis_lob[12])]) # newton_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[13], calc_xis_lob[13])]) # integratedgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[14], calc_xis_lob[14])]) # lensing_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[15], calc_xis_lob[15])]) # doppler_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[16], calc_xis_lob[16])]) # doppler_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[17], calc_xis_lob[17])]) # localgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[18], calc_xis_lob[18])]) # doppler_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[19], calc_xis_lob[19])]) # integratedgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[20], calc_xis_lob[20])]) # lensing_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[21], calc_xis_lob[21])]) # localgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[22], calc_xis_lob[22])]) # lensing_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[23], calc_xis_lob[23])]) # integratedgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[24], calc_xis_lob[24])]) # localgp_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[25], calc_xis_lob[25])]) # integratedgp_localgp

        rm(name)
        for effect in GaPSE.GR_EFFECTS_GNC
            rm("all_GNC_standalones_CFs/xi_GNC_" * effect * "_L$L" * ".txt")
        end
        rm("all_GNC_standalones_CFs")

    end

    @testset "L = 3 no_window no obs" begin
        L = 3


        #### trap ####
        ss_trap, res_sums_trap, res_xis_trap = GaPSE.readxyall(
            "datatest/GNC_SumXiMultipoles/xis_GNC_L$L" * "_trap_noF_noobs_specific_ss.txt", comments=true)

        name = "calc_xis_GNC_L$L" * "_trap_noF_noobs_specific_ss.txt"
        isfile(name) && rm(name)

        GaPSE.print_map_sum_ξ_GNC_multipole(COSMO, name, ss_trap;
            single=false, L=L, alg=:trap, kwargs...)
        calc_ss_trap, calc_sums_trap = GaPSE.readxy(name, comments=true)
        calc_xis_trap = [
            begin
                a_name = "all_GNC_standalones_CFs/xi_GNC_" * effect * "_L$L" * ".txt"
                _, xis = GaPSE.readxy(a_name, comments=true)
                xis
            end for effect in GaPSE.GR_EFFECTS_GNC
        ]

        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(ss_trap, calc_ss_trap)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_sums_trap, calc_sums_trap)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[1], calc_xis_trap[1])]) # auto_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[2], calc_xis_trap[2])]) # auto_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[3], calc_xis_trap[3])]) # auto_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[4], calc_xis_trap[4])]) # auto_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[5], calc_xis_trap[5])]) # auto_integrated
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[6], calc_xis_trap[6])]) # newton_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[7], calc_xis_trap[7])]) # doppler_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[8], calc_xis_trap[8])]) # newton_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[9], calc_xis_trap[9])]) # lensing_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[10], calc_xis_trap[10])]) # newton_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[11], calc_xis_trap[11])]) # localgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[12], calc_xis_trap[12])]) # newton_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[13], calc_xis_trap[13])]) # integratedgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[14], calc_xis_trap[14])]) # lensing_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[15], calc_xis_trap[15])]) # doppler_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[16], calc_xis_trap[16])]) # doppler_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[17], calc_xis_trap[17])]) # localgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[18], calc_xis_trap[18])]) # doppler_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[19], calc_xis_trap[19])]) # integratedgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[20], calc_xis_trap[20])]) # lensing_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[21], calc_xis_trap[21])]) # localgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[22], calc_xis_trap[22])]) # lensing_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[23], calc_xis_trap[23])]) # integratedgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[24], calc_xis_trap[24])]) # localgp_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[25], calc_xis_trap[25])]) # integratedgp_localgp

        rm(name)
        for effect in GaPSE.GR_EFFECTS_GNC
            rm("all_GNC_standalones_CFs/xi_GNC_" * effect * "_L$L" * ".txt")
        end
        rm("all_GNC_standalones_CFs")

        #### lob ####
        ss_lob, res_sums_lob, res_xis_lob = GaPSE.readxyall(
            "datatest/GNC_SumXiMultipoles/xis_GNC_L$L" * "_lob_noF_noobs_specific_ss.txt", comments=true)

        name = "calc_xis_GNC_L$L" * "_lob_noF_noobs_specific_ss.txt"
        isfile(name) && rm(name)

        GaPSE.print_map_sum_ξ_GNC_multipole(COSMO, name, ss_lob;
            single=false, L=L, alg=:lobatto, kwargs...)
        calc_ss_lob, calc_sums_lob = GaPSE.readxy(name, comments=true)
        calc_xis_lob = [
            begin
                a_name = "all_GNC_standalones_CFs/xi_GNC_" * effect * "_L$L" * ".txt"
                _, xis = GaPSE.readxy(a_name, comments=true)
                xis
            end for effect in GaPSE.GR_EFFECTS_GNC
        ]

        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(ss_lob, calc_ss_lob)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_sums_lob, calc_sums_lob)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[1], calc_xis_lob[1])]) # auto_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[2], calc_xis_lob[2])]) # auto_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[3], calc_xis_lob[3])]) # auto_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[4], calc_xis_lob[4])]) # auto_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[5], calc_xis_lob[5])]) # auto_integrated
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[6], calc_xis_lob[6])]) # newton_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[7], calc_xis_lob[7])]) # doppler_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[8], calc_xis_lob[8])]) # newton_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[9], calc_xis_lob[9])]) # lensing_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[10], calc_xis_lob[10])]) # newton_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[11], calc_xis_lob[11])]) # localgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[12], calc_xis_lob[12])]) # newton_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[13], calc_xis_lob[13])]) # integratedgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[14], calc_xis_lob[14])]) # lensing_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[15], calc_xis_lob[15])]) # doppler_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[16], calc_xis_lob[16])]) # doppler_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[17], calc_xis_lob[17])]) # localgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[18], calc_xis_lob[18])]) # doppler_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[19], calc_xis_lob[19])]) # integratedgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[20], calc_xis_lob[20])]) # lensing_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[21], calc_xis_lob[21])]) # localgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[22], calc_xis_lob[22])]) # lensing_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[23], calc_xis_lob[23])]) # integratedgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[24], calc_xis_lob[24])]) # localgp_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[25], calc_xis_lob[25])]) # integratedgp_localgp

        rm(name)
        for effect in GaPSE.GR_EFFECTS_GNC
            rm("all_GNC_standalones_CFs/xi_GNC_" * effect * "_L$L" * ".txt")
        end
        rm("all_GNC_standalones_CFs")

    end

    @testset "L = 4 no_window no obs" begin
        L = 4

        #### trap ####
        ss_trap, res_sums_trap, res_xis_trap = GaPSE.readxyall(
            "datatest/GNC_SumXiMultipoles/xis_GNC_L$L" * "_trap_noF_noobs_specific_ss.txt", comments=true)

        name = "calc_xis_GNC_L$L" * "_trap_noF_noobs_specific_ss.txt"
        isfile(name) && rm(name)

        GaPSE.print_map_sum_ξ_GNC_multipole(COSMO, name, ss_trap;
            single=false, L=L, alg=:trap, kwargs...)
        calc_ss_trap, calc_sums_trap = GaPSE.readxy(name, comments=true)
        calc_xis_trap = [
            begin
                a_name = "all_GNC_standalones_CFs/xi_GNC_" * effect * "_L$L" * ".txt"
                _, xis = GaPSE.readxy(a_name, comments=true)
                xis
            end for effect in GaPSE.GR_EFFECTS_GNC
        ]

        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(ss_trap, calc_ss_trap)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_sums_trap, calc_sums_trap)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[1], calc_xis_trap[1])]) # auto_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[2], calc_xis_trap[2])]) # auto_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[3], calc_xis_trap[3])]) # auto_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[4], calc_xis_trap[4])]) # auto_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[5], calc_xis_trap[5])]) # auto_integrated
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[6], calc_xis_trap[6])]) # newton_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[7], calc_xis_trap[7])]) # doppler_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[8], calc_xis_trap[8])]) # newton_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[9], calc_xis_trap[9])]) # lensing_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[10], calc_xis_trap[10])]) # newton_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[11], calc_xis_trap[11])]) # localgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[12], calc_xis_trap[12])]) # newton_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[13], calc_xis_trap[13])]) # integratedgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[14], calc_xis_trap[14])]) # lensing_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[15], calc_xis_trap[15])]) # doppler_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[16], calc_xis_trap[16])]) # doppler_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[17], calc_xis_trap[17])]) # localgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[18], calc_xis_trap[18])]) # doppler_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[19], calc_xis_trap[19])]) # integratedgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[20], calc_xis_trap[20])]) # lensing_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[21], calc_xis_trap[21])]) # localgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[22], calc_xis_trap[22])]) # lensing_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[23], calc_xis_trap[23])]) # integratedgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[24], calc_xis_trap[24])]) # localgp_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[25], calc_xis_trap[25])]) # integratedgp_localgp

        rm(name)
        for effect in GaPSE.GR_EFFECTS_GNC
            rm("all_GNC_standalones_CFs/xi_GNC_" * effect * "_L$L" * ".txt")
        end
        rm("all_GNC_standalones_CFs")

        #### lob ####
        ss_lob, res_sums_lob, res_xis_lob = GaPSE.readxyall(
            "datatest/GNC_SumXiMultipoles/xis_GNC_L$L" * "_lob_noF_noobs_specific_ss.txt", comments=true)

        name = "calc_xis_GNC_L$L" * "_lob_noF_noobs_specific_ss.txt"
        isfile(name) && rm(name)

        GaPSE.print_map_sum_ξ_GNC_multipole(COSMO, name, ss_lob;
            single=false, L=L, alg=:lobatto, kwargs...)
        calc_ss_lob, calc_sums_lob = GaPSE.readxy(name, comments=true)
        calc_xis_lob = [
            begin
                a_name = "all_GNC_standalones_CFs/xi_GNC_" * effect * "_L$L" * ".txt"
                _, xis = GaPSE.readxy(a_name, comments=true)
                xis
            end for effect in GaPSE.GR_EFFECTS_GNC
        ]

        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(ss_lob, calc_ss_lob)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_sums_lob, calc_sums_lob)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[1], calc_xis_lob[1])]) # auto_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[2], calc_xis_lob[2])]) # auto_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[3], calc_xis_lob[3])]) # auto_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[4], calc_xis_lob[4])]) # auto_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[5], calc_xis_lob[5])]) # auto_integrated
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[6], calc_xis_lob[6])]) # newton_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[7], calc_xis_lob[7])]) # doppler_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[8], calc_xis_lob[8])]) # newton_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[9], calc_xis_lob[9])]) # lensing_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[10], calc_xis_lob[10])]) # newton_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[11], calc_xis_lob[11])]) # localgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[12], calc_xis_lob[12])]) # newton_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[13], calc_xis_lob[13])]) # integratedgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[14], calc_xis_lob[14])]) # lensing_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[15], calc_xis_lob[15])]) # doppler_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[16], calc_xis_lob[16])]) # doppler_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[17], calc_xis_lob[17])]) # localgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[18], calc_xis_lob[18])]) # doppler_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[19], calc_xis_lob[19])]) # integratedgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[20], calc_xis_lob[20])]) # lensing_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[21], calc_xis_lob[21])]) # localgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[22], calc_xis_lob[22])]) # lensing_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[23], calc_xis_lob[23])]) # integratedgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[24], calc_xis_lob[24])]) # localgp_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[25], calc_xis_lob[25])]) # integratedgp_localgp

        rm(name)
        for effect in GaPSE.GR_EFFECTS_GNC
            rm("all_GNC_standalones_CFs/xi_GNC_" * effect * "_L$L" * ".txt")
        end
        rm("all_GNC_standalones_CFs")

    end
end

println("\n  5 / 8 ...")

@testset "test print_map_sum_ξ_GNC_multipole no_window no obs vel, second method" begin
    RTOL = 1.2e-2

    kwargs = Dict(
        :use_windows => false, :obs => :noobsvel,
        :enhancer => 1e8, :s1 => nothing,
        :N_trap => 30, :N_lob => 30,
        :atol_quad => 0.0, :rtol_quad => 1e-2,
        :N_χs => 40, :N_χs_2 => 20,
        :pr => false,
    )

    @testset "test_zeros" begin
        @test_throws AssertionError GaPSE.GaPSE.print_map_sum_ξ_GNC_multipole(COSMO, "NONEXISTINGDIR/file.txt", [1.0, 2.0, 3.0])
    end

    @testset "L = 0 no_window no obs vel" begin
        L = 0

        ##### quad ####
        ss_quad, res_sums_quad, res_xis_quad = GaPSE.readxyall(
            "datatest/GNC_SumXiMultipoles/xis_GNC_L$L" * "_quad_noF_noobsvel_specific_ss.txt", comments=true)

        name = "calc_xis_GNC_L$L" * "_quad_noF_noobsvel_specific_ss.txt"
        isfile(name) && rm(name)

        GaPSE.print_map_sum_ξ_GNC_multipole(COSMO, name, ss_quad;
            single=false, L=L, alg=:quad, kwargs...)

        calc_ss_quad, calc_sums_quad = GaPSE.readxy(name, comments=true)
        calc_xis_quad = [
            begin
                a_name = "all_GNC_standalones_CFs/xi_GNC_" * effect * "_L$L" * ".txt"
                _, xis = GaPSE.readxy(a_name, comments=true)
                xis
            end for effect in GaPSE.GR_EFFECTS_GNC
        ]

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
        for effect in GaPSE.GR_EFFECTS_GNC
            rm("all_GNC_standalones_CFs/xi_GNC_" * effect * "_L$L" * ".txt")
        end
        rm("all_GNC_standalones_CFs")

        #### trap ####
        ss_trap, res_sums_trap, res_xis_trap = GaPSE.readxyall(
            "datatest/GNC_SumXiMultipoles/xis_GNC_L$L" * "_trap_noF_noobsvel_specific_ss.txt", comments=true)

        name = "calc_xis_GNC_L$L" * "_trap_noF_noobsvel_specific_ss.txt"
        isfile(name) && rm(name)

        GaPSE.print_map_sum_ξ_GNC_multipole(COSMO, name, ss_trap;
            single=false, L=L, alg=:trap, kwargs...)
        calc_ss_trap, calc_sums_trap = GaPSE.readxy(name, comments=true)
        calc_xis_trap = [
            begin
                a_name = "all_GNC_standalones_CFs/xi_GNC_" * effect * "_L$L" * ".txt"
                _, xis = GaPSE.readxy(a_name, comments=true)
                xis
            end for effect in GaPSE.GR_EFFECTS_GNC
        ]

        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(ss_trap, calc_ss_trap)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_sums_trap, calc_sums_trap)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[1], calc_xis_trap[1])]) # auto_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[2], calc_xis_trap[2])]) # auto_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[3], calc_xis_trap[3])]) # auto_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[4], calc_xis_trap[4])]) # auto_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[5], calc_xis_trap[5])]) # auto_integrated
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[6], calc_xis_trap[6])]) # newton_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[7], calc_xis_trap[7])]) # doppler_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[8], calc_xis_trap[8])]) # newton_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[9], calc_xis_trap[9])]) # lensing_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[10], calc_xis_trap[10])]) # newton_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[11], calc_xis_trap[11])]) # localgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[12], calc_xis_trap[12])]) # newton_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[13], calc_xis_trap[13])]) # integratedgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[14], calc_xis_trap[14])]) # lensing_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[15], calc_xis_trap[15])]) # doppler_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[16], calc_xis_trap[16])]) # doppler_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[17], calc_xis_trap[17])]) # localgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[18], calc_xis_trap[18])]) # doppler_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[19], calc_xis_trap[19])]) # integratedgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[20], calc_xis_trap[20])]) # lensing_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[21], calc_xis_trap[21])]) # localgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[22], calc_xis_trap[22])]) # lensing_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[23], calc_xis_trap[23])]) # integratedgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[24], calc_xis_trap[24])]) # localgp_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[25], calc_xis_trap[25])]) # integratedgp_localgp

        rm(name)
        for effect in GaPSE.GR_EFFECTS_GNC
            rm("all_GNC_standalones_CFs/xi_GNC_" * effect * "_L$L" * ".txt")
        end
        rm("all_GNC_standalones_CFs")

        #### lob ####
        ss_lob, res_sums_lob, res_xis_lob = GaPSE.readxyall(
            "datatest/GNC_SumXiMultipoles/xis_GNC_L$L" * "_lob_noF_noobsvel_specific_ss.txt", comments=true)

        name = "calc_xis_GNC_L$L" * "_lob_noF_noobsvel_specific_ss.txt"
        isfile(name) && rm(name)

        GaPSE.print_map_sum_ξ_GNC_multipole(COSMO, name, ss_lob;
            single=false, L=L, alg=:lobatto, kwargs...)
        calc_ss_lob, calc_sums_lob, calc_xis_lob = GaPSE.readxyall(name, comments=true)
        calc_ss_lob, calc_sums_lob = GaPSE.readxy(name, comments=true)
        calc_xis_lob = [
            begin
                a_name = "all_GNC_standalones_CFs/xi_GNC_" * effect * "_L$L" * ".txt"
                _, xis = GaPSE.readxy(a_name, comments=true)
                xis
            end for effect in GaPSE.GR_EFFECTS_GNC
        ]

        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(ss_lob, calc_ss_lob)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_sums_lob, calc_sums_lob)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[1], calc_xis_lob[1])]) # auto_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[2], calc_xis_lob[2])]) # auto_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[3], calc_xis_lob[3])]) # auto_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[4], calc_xis_lob[4])]) # auto_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[5], calc_xis_lob[5])]) # auto_integrated
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[6], calc_xis_lob[6])]) # newton_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[7], calc_xis_lob[7])]) # doppler_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[8], calc_xis_lob[8])]) # newton_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[9], calc_xis_lob[9])]) # lensing_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[10], calc_xis_lob[10])]) # newton_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[11], calc_xis_lob[11])]) # localgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[12], calc_xis_lob[12])]) # newton_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[13], calc_xis_lob[13])]) # integratedgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[14], calc_xis_lob[14])]) # lensing_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[15], calc_xis_lob[15])]) # doppler_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[16], calc_xis_lob[16])]) # doppler_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[17], calc_xis_lob[17])]) # localgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[18], calc_xis_lob[18])]) # doppler_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[19], calc_xis_lob[19])]) # integratedgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[20], calc_xis_lob[20])]) # lensing_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[21], calc_xis_lob[21])]) # localgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[22], calc_xis_lob[22])]) # lensing_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[23], calc_xis_lob[23])]) # integratedgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[24], calc_xis_lob[24])]) # localgp_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[25], calc_xis_lob[25])]) # integratedgp_localgp

        rm(name)
        for effect in GaPSE.GR_EFFECTS_GNC
            rm("all_GNC_standalones_CFs/xi_GNC_" * effect * "_L$L" * ".txt")
        end
        rm("all_GNC_standalones_CFs")

    end

    @testset "L = 1 no_window no obs vel" begin
        L = 1

        #### trap ####
        ss_trap, res_sums_trap, res_xis_trap = GaPSE.readxyall(
            "datatest/GNC_SumXiMultipoles/xis_GNC_L$L" * "_trap_noF_noobsvel_specific_ss.txt", comments=true)

        name = "calc_xis_GNC_L$L" * "_trap_noF_noobsvel_specific_ss.txt"
        isfile(name) && rm(name)

        GaPSE.print_map_sum_ξ_GNC_multipole(COSMO, name, ss_trap;
            single=false, L=L, alg=:trap, kwargs...)
        calc_ss_trap, calc_sums_trap = GaPSE.readxy(name, comments=true)
        calc_xis_trap = [
            begin
                a_name = "all_GNC_standalones_CFs/xi_GNC_" * effect * "_L$L" * ".txt"
                _, xis = GaPSE.readxy(a_name, comments=true)
                xis
            end for effect in GaPSE.GR_EFFECTS_GNC
        ]

        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(ss_trap, calc_ss_trap)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_sums_trap, calc_sums_trap)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[1], calc_xis_trap[1])]) # auto_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[2], calc_xis_trap[2])]) # auto_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[3], calc_xis_trap[3])]) # auto_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[4], calc_xis_trap[4])]) # auto_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[5], calc_xis_trap[5])]) # auto_integrated
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[6], calc_xis_trap[6])]) # newton_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[7], calc_xis_trap[7])]) # doppler_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[8], calc_xis_trap[8])]) # newton_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[9], calc_xis_trap[9])]) # lensing_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[10], calc_xis_trap[10])]) # newton_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[11], calc_xis_trap[11])]) # localgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[12], calc_xis_trap[12])]) # newton_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[13], calc_xis_trap[13])]) # integratedgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[14], calc_xis_trap[14])]) # lensing_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[15], calc_xis_trap[15])]) # doppler_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[16], calc_xis_trap[16])]) # doppler_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[17], calc_xis_trap[17])]) # localgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[18], calc_xis_trap[18])]) # doppler_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[19], calc_xis_trap[19])]) # integratedgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[20], calc_xis_trap[20])]) # lensing_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[21], calc_xis_trap[21])]) # localgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[22], calc_xis_trap[22])]) # lensing_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[23], calc_xis_trap[23])]) # integratedgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[24], calc_xis_trap[24])]) # localgp_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[25], calc_xis_trap[25])]) # integratedgp_localgp

        rm(name)
        for effect in GaPSE.GR_EFFECTS_GNC
            rm("all_GNC_standalones_CFs/xi_GNC_" * effect * "_L$L" * ".txt")
        end
        rm("all_GNC_standalones_CFs")

        #### lob ####
        ss_lob, res_sums_lob, res_xis_lob = GaPSE.readxyall(
            "datatest/GNC_SumXiMultipoles/xis_GNC_L$L" * "_lob_noF_noobsvel_specific_ss.txt", comments=true)

        name = "calc_xis_GNC_L$L" * "_lob_noF_noobsvel_specific_ss.txt"
        isfile(name) && rm(name)

        GaPSE.print_map_sum_ξ_GNC_multipole(COSMO, name, ss_lob;
            single=false, L=L, alg=:lobatto, kwargs...)
        calc_ss_lob, calc_sums_lob = GaPSE.readxy(name, comments=true)
        calc_xis_lob = [
            begin
                a_name = "all_GNC_standalones_CFs/xi_GNC_" * effect * "_L$L" * ".txt"
                _, xis = GaPSE.readxy(a_name, comments=true)
                xis
            end for effect in GaPSE.GR_EFFECTS_GNC
        ]

        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(ss_lob, calc_ss_lob)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_sums_lob, calc_sums_lob)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[1], calc_xis_lob[1])]) # auto_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[2], calc_xis_lob[2])]) # auto_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[3], calc_xis_lob[3])]) # auto_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[4], calc_xis_lob[4])]) # auto_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[5], calc_xis_lob[5])]) # auto_integrated
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[6], calc_xis_lob[6])]) # newton_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[7], calc_xis_lob[7])]) # doppler_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[8], calc_xis_lob[8])]) # newton_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[9], calc_xis_lob[9])]) # lensing_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[10], calc_xis_lob[10])]) # newton_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[11], calc_xis_lob[11])]) # localgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[12], calc_xis_lob[12])]) # newton_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[13], calc_xis_lob[13])]) # integratedgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[14], calc_xis_lob[14])]) # lensing_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[15], calc_xis_lob[15])]) # doppler_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[16], calc_xis_lob[16])]) # doppler_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[17], calc_xis_lob[17])]) # localgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[18], calc_xis_lob[18])]) # doppler_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[19], calc_xis_lob[19])]) # integratedgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[20], calc_xis_lob[20])]) # lensing_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[21], calc_xis_lob[21])]) # localgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[22], calc_xis_lob[22])]) # lensing_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[23], calc_xis_lob[23])]) # integratedgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[24], calc_xis_lob[24])]) # localgp_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[25], calc_xis_lob[25])]) # integratedgp_localgp

        rm(name)
        for effect in GaPSE.GR_EFFECTS_GNC
            rm("all_GNC_standalones_CFs/xi_GNC_" * effect * "_L$L" * ".txt")
        end
        rm("all_GNC_standalones_CFs")

    end

    @testset "L = 2 no_window no obs vel" begin
        L = 2

        ##### quad ####
        ss_quad, res_sums_quad, res_xis_quad = GaPSE.readxyall(
            "datatest/GNC_SumXiMultipoles/xis_GNC_L$L" * "_quad_noF_noobsvel_specific_ss.txt", comments=true)

        name = "calc_xis_GNC_L$L" * "_quad_noF_noobsvel_specific_ss.txt"
        isfile(name) && rm(name)

        GaPSE.print_map_sum_ξ_GNC_multipole(COSMO, name, ss_quad;
            single=false, L=L, alg=:quad, kwargs...)
        calc_ss_quad, calc_sums_quad = GaPSE.readxy(name, comments=true)
        calc_xis_quad = [
            begin
                a_name = "all_GNC_standalones_CFs/xi_GNC_" * effect * "_L$L" * ".txt"
                _, xis = GaPSE.readxy(a_name, comments=true)
                xis
            end for effect in GaPSE.GR_EFFECTS_GNC
        ]

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
        for effect in GaPSE.GR_EFFECTS_GNC
            rm("all_GNC_standalones_CFs/xi_GNC_" * effect * "_L$L" * ".txt")
        end
        rm("all_GNC_standalones_CFs")

        #### trap ####
        ss_trap, res_sums_trap, res_xis_trap = GaPSE.readxyall(
            "datatest/GNC_SumXiMultipoles/xis_GNC_L$L" * "_trap_noF_noobsvel_specific_ss.txt", comments=true)

        name = "calc_xis_GNC_L$L" * "_trap_noF_noobsvel_specific_ss.txt"
        isfile(name) && rm(name)

        GaPSE.print_map_sum_ξ_GNC_multipole(COSMO, name, ss_trap;
            single=false, L=L, alg=:trap, kwargs...)
        calc_ss_trap, calc_sums_trap = GaPSE.readxy(name, comments=true)
        calc_xis_trap = [
            begin
                a_name = "all_GNC_standalones_CFs/xi_GNC_" * effect * "_L$L" * ".txt"
                _, xis = GaPSE.readxy(a_name, comments=true)
                xis
            end for effect in GaPSE.GR_EFFECTS_GNC
        ]

        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(ss_trap, calc_ss_trap)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_sums_trap, calc_sums_trap)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[1], calc_xis_trap[1])]) # auto_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[2], calc_xis_trap[2])]) # auto_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[3], calc_xis_trap[3])]) # auto_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[4], calc_xis_trap[4])]) # auto_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[5], calc_xis_trap[5])]) # auto_integrated
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[6], calc_xis_trap[6])]) # newton_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[7], calc_xis_trap[7])]) # doppler_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[8], calc_xis_trap[8])]) # newton_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[9], calc_xis_trap[9])]) # lensing_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[10], calc_xis_trap[10])]) # newton_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[11], calc_xis_trap[11])]) # localgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[12], calc_xis_trap[12])]) # newton_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[13], calc_xis_trap[13])]) # integratedgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[14], calc_xis_trap[14])]) # lensing_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[15], calc_xis_trap[15])]) # doppler_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[16], calc_xis_trap[16])]) # doppler_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[17], calc_xis_trap[17])]) # localgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[18], calc_xis_trap[18])]) # doppler_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[19], calc_xis_trap[19])]) # integratedgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[20], calc_xis_trap[20])]) # lensing_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[21], calc_xis_trap[21])]) # localgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[22], calc_xis_trap[22])]) # lensing_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[23], calc_xis_trap[23])]) # integratedgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[24], calc_xis_trap[24])]) # localgp_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[25], calc_xis_trap[25])]) # integratedgp_localgp

        rm(name)
        for effect in GaPSE.GR_EFFECTS_GNC
            rm("all_GNC_standalones_CFs/xi_GNC_" * effect * "_L$L" * ".txt")
        end
        rm("all_GNC_standalones_CFs")

        #### lob ####
        ss_lob, res_sums_lob, res_xis_lob = GaPSE.readxyall(
            "datatest/GNC_SumXiMultipoles/xis_GNC_L$L" * "_lob_noF_noobsvel_specific_ss.txt", comments=true)

        name = "calc_xis_GNC_L$L" * "_lob_noF_noobsvel_specific_ss.txt"
        isfile(name) && rm(name)

        GaPSE.print_map_sum_ξ_GNC_multipole(COSMO, name, ss_lob;
            single=false, L=L, alg=:lobatto, kwargs...)
        calc_ss_lob, calc_sums_lob = GaPSE.readxy(name, comments=true)
        calc_xis_lob = [
            begin
                a_name = "all_GNC_standalones_CFs/xi_GNC_" * effect * "_L$L" * ".txt"
                _, xis = GaPSE.readxy(a_name, comments=true)
                xis
            end for effect in GaPSE.GR_EFFECTS_GNC
        ]

        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(ss_lob, calc_ss_lob)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_sums_lob, calc_sums_lob)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[1], calc_xis_lob[1])]) # auto_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[2], calc_xis_lob[2])]) # auto_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[3], calc_xis_lob[3])]) # auto_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[4], calc_xis_lob[4])]) # auto_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[5], calc_xis_lob[5])]) # auto_integrated
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[6], calc_xis_lob[6])]) # newton_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[7], calc_xis_lob[7])]) # doppler_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[8], calc_xis_lob[8])]) # newton_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[9], calc_xis_lob[9])]) # lensing_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[10], calc_xis_lob[10])]) # newton_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[11], calc_xis_lob[11])]) # localgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[12], calc_xis_lob[12])]) # newton_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[13], calc_xis_lob[13])]) # integratedgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[14], calc_xis_lob[14])]) # lensing_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[15], calc_xis_lob[15])]) # doppler_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[16], calc_xis_lob[16])]) # doppler_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[17], calc_xis_lob[17])]) # localgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[18], calc_xis_lob[18])]) # doppler_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[19], calc_xis_lob[19])]) # integratedgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[20], calc_xis_lob[20])]) # lensing_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[21], calc_xis_lob[21])]) # localgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[22], calc_xis_lob[22])]) # lensing_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[23], calc_xis_lob[23])]) # integratedgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[24], calc_xis_lob[24])]) # localgp_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[25], calc_xis_lob[25])]) # integratedgp_localgp

        rm(name)
        for effect in GaPSE.GR_EFFECTS_GNC
            rm("all_GNC_standalones_CFs/xi_GNC_" * effect * "_L$L" * ".txt")
        end
        rm("all_GNC_standalones_CFs")

    end

    @testset "L = 3 no_window no obs vel" begin
        L = 3

        #### trap ####
        ss_trap, res_sums_trap, res_xis_trap = GaPSE.readxyall(
            "datatest/GNC_SumXiMultipoles/xis_GNC_L$L" * "_trap_noF_noobsvel_specific_ss.txt", comments=true)

        name = "calc_xis_GNC_L$L" * "_trap_noF_noobsvel_specific_ss.txt"
        isfile(name) && rm(name)

        GaPSE.print_map_sum_ξ_GNC_multipole(COSMO, name, ss_trap;
            single=false, L=L, alg=:trap, kwargs...)
        calc_ss_trap, calc_sums_trap = GaPSE.readxy(name, comments=true)
        calc_xis_trap = [
            begin
                a_name = "all_GNC_standalones_CFs/xi_GNC_" * effect * "_L$L" * ".txt"
                _, xis = GaPSE.readxy(a_name, comments=true)
                xis
            end for effect in GaPSE.GR_EFFECTS_GNC
        ]

        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(ss_trap, calc_ss_trap)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_sums_trap, calc_sums_trap)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[1], calc_xis_trap[1])]) # auto_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[2], calc_xis_trap[2])]) # auto_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[3], calc_xis_trap[3])]) # auto_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[4], calc_xis_trap[4])]) # auto_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[5], calc_xis_trap[5])]) # auto_integrated
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[6], calc_xis_trap[6])]) # newton_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[7], calc_xis_trap[7])]) # doppler_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[8], calc_xis_trap[8])]) # newton_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[9], calc_xis_trap[9])]) # lensing_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[10], calc_xis_trap[10])]) # newton_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[11], calc_xis_trap[11])]) # localgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[12], calc_xis_trap[12])]) # newton_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[13], calc_xis_trap[13])]) # integratedgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[14], calc_xis_trap[14])]) # lensing_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[15], calc_xis_trap[15])]) # doppler_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[16], calc_xis_trap[16])]) # doppler_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[17], calc_xis_trap[17])]) # localgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[18], calc_xis_trap[18])]) # doppler_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[19], calc_xis_trap[19])]) # integratedgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[20], calc_xis_trap[20])]) # lensing_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[21], calc_xis_trap[21])]) # localgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[22], calc_xis_trap[22])]) # lensing_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[23], calc_xis_trap[23])]) # integratedgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[24], calc_xis_trap[24])]) # localgp_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[25], calc_xis_trap[25])]) # integratedgp_localgp

        rm(name)
        for effect in GaPSE.GR_EFFECTS_GNC
            rm("all_GNC_standalones_CFs/xi_GNC_" * effect * "_L$L" * ".txt")
        end
        rm("all_GNC_standalones_CFs")

        #### lob ####
        ss_lob, res_sums_lob, res_xis_lob = GaPSE.readxyall(
            "datatest/GNC_SumXiMultipoles/xis_GNC_L$L" * "_lob_noF_noobsvel_specific_ss.txt", comments=true)

        name = "calc_xis_GNC_L$L" * "_lob_noF_noobsvel_specific_ss.txt"
        isfile(name) && rm(name)

        GaPSE.print_map_sum_ξ_GNC_multipole(COSMO, name, ss_lob;
            single=false, L=L, alg=:lobatto, kwargs...)
        calc_ss_lob, calc_sums_lob = GaPSE.readxy(name, comments=true)
        calc_xis_lob = [
            begin
                a_name = "all_GNC_standalones_CFs/xi_GNC_" * effect * "_L$L" * ".txt"
                _, xis = GaPSE.readxy(a_name, comments=true)
                xis
            end for effect in GaPSE.GR_EFFECTS_GNC
        ]

        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(ss_lob, calc_ss_lob)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_sums_lob, calc_sums_lob)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[1], calc_xis_lob[1])]) # auto_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[2], calc_xis_lob[2])]) # auto_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[3], calc_xis_lob[3])]) # auto_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[4], calc_xis_lob[4])]) # auto_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[5], calc_xis_lob[5])]) # auto_integrated
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[6], calc_xis_lob[6])]) # newton_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[7], calc_xis_lob[7])]) # doppler_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[8], calc_xis_lob[8])]) # newton_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[9], calc_xis_lob[9])]) # lensing_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[10], calc_xis_lob[10])]) # newton_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[11], calc_xis_lob[11])]) # localgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[12], calc_xis_lob[12])]) # newton_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[13], calc_xis_lob[13])]) # integratedgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[14], calc_xis_lob[14])]) # lensing_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[15], calc_xis_lob[15])]) # doppler_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[16], calc_xis_lob[16])]) # doppler_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[17], calc_xis_lob[17])]) # localgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[18], calc_xis_lob[18])]) # doppler_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[19], calc_xis_lob[19])]) # integratedgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[20], calc_xis_lob[20])]) # lensing_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[21], calc_xis_lob[21])]) # localgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[22], calc_xis_lob[22])]) # lensing_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[23], calc_xis_lob[23])]) # integratedgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[24], calc_xis_lob[24])]) # localgp_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[25], calc_xis_lob[25])]) # integratedgp_localgp

        rm(name)
        for effect in GaPSE.GR_EFFECTS_GNC
            rm("all_GNC_standalones_CFs/xi_GNC_" * effect * "_L$L" * ".txt")
        end
        rm("all_GNC_standalones_CFs")

    end

    @testset "L = 4 no_window no obs vel" begin
        L = 4

        #### trap ####
        ss_trap, res_sums_trap, res_xis_trap = GaPSE.readxyall(
            "datatest/GNC_SumXiMultipoles/xis_GNC_L$L" * "_trap_noF_noobsvel_specific_ss.txt", comments=true)

        name = "calc_xis_GNC_L$L" * "_trap_noF_noobsvel_specific_ss.txt"
        isfile(name) && rm(name)

        GaPSE.print_map_sum_ξ_GNC_multipole(COSMO, name, ss_trap;
            single=false, L=L, alg=:trap, kwargs...)
        calc_ss_trap, calc_sums_trap = GaPSE.readxy(name, comments=true)
        calc_xis_trap = [
            begin
                a_name = "all_GNC_standalones_CFs/xi_GNC_" * effect * "_L$L" * ".txt"
                _, xis = GaPSE.readxy(a_name, comments=true)
                xis
            end for effect in GaPSE.GR_EFFECTS_GNC
        ]

        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(ss_trap, calc_ss_trap)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_sums_trap, calc_sums_trap)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[1], calc_xis_trap[1])]) # auto_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[2], calc_xis_trap[2])]) # auto_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[3], calc_xis_trap[3])]) # auto_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[4], calc_xis_trap[4])]) # auto_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[5], calc_xis_trap[5])]) # auto_integrated
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[6], calc_xis_trap[6])]) # newton_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[7], calc_xis_trap[7])]) # doppler_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[8], calc_xis_trap[8])]) # newton_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[9], calc_xis_trap[9])]) # lensing_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[10], calc_xis_trap[10])]) # newton_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[11], calc_xis_trap[11])]) # localgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[12], calc_xis_trap[12])]) # newton_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[13], calc_xis_trap[13])]) # integratedgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[14], calc_xis_trap[14])]) # lensing_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[15], calc_xis_trap[15])]) # doppler_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[16], calc_xis_trap[16])]) # doppler_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[17], calc_xis_trap[17])]) # localgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[18], calc_xis_trap[18])]) # doppler_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[19], calc_xis_trap[19])]) # integratedgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[20], calc_xis_trap[20])]) # lensing_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[21], calc_xis_trap[21])]) # localgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[22], calc_xis_trap[22])]) # lensing_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[23], calc_xis_trap[23])]) # integratedgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[24], calc_xis_trap[24])]) # localgp_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[25], calc_xis_trap[25])]) # integratedgp_localgp

        rm(name)
        for effect in GaPSE.GR_EFFECTS_GNC
            rm("all_GNC_standalones_CFs/xi_GNC_" * effect * "_L$L" * ".txt")
        end
        rm("all_GNC_standalones_CFs")

        #### lob ####
        ss_lob, res_sums_lob, res_xis_lob = GaPSE.readxyall(
            "datatest/GNC_SumXiMultipoles/xis_GNC_L$L" * "_lob_noF_noobsvel_specific_ss.txt", comments=true)

        name = "calc_xis_GNC_L$L" * "_lob_noF_noobsvel_specific_ss.txt"
        isfile(name) && rm(name)

        GaPSE.print_map_sum_ξ_GNC_multipole(COSMO, name, ss_lob;
            single=false, L=L, alg=:lobatto, kwargs...)
        calc_ss_lob, calc_sums_lob = GaPSE.readxy(name, comments=true)
        calc_xis_lob = [
            begin
                a_name = "all_GNC_standalones_CFs/xi_GNC_" * effect * "_L$L" * ".txt"
                _, xis = GaPSE.readxy(a_name, comments=true)
                xis
            end for effect in GaPSE.GR_EFFECTS_GNC
        ]

        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(ss_lob, calc_ss_lob)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_sums_lob, calc_sums_lob)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[1], calc_xis_lob[1])]) # auto_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[2], calc_xis_lob[2])]) # auto_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[3], calc_xis_lob[3])]) # auto_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[4], calc_xis_lob[4])]) # auto_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[5], calc_xis_lob[5])]) # auto_integrated
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[6], calc_xis_lob[6])]) # newton_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[7], calc_xis_lob[7])]) # doppler_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[8], calc_xis_lob[8])]) # newton_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[9], calc_xis_lob[9])]) # lensing_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[10], calc_xis_lob[10])]) # newton_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[11], calc_xis_lob[11])]) # localgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[12], calc_xis_lob[12])]) # newton_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[13], calc_xis_lob[13])]) # integratedgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[14], calc_xis_lob[14])]) # lensing_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[15], calc_xis_lob[15])]) # doppler_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[16], calc_xis_lob[16])]) # doppler_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[17], calc_xis_lob[17])]) # localgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[18], calc_xis_lob[18])]) # doppler_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[19], calc_xis_lob[19])]) # integratedgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[20], calc_xis_lob[20])]) # lensing_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[21], calc_xis_lob[21])]) # localgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[22], calc_xis_lob[22])]) # lensing_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[23], calc_xis_lob[23])]) # integratedgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[24], calc_xis_lob[24])]) # localgp_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[25], calc_xis_lob[25])]) # integratedgp_localgp

        rm(name)
        for effect in GaPSE.GR_EFFECTS_GNC
            rm("all_GNC_standalones_CFs/xi_GNC_" * effect * "_L$L" * ".txt")
        end
        rm("all_GNC_standalones_CFs")

    end
end

println("\n  6 / 8 ...")

@testset "test print_map_sum_ξ_GNC_multipole no_window with obs, second method" begin
    RTOL = 1.2e-2

    kwargs = Dict(
        :use_windows => false, :obs => :yes,
        :enhancer => 1e8, :s1 => nothing,
        :N_trap => 30, :N_lob => 30,
        :atol_quad => 0.0, :rtol_quad => 1e-2,
        :N_χs => 40, :N_χs_2 => 20,
        :pr => false,
    )

    @testset "test_zeros" begin
        @test_throws AssertionError GaPSE.GaPSE.print_map_sum_ξ_GNC_multipole(COSMO, "NONEXISTINGDIR/file.txt", [1.0, 2.0, 3.0])
    end

    @testset "L = 0 no_window with obs" begin
        L = 0

        ##### quad ####
        ss_quad, res_sums_quad, res_xis_quad = GaPSE.readxyall(
            "datatest/GNC_SumXiMultipoles/xis_GNC_L$L" * "_quad_noF_withobs_specific_ss.txt", comments=true)

        name = "calc_xis_GNC_L$L" * "_quad_noF_withobs_specific_ss.txt"
        isfile(name) && rm(name)

        GaPSE.print_map_sum_ξ_GNC_multipole(COSMO, name, ss_quad;
            single=false, L=L, alg=:quad, kwargs...)

        calc_ss_quad, calc_sums_quad = GaPSE.readxy(name, comments=true)
        calc_xis_quad = [
            begin
                a_name = "all_GNC_standalones_CFs/xi_GNC_" * effect * "_L$L" * ".txt"
                _, xis = GaPSE.readxy(a_name, comments=true)
                xis
            end for effect in GaPSE.GR_EFFECTS_GNC
        ]

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
        for effect in GaPSE.GR_EFFECTS_GNC
            rm("all_GNC_standalones_CFs/xi_GNC_" * effect * "_L$L" * ".txt")
        end
        rm("all_GNC_standalones_CFs")

        #### trap ####
        ss_trap, res_sums_trap, res_xis_trap = GaPSE.readxyall(
            "datatest/GNC_SumXiMultipoles/xis_GNC_L$L" * "_trap_noF_withobs_specific_ss.txt", comments=true)

        name = "calc_xis_GNC_L$L" * "_trap_noF_withobs_specific_ss.txt"
        isfile(name) && rm(name)

        GaPSE.print_map_sum_ξ_GNC_multipole(COSMO, name, ss_trap;
            single=false, L=L, alg=:trap, kwargs...)
        calc_ss_trap, calc_sums_trap = GaPSE.readxy(name, comments=true)
        calc_xis_trap = [
            begin
                a_name = "all_GNC_standalones_CFs/xi_GNC_" * effect * "_L$L" * ".txt"
                _, xis = GaPSE.readxy(a_name, comments=true)
                xis
            end for effect in GaPSE.GR_EFFECTS_GNC
        ]

        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(ss_trap, calc_ss_trap)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_sums_trap, calc_sums_trap)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[1], calc_xis_trap[1])]) # auto_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[2], calc_xis_trap[2])]) # auto_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[3], calc_xis_trap[3])]) # auto_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[4], calc_xis_trap[4])]) # auto_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[5], calc_xis_trap[5])]) # auto_integrated
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[6], calc_xis_trap[6])]) # newton_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[7], calc_xis_trap[7])]) # doppler_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[8], calc_xis_trap[8])]) # newton_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[9], calc_xis_trap[9])]) # lensing_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[10], calc_xis_trap[10])]) # newton_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[11], calc_xis_trap[11])]) # localgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[12], calc_xis_trap[12])]) # newton_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[13], calc_xis_trap[13])]) # integratedgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[14], calc_xis_trap[14])]) # lensing_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[15], calc_xis_trap[15])]) # doppler_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[16], calc_xis_trap[16])]) # doppler_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[17], calc_xis_trap[17])]) # localgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[18], calc_xis_trap[18])]) # doppler_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[19], calc_xis_trap[19])]) # integratedgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[20], calc_xis_trap[20])]) # lensing_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[21], calc_xis_trap[21])]) # localgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[22], calc_xis_trap[22])]) # lensing_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[23], calc_xis_trap[23])]) # integratedgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[24], calc_xis_trap[24])]) # localgp_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[25], calc_xis_trap[25])]) # integratedgp_localgp

        rm(name)
        for effect in GaPSE.GR_EFFECTS_GNC
            rm("all_GNC_standalones_CFs/xi_GNC_" * effect * "_L$L" * ".txt")
        end
        rm("all_GNC_standalones_CFs")

        #### lob ####
        ss_lob, res_sums_lob, res_xis_lob = GaPSE.readxyall(
            "datatest/GNC_SumXiMultipoles/xis_GNC_L$L" * "_lob_noF_withobs_specific_ss.txt", comments=true)

        name = "calc_xis_GNC_L$L" * "_lob_noF_withobs_specific_ss.txt"
        isfile(name) && rm(name)

        GaPSE.print_map_sum_ξ_GNC_multipole(COSMO, name, ss_lob;
            single=false, L=L, alg=:lobatto, kwargs...)
        calc_ss_lob, calc_sums_lob, calc_xis_lob = GaPSE.readxyall(name, comments=true)
        calc_ss_lob, calc_sums_lob = GaPSE.readxy(name, comments=true)
        calc_xis_lob = [
            begin
                a_name = "all_GNC_standalones_CFs/xi_GNC_" * effect * "_L$L" * ".txt"
                _, xis = GaPSE.readxy(a_name, comments=true)
                xis
            end for effect in GaPSE.GR_EFFECTS_GNC
        ]

        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(ss_lob, calc_ss_lob)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_sums_lob, calc_sums_lob)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[1], calc_xis_lob[1])]) # auto_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[2], calc_xis_lob[2])]) # auto_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[3], calc_xis_lob[3])]) # auto_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[4], calc_xis_lob[4])]) # auto_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[5], calc_xis_lob[5])]) # auto_integrated
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[6], calc_xis_lob[6])]) # newton_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[7], calc_xis_lob[7])]) # doppler_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[8], calc_xis_lob[8])]) # newton_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[9], calc_xis_lob[9])]) # lensing_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[10], calc_xis_lob[10])]) # newton_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[11], calc_xis_lob[11])]) # localgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[12], calc_xis_lob[12])]) # newton_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[13], calc_xis_lob[13])]) # integratedgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[14], calc_xis_lob[14])]) # lensing_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[15], calc_xis_lob[15])]) # doppler_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[16], calc_xis_lob[16])]) # doppler_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[17], calc_xis_lob[17])]) # localgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[18], calc_xis_lob[18])]) # doppler_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[19], calc_xis_lob[19])]) # integratedgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[20], calc_xis_lob[20])]) # lensing_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[21], calc_xis_lob[21])]) # localgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[22], calc_xis_lob[22])]) # lensing_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[23], calc_xis_lob[23])]) # integratedgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[24], calc_xis_lob[24])]) # localgp_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[25], calc_xis_lob[25])]) # integratedgp_localgp

        rm(name)
        for effect in GaPSE.GR_EFFECTS_GNC
            rm("all_GNC_standalones_CFs/xi_GNC_" * effect * "_L$L" * ".txt")
        end
        rm("all_GNC_standalones_CFs")

    end

    @testset "L = 1 no_window with obs" begin
        L = 1

        #### trap ####
        ss_trap, res_sums_trap, res_xis_trap = GaPSE.readxyall(
            "datatest/GNC_SumXiMultipoles/xis_GNC_L$L" * "_trap_noF_withobs_specific_ss.txt", comments=true)

        name = "calc_xis_GNC_L$L" * "_trap_noF_withobs_specific_ss.txt"
        isfile(name) && rm(name)

        GaPSE.print_map_sum_ξ_GNC_multipole(COSMO, name, ss_trap;
            single=false, L=L, alg=:trap, kwargs...)
        calc_ss_trap, calc_sums_trap = GaPSE.readxy(name, comments=true)
        calc_xis_trap = [
            begin
                a_name = "all_GNC_standalones_CFs/xi_GNC_" * effect * "_L$L" * ".txt"
                _, xis = GaPSE.readxy(a_name, comments=true)
                xis
            end for effect in GaPSE.GR_EFFECTS_GNC
        ]

        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(ss_trap, calc_ss_trap)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_sums_trap, calc_sums_trap)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[1], calc_xis_trap[1])]) # auto_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[2], calc_xis_trap[2])]) # auto_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[3], calc_xis_trap[3])]) # auto_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[4], calc_xis_trap[4])]) # auto_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[5], calc_xis_trap[5])]) # auto_integrated
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[6], calc_xis_trap[6])]) # newton_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[7], calc_xis_trap[7])]) # doppler_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[8], calc_xis_trap[8])]) # newton_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[9], calc_xis_trap[9])]) # lensing_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[10], calc_xis_trap[10])]) # newton_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[11], calc_xis_trap[11])]) # localgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[12], calc_xis_trap[12])]) # newton_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[13], calc_xis_trap[13])]) # integratedgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[14], calc_xis_trap[14])]) # lensing_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[15], calc_xis_trap[15])]) # doppler_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[16], calc_xis_trap[16])]) # doppler_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[17], calc_xis_trap[17])]) # localgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[18], calc_xis_trap[18])]) # doppler_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[19], calc_xis_trap[19])]) # integratedgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[20], calc_xis_trap[20])]) # lensing_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[21], calc_xis_trap[21])]) # localgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[22], calc_xis_trap[22])]) # lensing_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[23], calc_xis_trap[23])]) # integratedgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[24], calc_xis_trap[24])]) # localgp_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[25], calc_xis_trap[25])]) # integratedgp_localgp

        rm(name)
        for effect in GaPSE.GR_EFFECTS_GNC
            rm("all_GNC_standalones_CFs/xi_GNC_" * effect * "_L$L" * ".txt")
        end
        rm("all_GNC_standalones_CFs")

        #### lob ####
        ss_lob, res_sums_lob, res_xis_lob = GaPSE.readxyall(
            "datatest/GNC_SumXiMultipoles/xis_GNC_L$L" * "_lob_noF_withobs_specific_ss.txt", comments=true)

        name = "calc_xis_GNC_L$L" * "_lob_noF_withobs_specific_ss.txt"
        isfile(name) && rm(name)

        GaPSE.print_map_sum_ξ_GNC_multipole(COSMO, name, ss_lob;
            single=false, L=L, alg=:lobatto, kwargs...)
        calc_ss_lob, calc_sums_lob = GaPSE.readxy(name, comments=true)
        calc_xis_lob = [
            begin
                a_name = "all_GNC_standalones_CFs/xi_GNC_" * effect * "_L$L" * ".txt"
                _, xis = GaPSE.readxy(a_name, comments=true)
                xis
            end for effect in GaPSE.GR_EFFECTS_GNC
        ]

        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(ss_lob, calc_ss_lob)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_sums_lob, calc_sums_lob)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[1], calc_xis_lob[1])]) # auto_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[2], calc_xis_lob[2])]) # auto_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[3], calc_xis_lob[3])]) # auto_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[4], calc_xis_lob[4])]) # auto_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[5], calc_xis_lob[5])]) # auto_integrated
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[6], calc_xis_lob[6])]) # newton_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[7], calc_xis_lob[7])]) # doppler_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[8], calc_xis_lob[8])]) # newton_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[9], calc_xis_lob[9])]) # lensing_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[10], calc_xis_lob[10])]) # newton_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[11], calc_xis_lob[11])]) # localgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[12], calc_xis_lob[12])]) # newton_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[13], calc_xis_lob[13])]) # integratedgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[14], calc_xis_lob[14])]) # lensing_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[15], calc_xis_lob[15])]) # doppler_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[16], calc_xis_lob[16])]) # doppler_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[17], calc_xis_lob[17])]) # localgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[18], calc_xis_lob[18])]) # doppler_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[19], calc_xis_lob[19])]) # integratedgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[20], calc_xis_lob[20])]) # lensing_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[21], calc_xis_lob[21])]) # localgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[22], calc_xis_lob[22])]) # lensing_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[23], calc_xis_lob[23])]) # integratedgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[24], calc_xis_lob[24])]) # localgp_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[25], calc_xis_lob[25])]) # integratedgp_localgp

        rm(name)
        for effect in GaPSE.GR_EFFECTS_GNC
            rm("all_GNC_standalones_CFs/xi_GNC_" * effect * "_L$L" * ".txt")
        end
        rm("all_GNC_standalones_CFs")

    end

    @testset "L = 2 no_window with obs" begin
        L = 2

        ##### quad ####
        ss_quad, res_sums_quad, res_xis_quad = GaPSE.readxyall(
            "datatest/GNC_SumXiMultipoles/xis_GNC_L$L" * "_quad_noF_withobs_specific_ss.txt", comments=true)

        name = "calc_xis_GNC_L$L" * "_quad_noF_withobs_specific_ss.txt"
        isfile(name) && rm(name)

        GaPSE.print_map_sum_ξ_GNC_multipole(COSMO, name, ss_quad;
            single=false, L=L, alg=:quad, kwargs...)
        calc_ss_quad, calc_sums_quad = GaPSE.readxy(name, comments=true)
        calc_xis_quad = [
            begin
                a_name = "all_GNC_standalones_CFs/xi_GNC_" * effect * "_L$L" * ".txt"
                _, xis = GaPSE.readxy(a_name, comments=true)
                xis
            end for effect in GaPSE.GR_EFFECTS_GNC
        ]

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
        for effect in GaPSE.GR_EFFECTS_GNC
            rm("all_GNC_standalones_CFs/xi_GNC_" * effect * "_L$L" * ".txt")
        end
        rm("all_GNC_standalones_CFs")

        #### trap ####
        ss_trap, res_sums_trap, res_xis_trap = GaPSE.readxyall(
            "datatest/GNC_SumXiMultipoles/xis_GNC_L$L" * "_trap_noF_withobs_specific_ss.txt", comments=true)

        name = "calc_xis_GNC_L$L" * "_trap_noF_withobs_specific_ss.txt"
        isfile(name) && rm(name)

        GaPSE.print_map_sum_ξ_GNC_multipole(COSMO, name, ss_trap;
            single=false, L=L, alg=:trap, kwargs...)
        calc_ss_trap, calc_sums_trap = GaPSE.readxy(name, comments=true)
        calc_xis_trap = [
            begin
                a_name = "all_GNC_standalones_CFs/xi_GNC_" * effect * "_L$L" * ".txt"
                _, xis = GaPSE.readxy(a_name, comments=true)
                xis
            end for effect in GaPSE.GR_EFFECTS_GNC
        ]

        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(ss_trap, calc_ss_trap)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_sums_trap, calc_sums_trap)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[1], calc_xis_trap[1])]) # auto_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[2], calc_xis_trap[2])]) # auto_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[3], calc_xis_trap[3])]) # auto_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[4], calc_xis_trap[4])]) # auto_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[5], calc_xis_trap[5])]) # auto_integrated
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[6], calc_xis_trap[6])]) # newton_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[7], calc_xis_trap[7])]) # doppler_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[8], calc_xis_trap[8])]) # newton_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[9], calc_xis_trap[9])]) # lensing_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[10], calc_xis_trap[10])]) # newton_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[11], calc_xis_trap[11])]) # localgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[12], calc_xis_trap[12])]) # newton_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[13], calc_xis_trap[13])]) # integratedgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[14], calc_xis_trap[14])]) # lensing_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[15], calc_xis_trap[15])]) # doppler_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[16], calc_xis_trap[16])]) # doppler_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[17], calc_xis_trap[17])]) # localgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[18], calc_xis_trap[18])]) # doppler_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[19], calc_xis_trap[19])]) # integratedgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[20], calc_xis_trap[20])]) # lensing_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[21], calc_xis_trap[21])]) # localgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[22], calc_xis_trap[22])]) # lensing_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[23], calc_xis_trap[23])]) # integratedgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[24], calc_xis_trap[24])]) # localgp_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[25], calc_xis_trap[25])]) # integratedgp_localgp

        rm(name)
        for effect in GaPSE.GR_EFFECTS_GNC
            rm("all_GNC_standalones_CFs/xi_GNC_" * effect * "_L$L" * ".txt")
        end
        rm("all_GNC_standalones_CFs")

        #### lob ####
        ss_lob, res_sums_lob, res_xis_lob = GaPSE.readxyall(
            "datatest/GNC_SumXiMultipoles/xis_GNC_L$L" * "_lob_noF_withobs_specific_ss.txt", comments=true)

        name = "calc_xis_GNC_L$L" * "_lob_noF_withobs_specific_ss.txt"
        isfile(name) && rm(name)

        GaPSE.print_map_sum_ξ_GNC_multipole(COSMO, name, ss_lob;
            single=false, L=L, alg=:lobatto, kwargs...)
        calc_ss_lob, calc_sums_lob = GaPSE.readxy(name, comments=true)
        calc_xis_lob = [
            begin
                a_name = "all_GNC_standalones_CFs/xi_GNC_" * effect * "_L$L" * ".txt"
                _, xis = GaPSE.readxy(a_name, comments=true)
                xis
            end for effect in GaPSE.GR_EFFECTS_GNC
        ]

        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(ss_lob, calc_ss_lob)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_sums_lob, calc_sums_lob)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[1], calc_xis_lob[1])]) # auto_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[2], calc_xis_lob[2])]) # auto_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[3], calc_xis_lob[3])]) # auto_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[4], calc_xis_lob[4])]) # auto_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[5], calc_xis_lob[5])]) # auto_integrated
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[6], calc_xis_lob[6])]) # newton_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[7], calc_xis_lob[7])]) # doppler_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[8], calc_xis_lob[8])]) # newton_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[9], calc_xis_lob[9])]) # lensing_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[10], calc_xis_lob[10])]) # newton_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[11], calc_xis_lob[11])]) # localgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[12], calc_xis_lob[12])]) # newton_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[13], calc_xis_lob[13])]) # integratedgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[14], calc_xis_lob[14])]) # lensing_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[15], calc_xis_lob[15])]) # doppler_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[16], calc_xis_lob[16])]) # doppler_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[17], calc_xis_lob[17])]) # localgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[18], calc_xis_lob[18])]) # doppler_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[19], calc_xis_lob[19])]) # integratedgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[20], calc_xis_lob[20])]) # lensing_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[21], calc_xis_lob[21])]) # localgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[22], calc_xis_lob[22])]) # lensing_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[23], calc_xis_lob[23])]) # integratedgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[24], calc_xis_lob[24])]) # localgp_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[25], calc_xis_lob[25])]) # integratedgp_localgp

        rm(name)
        for effect in GaPSE.GR_EFFECTS_GNC
            rm("all_GNC_standalones_CFs/xi_GNC_" * effect * "_L$L" * ".txt")
        end
        rm("all_GNC_standalones_CFs")

    end

    @testset "L = 3 no_window with obs" begin
        L = 3

        #### trap ####
        ss_trap, res_sums_trap, res_xis_trap = GaPSE.readxyall(
            "datatest/GNC_SumXiMultipoles/xis_GNC_L$L" * "_trap_noF_withobs_specific_ss.txt", comments=true)

        name = "calc_xis_GNC_L$L" * "_trap_noF_withobs_specific_ss.txt"
        isfile(name) && rm(name)

        GaPSE.print_map_sum_ξ_GNC_multipole(COSMO, name, ss_trap;
            single=false, L=L, alg=:trap, kwargs...)
        calc_ss_trap, calc_sums_trap = GaPSE.readxy(name, comments=true)
        calc_xis_trap = [
            begin
                a_name = "all_GNC_standalones_CFs/xi_GNC_" * effect * "_L$L" * ".txt"
                _, xis = GaPSE.readxy(a_name, comments=true)
                xis
            end for effect in GaPSE.GR_EFFECTS_GNC
        ]

        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(ss_trap, calc_ss_trap)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_sums_trap, calc_sums_trap)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[1], calc_xis_trap[1])]) # auto_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[2], calc_xis_trap[2])]) # auto_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[3], calc_xis_trap[3])]) # auto_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[4], calc_xis_trap[4])]) # auto_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[5], calc_xis_trap[5])]) # auto_integrated
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[6], calc_xis_trap[6])]) # newton_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[7], calc_xis_trap[7])]) # doppler_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[8], calc_xis_trap[8])]) # newton_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[9], calc_xis_trap[9])]) # lensing_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[10], calc_xis_trap[10])]) # newton_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[11], calc_xis_trap[11])]) # localgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[12], calc_xis_trap[12])]) # newton_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[13], calc_xis_trap[13])]) # integratedgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[14], calc_xis_trap[14])]) # lensing_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[15], calc_xis_trap[15])]) # doppler_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[16], calc_xis_trap[16])]) # doppler_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[17], calc_xis_trap[17])]) # localgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[18], calc_xis_trap[18])]) # doppler_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[19], calc_xis_trap[19])]) # integratedgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[20], calc_xis_trap[20])]) # lensing_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[21], calc_xis_trap[21])]) # localgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[22], calc_xis_trap[22])]) # lensing_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[23], calc_xis_trap[23])]) # integratedgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[24], calc_xis_trap[24])]) # localgp_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[25], calc_xis_trap[25])]) # integratedgp_localgp

        rm(name)
        for effect in GaPSE.GR_EFFECTS_GNC
            rm("all_GNC_standalones_CFs/xi_GNC_" * effect * "_L$L" * ".txt")
        end
        rm("all_GNC_standalones_CFs")

        #### lob ####
        ss_lob, res_sums_lob, res_xis_lob = GaPSE.readxyall(
            "datatest/GNC_SumXiMultipoles/xis_GNC_L$L" * "_lob_noF_withobs_specific_ss.txt", comments=true)

        name = "calc_xis_GNC_L$L" * "_lob_noF_withobs_specific_ss.txt"
        isfile(name) && rm(name)

        GaPSE.print_map_sum_ξ_GNC_multipole(COSMO, name, ss_lob;
            single=false, L=L, alg=:lobatto, kwargs...)
        calc_ss_lob, calc_sums_lob = GaPSE.readxy(name, comments=true)
        calc_xis_lob = [
            begin
                a_name = "all_GNC_standalones_CFs/xi_GNC_" * effect * "_L$L" * ".txt"
                _, xis = GaPSE.readxy(a_name, comments=true)
                xis
            end for effect in GaPSE.GR_EFFECTS_GNC
        ]

        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(ss_lob, calc_ss_lob)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_sums_lob, calc_sums_lob)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[1], calc_xis_lob[1])]) # auto_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[2], calc_xis_lob[2])]) # auto_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[3], calc_xis_lob[3])]) # auto_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[4], calc_xis_lob[4])]) # auto_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[5], calc_xis_lob[5])]) # auto_integrated
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[6], calc_xis_lob[6])]) # newton_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[7], calc_xis_lob[7])]) # doppler_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[8], calc_xis_lob[8])]) # newton_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[9], calc_xis_lob[9])]) # lensing_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[10], calc_xis_lob[10])]) # newton_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[11], calc_xis_lob[11])]) # localgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[12], calc_xis_lob[12])]) # newton_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[13], calc_xis_lob[13])]) # integratedgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[14], calc_xis_lob[14])]) # lensing_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[15], calc_xis_lob[15])]) # doppler_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[16], calc_xis_lob[16])]) # doppler_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[17], calc_xis_lob[17])]) # localgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[18], calc_xis_lob[18])]) # doppler_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[19], calc_xis_lob[19])]) # integratedgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[20], calc_xis_lob[20])]) # lensing_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[21], calc_xis_lob[21])]) # localgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[22], calc_xis_lob[22])]) # lensing_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[23], calc_xis_lob[23])]) # integratedgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[24], calc_xis_lob[24])]) # localgp_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[25], calc_xis_lob[25])]) # integratedgp_localgp

        rm(name)
        for effect in GaPSE.GR_EFFECTS_GNC
            rm("all_GNC_standalones_CFs/xi_GNC_" * effect * "_L$L" * ".txt")
        end
        rm("all_GNC_standalones_CFs")

    end

    @testset "L = 4 no_window with obs" begin
        L = 4

        #### trap ####
        ss_trap, res_sums_trap, res_xis_trap = GaPSE.readxyall(
            "datatest/GNC_SumXiMultipoles/xis_GNC_L$L" * "_trap_noF_withobs_specific_ss.txt", comments=true)

        name = "calc_xis_GNC_L$L" * "_trap_noF_withobs_specific_ss.txt"
        isfile(name) && rm(name)

        GaPSE.print_map_sum_ξ_GNC_multipole(COSMO, name, ss_trap;
            single=false, L=L, alg=:trap, kwargs...)
        calc_ss_trap, calc_sums_trap = GaPSE.readxy(name, comments=true)
        calc_xis_trap = [
            begin
                a_name = "all_GNC_standalones_CFs/xi_GNC_" * effect * "_L$L" * ".txt"
                _, xis = GaPSE.readxy(a_name, comments=true)
                xis
            end for effect in GaPSE.GR_EFFECTS_GNC
        ]

        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(ss_trap, calc_ss_trap)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_sums_trap, calc_sums_trap)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[1], calc_xis_trap[1])]) # auto_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[2], calc_xis_trap[2])]) # auto_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[3], calc_xis_trap[3])]) # auto_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[4], calc_xis_trap[4])]) # auto_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[5], calc_xis_trap[5])]) # auto_integrated
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[6], calc_xis_trap[6])]) # newton_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[7], calc_xis_trap[7])]) # doppler_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[8], calc_xis_trap[8])]) # newton_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[9], calc_xis_trap[9])]) # lensing_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[10], calc_xis_trap[10])]) # newton_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[11], calc_xis_trap[11])]) # localgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[12], calc_xis_trap[12])]) # newton_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[13], calc_xis_trap[13])]) # integratedgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[14], calc_xis_trap[14])]) # lensing_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[15], calc_xis_trap[15])]) # doppler_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[16], calc_xis_trap[16])]) # doppler_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[17], calc_xis_trap[17])]) # localgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[18], calc_xis_trap[18])]) # doppler_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[19], calc_xis_trap[19])]) # integratedgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[20], calc_xis_trap[20])]) # lensing_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[21], calc_xis_trap[21])]) # localgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[22], calc_xis_trap[22])]) # lensing_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[23], calc_xis_trap[23])]) # integratedgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[24], calc_xis_trap[24])]) # localgp_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[25], calc_xis_trap[25])]) # integratedgp_localgp

        rm(name)
        for effect in GaPSE.GR_EFFECTS_GNC
            rm("all_GNC_standalones_CFs/xi_GNC_" * effect * "_L$L" * ".txt")
        end
        rm("all_GNC_standalones_CFs")

        #### lob ####
        ss_lob, res_sums_lob, res_xis_lob = GaPSE.readxyall(
            "datatest/GNC_SumXiMultipoles/xis_GNC_L$L" * "_lob_noF_withobs_specific_ss.txt", comments=true)

        name = "calc_xis_GNC_L$L" * "_lob_noF_withobs_specific_ss.txt"
        isfile(name) && rm(name)

        GaPSE.print_map_sum_ξ_GNC_multipole(COSMO, name, ss_lob;
            single=false, L=L, alg=:lobatto, kwargs...)
        calc_ss_lob, calc_sums_lob = GaPSE.readxy(name, comments=true)
        calc_xis_lob = [
            begin
                a_name = "all_GNC_standalones_CFs/xi_GNC_" * effect * "_L$L" * ".txt"
                _, xis = GaPSE.readxy(a_name, comments=true)
                xis
            end for effect in GaPSE.GR_EFFECTS_GNC
        ]

        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(ss_lob, calc_ss_lob)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_sums_lob, calc_sums_lob)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[1], calc_xis_lob[1])]) # auto_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[2], calc_xis_lob[2])]) # auto_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[3], calc_xis_lob[3])]) # auto_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[4], calc_xis_lob[4])]) # auto_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[5], calc_xis_lob[5])]) # auto_integrated
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[6], calc_xis_lob[6])]) # newton_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[7], calc_xis_lob[7])]) # doppler_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[8], calc_xis_lob[8])]) # newton_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[9], calc_xis_lob[9])]) # lensing_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[10], calc_xis_lob[10])]) # newton_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[11], calc_xis_lob[11])]) # localgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[12], calc_xis_lob[12])]) # newton_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[13], calc_xis_lob[13])]) # integratedgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[14], calc_xis_lob[14])]) # lensing_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[15], calc_xis_lob[15])]) # doppler_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[16], calc_xis_lob[16])]) # doppler_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[17], calc_xis_lob[17])]) # localgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[18], calc_xis_lob[18])]) # doppler_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[19], calc_xis_lob[19])]) # integratedgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[20], calc_xis_lob[20])]) # lensing_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[21], calc_xis_lob[21])]) # localgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[22], calc_xis_lob[22])]) # lensing_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[23], calc_xis_lob[23])]) # integratedgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[24], calc_xis_lob[24])]) # localgp_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[25], calc_xis_lob[25])]) # integratedgp_localgp

        rm(name)
        for effect in GaPSE.GR_EFFECTS_GNC
            rm("all_GNC_standalones_CFs/xi_GNC_" * effect * "_L$L" * ".txt")
        end
        rm("all_GNC_standalones_CFs")

    end
end

println("\n  7 / 8 ...")

@testset "test print_map_sum_ξ_GNC_multipole with_window no obs, second method" begin
    RTOL = 1.2e-2

    kwargs = Dict(
        :use_windows => true, :obs => :no,
        :enhancer => 1e8, :s1 => nothing,
        :N_trap => 30, :N_lob => 30,
        :atol_quad => 0.0, :rtol_quad => 1e-2,
        :N_χs => 40, :N_χs_2 => 20,
        :pr => false,
    )

    @testset "test_zeros" begin
        @test_throws AssertionError GaPSE.GaPSE.print_map_sum_ξ_GNC_multipole(COSMO, "NONEXISTINGDIR/file.txt", [1.0, 2.0, 3.0])
    end

    @testset "L = 0 with_window no obs" begin
        L = 0

        ##### quad ####
        ss_quad, res_sums_quad, res_xis_quad = GaPSE.readxyall(
            "datatest/GNC_SumXiMultipoles/xis_GNC_L$L" * "_quad_withF_noobs_specific_ss.txt", comments=true)

        name = "calc_xis_GNC_L$L" * "_quad_withF_noobs_specific_ss.txt"
        isfile(name) && rm(name)

        GaPSE.print_map_sum_ξ_GNC_multipole(COSMO, name, ss_quad;
            single=false, L=L, alg=:quad, kwargs...)

        calc_ss_quad, calc_sums_quad = GaPSE.readxy(name, comments=true)
        calc_xis_quad = [
            begin
                a_name = "all_GNC_standalones_CFs/xi_GNC_" * effect * "_L$L" * ".txt"
                _, xis = GaPSE.readxy(a_name, comments=true)
                xis
            end for effect in GaPSE.GR_EFFECTS_GNC
        ]

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
        for effect in GaPSE.GR_EFFECTS_GNC
            rm("all_GNC_standalones_CFs/xi_GNC_" * effect * "_L$L" * ".txt")
        end
        rm("all_GNC_standalones_CFs")

        #### trap ####
        ss_trap, res_sums_trap, res_xis_trap = GaPSE.readxyall(
            "datatest/GNC_SumXiMultipoles/xis_GNC_L$L" * "_trap_withF_noobs_specific_ss.txt", comments=true)

        name = "calc_xis_GNC_L$L" * "_trap_withF_noobs_specific_ss.txt"
        isfile(name) && rm(name)

        GaPSE.print_map_sum_ξ_GNC_multipole(COSMO, name, ss_trap;
            single=false, L=L, alg=:trap, kwargs...)
        calc_ss_trap, calc_sums_trap = GaPSE.readxy(name, comments=true)
        calc_xis_trap = [
            begin
                a_name = "all_GNC_standalones_CFs/xi_GNC_" * effect * "_L$L" * ".txt"
                _, xis = GaPSE.readxy(a_name, comments=true)
                xis
            end for effect in GaPSE.GR_EFFECTS_GNC
        ]

        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(ss_trap, calc_ss_trap)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_sums_trap, calc_sums_trap)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[1], calc_xis_trap[1])]) # auto_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[2], calc_xis_trap[2])]) # auto_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[3], calc_xis_trap[3])]) # auto_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[4], calc_xis_trap[4])]) # auto_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[5], calc_xis_trap[5])]) # auto_integrated
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[6], calc_xis_trap[6])]) # newton_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[7], calc_xis_trap[7])]) # doppler_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[8], calc_xis_trap[8])]) # newton_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[9], calc_xis_trap[9])]) # lensing_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[10], calc_xis_trap[10])]) # newton_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[11], calc_xis_trap[11])]) # localgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[12], calc_xis_trap[12])]) # newton_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[13], calc_xis_trap[13])]) # integratedgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[14], calc_xis_trap[14])]) # lensing_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[15], calc_xis_trap[15])]) # doppler_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[16], calc_xis_trap[16])]) # doppler_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[17], calc_xis_trap[17])]) # localgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[18], calc_xis_trap[18])]) # doppler_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[19], calc_xis_trap[19])]) # integratedgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[20], calc_xis_trap[20])]) # lensing_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[21], calc_xis_trap[21])]) # localgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[22], calc_xis_trap[22])]) # lensing_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[23], calc_xis_trap[23])]) # integratedgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[24], calc_xis_trap[24])]) # localgp_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[25], calc_xis_trap[25])]) # integratedgp_localgp

        rm(name)
        for effect in GaPSE.GR_EFFECTS_GNC
            rm("all_GNC_standalones_CFs/xi_GNC_" * effect * "_L$L" * ".txt")
        end
        rm("all_GNC_standalones_CFs")

        #### lob ####
        ss_lob, res_sums_lob, res_xis_lob = GaPSE.readxyall(
            "datatest/GNC_SumXiMultipoles/xis_GNC_L$L" * "_lob_withF_noobs_specific_ss.txt", comments=true)

        name = "calc_xis_GNC_L$L" * "_lob_withF_noobs_specific_ss.txt"
        isfile(name) && rm(name)

        GaPSE.print_map_sum_ξ_GNC_multipole(COSMO, name, ss_lob;
            single=false, L=L, alg=:lobatto, kwargs...)
        calc_ss_lob, calc_sums_lob, calc_xis_lob = GaPSE.readxyall(name, comments=true)
        calc_ss_lob, calc_sums_lob = GaPSE.readxy(name, comments=true)
        calc_xis_lob = [
            begin
                a_name = "all_GNC_standalones_CFs/xi_GNC_" * effect * "_L$L" * ".txt"
                _, xis = GaPSE.readxy(a_name, comments=true)
                xis
            end for effect in GaPSE.GR_EFFECTS_GNC
        ]

        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(ss_lob, calc_ss_lob)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_sums_lob, calc_sums_lob)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[1], calc_xis_lob[1])]) # auto_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[2], calc_xis_lob[2])]) # auto_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[3], calc_xis_lob[3])]) # auto_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[4], calc_xis_lob[4])]) # auto_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[5], calc_xis_lob[5])]) # auto_integrated
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[6], calc_xis_lob[6])]) # newton_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[7], calc_xis_lob[7])]) # doppler_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[8], calc_xis_lob[8])]) # newton_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[9], calc_xis_lob[9])]) # lensing_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[10], calc_xis_lob[10])]) # newton_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[11], calc_xis_lob[11])]) # localgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[12], calc_xis_lob[12])]) # newton_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[13], calc_xis_lob[13])]) # integratedgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[14], calc_xis_lob[14])]) # lensing_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[15], calc_xis_lob[15])]) # doppler_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[16], calc_xis_lob[16])]) # doppler_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[17], calc_xis_lob[17])]) # localgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[18], calc_xis_lob[18])]) # doppler_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[19], calc_xis_lob[19])]) # integratedgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[20], calc_xis_lob[20])]) # lensing_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[21], calc_xis_lob[21])]) # localgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[22], calc_xis_lob[22])]) # lensing_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[23], calc_xis_lob[23])]) # integratedgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[24], calc_xis_lob[24])]) # localgp_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[25], calc_xis_lob[25])]) # integratedgp_localgp

        rm(name)
        for effect in GaPSE.GR_EFFECTS_GNC
            rm("all_GNC_standalones_CFs/xi_GNC_" * effect * "_L$L" * ".txt")
        end
        rm("all_GNC_standalones_CFs")

    end

    @testset "L = 1 with_window no obs" begin
        L = 1

        #### trap ####
        ss_trap, res_sums_trap, res_xis_trap = GaPSE.readxyall(
            "datatest/GNC_SumXiMultipoles/xis_GNC_L$L" * "_trap_withF_noobs_specific_ss.txt", comments=true)

        name = "calc_xis_GNC_L$L" * "_trap_withF_noobs_specific_ss.txt"
        isfile(name) && rm(name)

        GaPSE.print_map_sum_ξ_GNC_multipole(COSMO, name, ss_trap;
            single=false, L=L, alg=:trap, kwargs...)
        calc_ss_trap, calc_sums_trap = GaPSE.readxy(name, comments=true)
        calc_xis_trap = [
            begin
                a_name = "all_GNC_standalones_CFs/xi_GNC_" * effect * "_L$L" * ".txt"
                _, xis = GaPSE.readxy(a_name, comments=true)
                xis
            end for effect in GaPSE.GR_EFFECTS_GNC
        ]

        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(ss_trap, calc_ss_trap)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_sums_trap, calc_sums_trap)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[1], calc_xis_trap[1])]) # auto_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[2], calc_xis_trap[2])]) # auto_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[3], calc_xis_trap[3])]) # auto_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[4], calc_xis_trap[4])]) # auto_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[5], calc_xis_trap[5])]) # auto_integrated
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[6], calc_xis_trap[6])]) # newton_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[7], calc_xis_trap[7])]) # doppler_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[8], calc_xis_trap[8])]) # newton_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[9], calc_xis_trap[9])]) # lensing_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[10], calc_xis_trap[10])]) # newton_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[11], calc_xis_trap[11])]) # localgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[12], calc_xis_trap[12])]) # newton_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[13], calc_xis_trap[13])]) # integratedgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[14], calc_xis_trap[14])]) # lensing_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[15], calc_xis_trap[15])]) # doppler_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[16], calc_xis_trap[16])]) # doppler_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[17], calc_xis_trap[17])]) # localgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[18], calc_xis_trap[18])]) # doppler_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[19], calc_xis_trap[19])]) # integratedgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[20], calc_xis_trap[20])]) # lensing_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[21], calc_xis_trap[21])]) # localgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[22], calc_xis_trap[22])]) # lensing_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[23], calc_xis_trap[23])]) # integratedgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[24], calc_xis_trap[24])]) # localgp_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[25], calc_xis_trap[25])]) # integratedgp_localgp

        rm(name)
        for effect in GaPSE.GR_EFFECTS_GNC
            rm("all_GNC_standalones_CFs/xi_GNC_" * effect * "_L$L" * ".txt")
        end
        rm("all_GNC_standalones_CFs")

        #### lob ####
        ss_lob, res_sums_lob, res_xis_lob = GaPSE.readxyall(
            "datatest/GNC_SumXiMultipoles/xis_GNC_L$L" * "_lob_withF_noobs_specific_ss.txt", comments=true)

        name = "calc_xis_GNC_L$L" * "_lob_withF_noobs_specific_ss.txt"
        isfile(name) && rm(name)

        GaPSE.print_map_sum_ξ_GNC_multipole(COSMO, name, ss_lob;
            single=false, L=L, alg=:lobatto, kwargs...)
        calc_ss_lob, calc_sums_lob = GaPSE.readxy(name, comments=true)
        calc_xis_lob = [
            begin
                a_name = "all_GNC_standalones_CFs/xi_GNC_" * effect * "_L$L" * ".txt"
                _, xis = GaPSE.readxy(a_name, comments=true)
                xis
            end for effect in GaPSE.GR_EFFECTS_GNC
        ]

        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(ss_lob, calc_ss_lob)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_sums_lob, calc_sums_lob)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[1], calc_xis_lob[1])]) # auto_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[2], calc_xis_lob[2])]) # auto_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[3], calc_xis_lob[3])]) # auto_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[4], calc_xis_lob[4])]) # auto_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[5], calc_xis_lob[5])]) # auto_integrated
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[6], calc_xis_lob[6])]) # newton_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[7], calc_xis_lob[7])]) # doppler_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[8], calc_xis_lob[8])]) # newton_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[9], calc_xis_lob[9])]) # lensing_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[10], calc_xis_lob[10])]) # newton_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[11], calc_xis_lob[11])]) # localgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[12], calc_xis_lob[12])]) # newton_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[13], calc_xis_lob[13])]) # integratedgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[14], calc_xis_lob[14])]) # lensing_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[15], calc_xis_lob[15])]) # doppler_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[16], calc_xis_lob[16])]) # doppler_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[17], calc_xis_lob[17])]) # localgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[18], calc_xis_lob[18])]) # doppler_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[19], calc_xis_lob[19])]) # integratedgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[20], calc_xis_lob[20])]) # lensing_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[21], calc_xis_lob[21])]) # localgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[22], calc_xis_lob[22])]) # lensing_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[23], calc_xis_lob[23])]) # integratedgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[24], calc_xis_lob[24])]) # localgp_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[25], calc_xis_lob[25])]) # integratedgp_localgp

        rm(name)
        for effect in GaPSE.GR_EFFECTS_GNC
            rm("all_GNC_standalones_CFs/xi_GNC_" * effect * "_L$L" * ".txt")
        end
        rm("all_GNC_standalones_CFs")

    end

    @testset "L = 2 with_window no obs" begin
        L = 2

        ##### quad ####
        ss_quad, res_sums_quad, res_xis_quad = GaPSE.readxyall(
            "datatest/GNC_SumXiMultipoles/xis_GNC_L$L" * "_quad_withF_noobs_specific_ss.txt", comments=true)

        name = "calc_xis_GNC_L$L" * "_quad_withF_noobs_specific_ss.txt"
        isfile(name) && rm(name)

        GaPSE.print_map_sum_ξ_GNC_multipole(COSMO, name, ss_quad;
            single=false, L=L, alg=:quad, kwargs...)
        calc_ss_quad, calc_sums_quad = GaPSE.readxy(name, comments=true)
        calc_xis_quad = [
            begin
                a_name = "all_GNC_standalones_CFs/xi_GNC_" * effect * "_L$L" * ".txt"
                _, xis = GaPSE.readxy(a_name, comments=true)
                xis
            end for effect in GaPSE.GR_EFFECTS_GNC
        ]

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
        for effect in GaPSE.GR_EFFECTS_GNC
            rm("all_GNC_standalones_CFs/xi_GNC_" * effect * "_L$L" * ".txt")
        end
        rm("all_GNC_standalones_CFs")

        #### trap ####
        ss_trap, res_sums_trap, res_xis_trap = GaPSE.readxyall(
            "datatest/GNC_SumXiMultipoles/xis_GNC_L$L" * "_trap_withF_noobs_specific_ss.txt", comments=true)

        name = "calc_xis_GNC_L$L" * "_trap_withF_noobs_specific_ss.txt"
        isfile(name) && rm(name)

        GaPSE.print_map_sum_ξ_GNC_multipole(COSMO, name, ss_trap;
            single=false, L=L, alg=:trap, kwargs...)
        calc_ss_trap, calc_sums_trap = GaPSE.readxy(name, comments=true)
        calc_xis_trap = [
            begin
                a_name = "all_GNC_standalones_CFs/xi_GNC_" * effect * "_L$L" * ".txt"
                _, xis = GaPSE.readxy(a_name, comments=true)
                xis
            end for effect in GaPSE.GR_EFFECTS_GNC
        ]

        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(ss_trap, calc_ss_trap)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_sums_trap, calc_sums_trap)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[1], calc_xis_trap[1])]) # auto_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[2], calc_xis_trap[2])]) # auto_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[3], calc_xis_trap[3])]) # auto_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[4], calc_xis_trap[4])]) # auto_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[5], calc_xis_trap[5])]) # auto_integrated
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[6], calc_xis_trap[6])]) # newton_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[7], calc_xis_trap[7])]) # doppler_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[8], calc_xis_trap[8])]) # newton_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[9], calc_xis_trap[9])]) # lensing_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[10], calc_xis_trap[10])]) # newton_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[11], calc_xis_trap[11])]) # localgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[12], calc_xis_trap[12])]) # newton_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[13], calc_xis_trap[13])]) # integratedgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[14], calc_xis_trap[14])]) # lensing_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[15], calc_xis_trap[15])]) # doppler_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[16], calc_xis_trap[16])]) # doppler_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[17], calc_xis_trap[17])]) # localgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[18], calc_xis_trap[18])]) # doppler_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[19], calc_xis_trap[19])]) # integratedgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[20], calc_xis_trap[20])]) # lensing_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[21], calc_xis_trap[21])]) # localgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[22], calc_xis_trap[22])]) # lensing_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[23], calc_xis_trap[23])]) # integratedgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[24], calc_xis_trap[24])]) # localgp_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[25], calc_xis_trap[25])]) # integratedgp_localgp

        rm(name)
        for effect in GaPSE.GR_EFFECTS_GNC
            rm("all_GNC_standalones_CFs/xi_GNC_" * effect * "_L$L" * ".txt")
        end
        rm("all_GNC_standalones_CFs")

        #### lob ####
        ss_lob, res_sums_lob, res_xis_lob = GaPSE.readxyall(
            "datatest/GNC_SumXiMultipoles/xis_GNC_L$L" * "_lob_withF_noobs_specific_ss.txt", comments=true)

        name = "calc_xis_GNC_L$L" * "_lob_withF_noobs_specific_ss.txt"
        isfile(name) && rm(name)

        GaPSE.print_map_sum_ξ_GNC_multipole(COSMO, name, ss_lob;
            single=false, L=L, alg=:lobatto, kwargs...)
        calc_ss_lob, calc_sums_lob = GaPSE.readxy(name, comments=true)
        calc_xis_lob = [
            begin
                a_name = "all_GNC_standalones_CFs/xi_GNC_" * effect * "_L$L" * ".txt"
                _, xis = GaPSE.readxy(a_name, comments=true)
                xis
            end for effect in GaPSE.GR_EFFECTS_GNC
        ]

        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(ss_lob, calc_ss_lob)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_sums_lob, calc_sums_lob)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[1], calc_xis_lob[1])]) # auto_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[2], calc_xis_lob[2])]) # auto_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[3], calc_xis_lob[3])]) # auto_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[4], calc_xis_lob[4])]) # auto_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[5], calc_xis_lob[5])]) # auto_integrated
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[6], calc_xis_lob[6])]) # newton_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[7], calc_xis_lob[7])]) # doppler_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[8], calc_xis_lob[8])]) # newton_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[9], calc_xis_lob[9])]) # lensing_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[10], calc_xis_lob[10])]) # newton_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[11], calc_xis_lob[11])]) # localgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[12], calc_xis_lob[12])]) # newton_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[13], calc_xis_lob[13])]) # integratedgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[14], calc_xis_lob[14])]) # lensing_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[15], calc_xis_lob[15])]) # doppler_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[16], calc_xis_lob[16])]) # doppler_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[17], calc_xis_lob[17])]) # localgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[18], calc_xis_lob[18])]) # doppler_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[19], calc_xis_lob[19])]) # integratedgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[20], calc_xis_lob[20])]) # lensing_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[21], calc_xis_lob[21])]) # localgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[22], calc_xis_lob[22])]) # lensing_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[23], calc_xis_lob[23])]) # integratedgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[24], calc_xis_lob[24])]) # localgp_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[25], calc_xis_lob[25])]) # integratedgp_localgp

        rm(name)
        for effect in GaPSE.GR_EFFECTS_GNC
            rm("all_GNC_standalones_CFs/xi_GNC_" * effect * "_L$L" * ".txt")
        end
        rm("all_GNC_standalones_CFs")

    end

    @testset "L = 3 with_window no obs" begin
        L = 3
        
        #### trap ####
        ss_trap, res_sums_trap, res_xis_trap = GaPSE.readxyall(
            "datatest/GNC_SumXiMultipoles/xis_GNC_L$L" * "_trap_withF_noobs_specific_ss.txt", comments=true)

        name = "calc_xis_GNC_L$L" * "_trap_withF_noobs_specific_ss.txt"
        isfile(name) && rm(name)

        GaPSE.print_map_sum_ξ_GNC_multipole(COSMO, name, ss_trap;
            single=false, L=L, alg=:trap, kwargs...)
        calc_ss_trap, calc_sums_trap = GaPSE.readxy(name, comments=true)
        calc_xis_trap = [
            begin
                a_name = "all_GNC_standalones_CFs/xi_GNC_" * effect * "_L$L" * ".txt"
                _, xis = GaPSE.readxy(a_name, comments=true)
                xis
            end for effect in GaPSE.GR_EFFECTS_GNC
        ]

        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(ss_trap, calc_ss_trap)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_sums_trap, calc_sums_trap)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[1], calc_xis_trap[1])]) # auto_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[2], calc_xis_trap[2])]) # auto_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[3], calc_xis_trap[3])]) # auto_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[4], calc_xis_trap[4])]) # auto_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[5], calc_xis_trap[5])]) # auto_integrated
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[6], calc_xis_trap[6])]) # newton_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[7], calc_xis_trap[7])]) # doppler_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[8], calc_xis_trap[8])]) # newton_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[9], calc_xis_trap[9])]) # lensing_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[10], calc_xis_trap[10])]) # newton_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[11], calc_xis_trap[11])]) # localgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[12], calc_xis_trap[12])]) # newton_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[13], calc_xis_trap[13])]) # integratedgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[14], calc_xis_trap[14])]) # lensing_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[15], calc_xis_trap[15])]) # doppler_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[16], calc_xis_trap[16])]) # doppler_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[17], calc_xis_trap[17])]) # localgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[18], calc_xis_trap[18])]) # doppler_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[19], calc_xis_trap[19])]) # integratedgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[20], calc_xis_trap[20])]) # lensing_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[21], calc_xis_trap[21])]) # localgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[22], calc_xis_trap[22])]) # lensing_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[23], calc_xis_trap[23])]) # integratedgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[24], calc_xis_trap[24])]) # localgp_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[25], calc_xis_trap[25])]) # integratedgp_localgp

        rm(name)
        for effect in GaPSE.GR_EFFECTS_GNC
            rm("all_GNC_standalones_CFs/xi_GNC_" * effect * "_L$L" * ".txt")
        end
        rm("all_GNC_standalones_CFs")

        #### lob ####
        ss_lob, res_sums_lob, res_xis_lob = GaPSE.readxyall(
            "datatest/GNC_SumXiMultipoles/xis_GNC_L$L" * "_lob_withF_noobs_specific_ss.txt", comments=true)

        name = "calc_xis_GNC_L$L" * "_lob_withF_noobs_specific_ss.txt"
        isfile(name) && rm(name)

        GaPSE.print_map_sum_ξ_GNC_multipole(COSMO, name, ss_lob;
            single=false, L=L, alg=:lobatto, kwargs...)
        calc_ss_lob, calc_sums_lob = GaPSE.readxy(name, comments=true)
        calc_xis_lob = [
            begin
                a_name = "all_GNC_standalones_CFs/xi_GNC_" * effect * "_L$L" * ".txt"
                _, xis = GaPSE.readxy(a_name, comments=true)
                xis
            end for effect in GaPSE.GR_EFFECTS_GNC
        ]

        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(ss_lob, calc_ss_lob)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_sums_lob, calc_sums_lob)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[1], calc_xis_lob[1])]) # auto_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[2], calc_xis_lob[2])]) # auto_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[3], calc_xis_lob[3])]) # auto_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[4], calc_xis_lob[4])]) # auto_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[5], calc_xis_lob[5])]) # auto_integrated
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[6], calc_xis_lob[6])]) # newton_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[7], calc_xis_lob[7])]) # doppler_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[8], calc_xis_lob[8])]) # newton_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[9], calc_xis_lob[9])]) # lensing_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[10], calc_xis_lob[10])]) # newton_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[11], calc_xis_lob[11])]) # localgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[12], calc_xis_lob[12])]) # newton_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[13], calc_xis_lob[13])]) # integratedgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[14], calc_xis_lob[14])]) # lensing_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[15], calc_xis_lob[15])]) # doppler_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[16], calc_xis_lob[16])]) # doppler_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[17], calc_xis_lob[17])]) # localgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[18], calc_xis_lob[18])]) # doppler_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[19], calc_xis_lob[19])]) # integratedgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[20], calc_xis_lob[20])]) # lensing_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[21], calc_xis_lob[21])]) # localgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[22], calc_xis_lob[22])]) # lensing_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[23], calc_xis_lob[23])]) # integratedgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[24], calc_xis_lob[24])]) # localgp_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[25], calc_xis_lob[25])]) # integratedgp_localgp

        rm(name)
        for effect in GaPSE.GR_EFFECTS_GNC
            rm("all_GNC_standalones_CFs/xi_GNC_" * effect * "_L$L" * ".txt")
        end
        rm("all_GNC_standalones_CFs")

    end

    @testset "L = 4 with_window no obs" begin
        L = 4

        #### trap ####
        ss_trap, res_sums_trap, res_xis_trap = GaPSE.readxyall(
            "datatest/GNC_SumXiMultipoles/xis_GNC_L$L" * "_trap_withF_noobs_specific_ss.txt", comments=true)

        name = "calc_xis_GNC_L$L" * "_trap_withF_noobs_specific_ss.txt"
        isfile(name) && rm(name)

        GaPSE.print_map_sum_ξ_GNC_multipole(COSMO, name, ss_trap;
            single=false, L=L, alg=:trap, kwargs...)
        calc_ss_trap, calc_sums_trap = GaPSE.readxy(name, comments=true)
        calc_xis_trap = [
            begin
                a_name = "all_GNC_standalones_CFs/xi_GNC_" * effect * "_L$L" * ".txt"
                _, xis = GaPSE.readxy(a_name, comments=true)
                xis
            end for effect in GaPSE.GR_EFFECTS_GNC
        ]

        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(ss_trap, calc_ss_trap)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_sums_trap, calc_sums_trap)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[1], calc_xis_trap[1])]) # auto_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[2], calc_xis_trap[2])]) # auto_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[3], calc_xis_trap[3])]) # auto_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[4], calc_xis_trap[4])]) # auto_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[5], calc_xis_trap[5])]) # auto_integrated
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[6], calc_xis_trap[6])]) # newton_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[7], calc_xis_trap[7])]) # doppler_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[8], calc_xis_trap[8])]) # newton_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[9], calc_xis_trap[9])]) # lensing_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[10], calc_xis_trap[10])]) # newton_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[11], calc_xis_trap[11])]) # localgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[12], calc_xis_trap[12])]) # newton_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[13], calc_xis_trap[13])]) # integratedgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[14], calc_xis_trap[14])]) # lensing_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[15], calc_xis_trap[15])]) # doppler_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[16], calc_xis_trap[16])]) # doppler_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[17], calc_xis_trap[17])]) # localgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[18], calc_xis_trap[18])]) # doppler_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[19], calc_xis_trap[19])]) # integratedgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[20], calc_xis_trap[20])]) # lensing_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[21], calc_xis_trap[21])]) # localgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[22], calc_xis_trap[22])]) # lensing_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[23], calc_xis_trap[23])]) # integratedgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[24], calc_xis_trap[24])]) # localgp_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[25], calc_xis_trap[25])]) # integratedgp_localgp

        rm(name)
        for effect in GaPSE.GR_EFFECTS_GNC
            rm("all_GNC_standalones_CFs/xi_GNC_" * effect * "_L$L" * ".txt")
        end
        rm("all_GNC_standalones_CFs")

        #### lob ####
        ss_lob, res_sums_lob, res_xis_lob = GaPSE.readxyall(
            "datatest/GNC_SumXiMultipoles/xis_GNC_L$L" * "_lob_withF_noobs_specific_ss.txt", comments=true)

        name = "calc_xis_GNC_L$L" * "_lob_withF_noobs_specific_ss.txt"
        isfile(name) && rm(name)

        GaPSE.print_map_sum_ξ_GNC_multipole(COSMO, name, ss_lob;
            single=false, L=L, alg=:lobatto, kwargs...)
        calc_ss_lob, calc_sums_lob = GaPSE.readxy(name, comments=true)
        calc_xis_lob = [
            begin
                a_name = "all_GNC_standalones_CFs/xi_GNC_" * effect * "_L$L" * ".txt"
                _, xis = GaPSE.readxy(a_name, comments=true)
                xis
            end for effect in GaPSE.GR_EFFECTS_GNC
        ]

        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(ss_lob, calc_ss_lob)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_sums_lob, calc_sums_lob)])
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[1], calc_xis_lob[1])]) # auto_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[2], calc_xis_lob[2])]) # auto_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[3], calc_xis_lob[3])]) # auto_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[4], calc_xis_lob[4])]) # auto_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[5], calc_xis_lob[5])]) # auto_integrated
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[6], calc_xis_lob[6])]) # newton_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[7], calc_xis_lob[7])]) # doppler_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[8], calc_xis_lob[8])]) # newton_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[9], calc_xis_lob[9])]) # lensing_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[10], calc_xis_lob[10])]) # newton_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[11], calc_xis_lob[11])]) # localgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[12], calc_xis_lob[12])]) # newton_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[13], calc_xis_lob[13])]) # integratedgp_newton
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[14], calc_xis_lob[14])]) # lensing_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[15], calc_xis_lob[15])]) # doppler_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[16], calc_xis_lob[16])]) # doppler_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[17], calc_xis_lob[17])]) # localgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[18], calc_xis_lob[18])]) # doppler_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[19], calc_xis_lob[19])]) # integratedgp_doppler
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[20], calc_xis_lob[20])]) # lensing_localgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[21], calc_xis_lob[21])]) # localgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[22], calc_xis_lob[22])]) # lensing_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[23], calc_xis_lob[23])]) # integratedgp_lensing
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[24], calc_xis_lob[24])]) # localgp_integratedgp
        @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[25], calc_xis_lob[25])]) # integratedgp_localgp

        rm(name)
        for effect in GaPSE.GR_EFFECTS_GNC
            rm("all_GNC_standalones_CFs/xi_GNC_" * effect * "_L$L" * ".txt")
        end
        rm("all_GNC_standalones_CFs")

    end
end

println("\nEnded tests on print_map_sum_ξ_GNC_multipole!\n")
