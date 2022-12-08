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


println("\nNow I work on print_map_sum_ξ_LDxGNC_multipole.")
println("This will be longer...")


@testset "test print_map_sum_ξ_LDxGNC_multipole no_window, first method" begin
     RTOL = 2e-2

     kwargs = Dict(
          :use_windows => false,
          :enhancer => 1e8, :s1 => nothing,
          :N_trap => 30, :N_lob => 30,
          :atol_quad => 0.0, :rtol_quad => 1e-2,
          :N_χs => 40, :N_χs_2 => 20,
          :pr => false,
     )

     @testset "test_zeros" begin
          @test_throws AssertionError GaPSE.GaPSE.print_map_sum_ξ_LDxGNC_multipole(COSMO, "NONEXISTINGDIR/file.txt", [1.0, 2.0, 3.0])
     end

     @testset "L = 0 no_window" begin
          L = 0

          ##### quad ####
          ss_quad, res_sums_quad, res_xis_quad = GaPSE.readxyall(
               "datatest/LDxGNC_SumXiMultipoles/xis_LDxGNC_L$L" * "_quad_noF_specific_ss.txt", comments=true)

          name = "calc_xis_LDxGNC_L$L" * "_quad_noF_specific_ss.txt"
          isfile(name) && rm(name)

          GaPSE.print_map_sum_ξ_LDxGNC_multipole(COSMO, name, ss_quad; 
               single=true, L=L, alg=:quad, kwargs...)
          calc_ss_quad, calc_sums_quad, calc_xis_quad =  GaPSE.readxyall(name, comments=true)

          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(ss_quad, calc_ss_quad)])
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_sums_quad, calc_sums_quad)])
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[1], calc_xis_quad[1])]) # auto_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[2], calc_xis_quad[2])]) # auto_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[3], calc_xis_quad[3])]) # auto_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[4], calc_xis_quad[4])]) # auto_integrated
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[5], calc_xis_quad[5])]) # lensing_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[6], calc_xis_quad[6])]) # doppler_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[7], calc_xis_quad[7])]) # doppler_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[8], calc_xis_quad[8])]) # localgp_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[9], calc_xis_quad[9])]) # doppler_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[10], calc_xis_quad[10])]) # integratedgp_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[11], calc_xis_quad[11])]) # lensing_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[12], calc_xis_quad[12])]) # localgp_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[13], calc_xis_quad[13])]) # lensing_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[14], calc_xis_quad[14])]) # integratedgp_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[15], calc_xis_quad[15])]) # localgp_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[16], calc_xis_quad[16])]) # integratedgp_localgp

          rm(name)

          #### trap ####
          ss_trap, res_sums_trap, res_xis_trap = GaPSE.readxyall(
               "datatest/LDxGNC_SumXiMultipoles/xis_LDxGNC_L$L" * "_trap_noF_specific_ss.txt", comments=true)

          name = "calc_xis_LDxGNC_L$L" * "_trap_noF_specific_ss.txt"
          isfile(name) && rm(name)

          GaPSE.print_map_sum_ξ_LDxGNC_multipole(COSMO, name, ss_trap; 
               single=true, L=L, alg=:trap, kwargs...)
          calc_ss_trap, calc_sums_trap, calc_xis_trap =  GaPSE.readxyall(name, comments=true)

          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(ss_trap, calc_ss_trap)])
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_sums_trap, calc_sums_trap)])
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[1], calc_xis_trap[1])]) # auto_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[2], calc_xis_trap[2])]) # auto_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[3], calc_xis_trap[3])]) # auto_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[4], calc_xis_trap[4])]) # auto_integrated
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[5], calc_xis_trap[5])]) # lensing_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[6], calc_xis_trap[6])]) # doppler_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[7], calc_xis_trap[7])]) # doppler_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[8], calc_xis_trap[8])]) # localgp_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[9], calc_xis_trap[9])]) # doppler_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[10], calc_xis_trap[10])]) # integratedgp_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[11], calc_xis_trap[11])]) # lensing_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[12], calc_xis_trap[12])]) # localgp_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[13], calc_xis_trap[13])]) # lensing_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[14], calc_xis_trap[14])]) # integratedgp_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[15], calc_xis_trap[15])]) # localgp_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[16], calc_xis_trap[16])]) # integratedgp_localgp

          rm(name)

          #### lob ####
          ss_lob, res_sums_lob, res_xis_lob = GaPSE.readxyall(
               "datatest/LDxGNC_SumXiMultipoles/xis_LDxGNC_L$L" * "_lob_noF_specific_ss.txt", comments=true)

          name = "calc_xis_LDxGNC_L$L" * "_lob_noF_specific_ss.txt"
          isfile(name) && rm(name)

          GaPSE.print_map_sum_ξ_LDxGNC_multipole(COSMO, name, ss_lob; 
               single=true, L=L, alg=:lobatto, kwargs...)
          calc_ss_lob, calc_sums_lob, calc_xis_lob =  GaPSE.readxyall(name,  comments=true)

          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(ss_lob, calc_ss_lob)])
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_sums_lob, calc_sums_lob)])
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[1], calc_xis_lob[1])]) # auto_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[2], calc_xis_lob[2])]) # auto_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[3], calc_xis_lob[3])]) # auto_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[4], calc_xis_lob[4])]) # auto_integrated
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[5], calc_xis_lob[5])]) # lensing_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[6], calc_xis_lob[6])]) # doppler_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[7], calc_xis_lob[7])]) # doppler_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[8], calc_xis_lob[8])]) # localgp_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[9], calc_xis_lob[9])]) # doppler_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[10], calc_xis_lob[10])]) # integratedgp_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[11], calc_xis_lob[11])]) # lensing_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[12], calc_xis_lob[12])]) # localgp_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[13], calc_xis_lob[13])]) # lensing_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[14], calc_xis_lob[14])]) # integratedgp_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[15], calc_xis_lob[15])]) # localgp_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[16], calc_xis_lob[16])]) # integratedgp_localgp

          rm(name)

     end

     @testset "L = 1 no_window" begin
          L = 1

          #### trap ####
          ss_trap, res_sums_trap, res_xis_trap = GaPSE.readxyall(
               "datatest/LDxGNC_SumXiMultipoles/xis_LDxGNC_L$L" * "_trap_noF_specific_ss.txt", comments=true)

          name = "calc_xis_LDxGNC_L$L" * "_trap_noF_specific_ss.txt"
          isfile(name) && rm(name)

          GaPSE.print_map_sum_ξ_LDxGNC_multipole(COSMO, name, ss_trap; 
               single=true, L=L, alg=:trap, kwargs...)
          calc_ss_trap, calc_sums_trap, calc_xis_trap =  GaPSE.readxyall(name, comments=true)

          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(ss_trap, calc_ss_trap)])
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_sums_trap, calc_sums_trap)])
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[1], calc_xis_trap[1])]) # auto_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[2], calc_xis_trap[2])]) # auto_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[3], calc_xis_trap[3])]) # auto_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[4], calc_xis_trap[4])]) # auto_integrated
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[5], calc_xis_trap[5])]) # lensing_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[6], calc_xis_trap[6])]) # doppler_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[7], calc_xis_trap[7])]) # doppler_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[8], calc_xis_trap[8])]) # localgp_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[9], calc_xis_trap[9])]) # doppler_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[10], calc_xis_trap[10])]) # integratedgp_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[11], calc_xis_trap[11])]) # lensing_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[12], calc_xis_trap[12])]) # localgp_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[13], calc_xis_trap[13])]) # lensing_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[14], calc_xis_trap[14])]) # integratedgp_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[15], calc_xis_trap[15])]) # localgp_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[16], calc_xis_trap[16])]) # integratedgp_localgp

          rm(name)

          #### lob ####
          ss_lob, res_sums_lob, res_xis_lob = GaPSE.readxyall(
               "datatest/LDxGNC_SumXiMultipoles/xis_LDxGNC_L$L" * "_lob_noF_specific_ss.txt", comments=true)

          name = "calc_xis_LDxGNC_L$L" * "_lob_noF_specific_ss.txt"
          isfile(name) && rm(name)

          GaPSE.print_map_sum_ξ_LDxGNC_multipole(COSMO, name, ss_lob; 
               single=true, L=L, alg=:lobatto, kwargs...)
          calc_ss_lob, calc_sums_lob, calc_xis_lob =  GaPSE.readxyall(name,  comments=true)

          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(ss_lob, calc_ss_lob)])
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_sums_lob, calc_sums_lob)])
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[1], calc_xis_lob[1])]) # auto_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[2], calc_xis_lob[2])]) # auto_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[3], calc_xis_lob[3])]) # auto_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[4], calc_xis_lob[4])]) # auto_integrated
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[5], calc_xis_lob[5])]) # lensing_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[6], calc_xis_lob[6])]) # doppler_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[7], calc_xis_lob[7])]) # doppler_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[8], calc_xis_lob[8])]) # localgp_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[9], calc_xis_lob[9])]) # doppler_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[10], calc_xis_lob[10])]) # integratedgp_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[11], calc_xis_lob[11])]) # lensing_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[12], calc_xis_lob[12])]) # localgp_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[13], calc_xis_lob[13])]) # lensing_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[14], calc_xis_lob[14])]) # integratedgp_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[15], calc_xis_lob[15])]) # localgp_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[16], calc_xis_lob[16])]) # integratedgp_localgp

          rm(name)

     end

     @testset "L = 2 no_window" begin
          L = 2

          ##### quad ####
          ss_quad, res_sums_quad, res_xis_quad = GaPSE.readxyall(
               "datatest/LDxGNC_SumXiMultipoles/xis_LDxGNC_L$L" * "_quad_noF_specific_ss.txt", comments=true)

          name = "calc_xis_LDxGNC_L$L" * "_quad_noF_specific_ss.txt"
          isfile(name) && rm(name)

          GaPSE.print_map_sum_ξ_LDxGNC_multipole(COSMO, name, ss_quad; 
               single=true, L=L, alg=:quad, kwargs...)
          calc_ss_quad, calc_sums_quad, calc_xis_quad =  GaPSE.readxyall(name, comments=true)

          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(ss_quad, calc_ss_quad)])
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_sums_quad, calc_sums_quad)])
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[1], calc_xis_quad[1])]) # auto_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[2], calc_xis_quad[2])]) # auto_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[3], calc_xis_quad[3])]) # auto_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[4], calc_xis_quad[4])]) # auto_integrated
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[5], calc_xis_quad[5])]) # lensing_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[6], calc_xis_quad[6])]) # doppler_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[7], calc_xis_quad[7])]) # doppler_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[8], calc_xis_quad[8])]) # localgp_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[9], calc_xis_quad[9])]) # doppler_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[10], calc_xis_quad[10])]) # integratedgp_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[11], calc_xis_quad[11])]) # lensing_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[12], calc_xis_quad[12])]) # localgp_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[13], calc_xis_quad[13])]) # lensing_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[14], calc_xis_quad[14])]) # integratedgp_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[15], calc_xis_quad[15])]) # localgp_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[16], calc_xis_quad[16])]) # integratedgp_localgp

          rm(name)

          #### trap ####
          ss_trap, res_sums_trap, res_xis_trap = GaPSE.readxyall(
               "datatest/LDxGNC_SumXiMultipoles/xis_LDxGNC_L$L" * "_trap_noF_specific_ss.txt", comments=true)

          name = "calc_xis_LDxGNC_L$L" * "_trap_noF_specific_ss.txt"
          isfile(name) && rm(name)

          GaPSE.print_map_sum_ξ_LDxGNC_multipole(COSMO, name, ss_trap; 
               single=true, L=L, alg=:trap, kwargs...)
          calc_ss_trap, calc_sums_trap, calc_xis_trap =  GaPSE.readxyall(name, comments=true)

          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(ss_trap, calc_ss_trap)])
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_sums_trap, calc_sums_trap)])
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[1], calc_xis_trap[1])]) # auto_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[2], calc_xis_trap[2])]) # auto_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[3], calc_xis_trap[3])]) # auto_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[4], calc_xis_trap[4])]) # auto_integrated
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[5], calc_xis_trap[5])]) # lensing_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[6], calc_xis_trap[6])]) # doppler_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[7], calc_xis_trap[7])]) # doppler_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[8], calc_xis_trap[8])]) # localgp_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[9], calc_xis_trap[9])]) # doppler_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[10], calc_xis_trap[10])]) # integratedgp_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[11], calc_xis_trap[11])]) # lensing_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[12], calc_xis_trap[12])]) # localgp_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[13], calc_xis_trap[13])]) # lensing_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[14], calc_xis_trap[14])]) # integratedgp_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[15], calc_xis_trap[15])]) # localgp_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[16], calc_xis_trap[16])]) # integratedgp_localgp

          rm(name)

          #### lob ####
          ss_lob, res_sums_lob, res_xis_lob = GaPSE.readxyall(
               "datatest/LDxGNC_SumXiMultipoles/xis_LDxGNC_L$L" * "_lob_noF_specific_ss.txt", comments=true)

          name = "calc_xis_LDxGNC_L$L" * "_lob_noF_specific_ss.txt"
          isfile(name) && rm(name)

          GaPSE.print_map_sum_ξ_LDxGNC_multipole(COSMO, name, ss_lob; 
               single=true, L=L, alg=:lobatto, kwargs...)
          calc_ss_lob, calc_sums_lob, calc_xis_lob =  GaPSE.readxyall(name,  comments=true)

          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(ss_lob, calc_ss_lob)])
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_sums_lob, calc_sums_lob)])
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[1], calc_xis_lob[1])]) # auto_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[2], calc_xis_lob[2])]) # auto_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[3], calc_xis_lob[3])]) # auto_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[4], calc_xis_lob[4])]) # auto_integrated
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[5], calc_xis_lob[5])]) # lensing_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[6], calc_xis_lob[6])]) # doppler_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[7], calc_xis_lob[7])]) # doppler_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[8], calc_xis_lob[8])]) # localgp_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[9], calc_xis_lob[9])]) # doppler_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[10], calc_xis_lob[10])]) # integratedgp_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[11], calc_xis_lob[11])]) # lensing_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[12], calc_xis_lob[12])]) # localgp_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[13], calc_xis_lob[13])]) # lensing_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[14], calc_xis_lob[14])]) # integratedgp_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[15], calc_xis_lob[15])]) # localgp_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[16], calc_xis_lob[16])]) # integratedgp_localgp

          rm(name)

     end

     @testset "L = 3 no_window" begin
          L = 3


          #### trap ####
          ss_trap, res_sums_trap, res_xis_trap = GaPSE.readxyall(
               "datatest/LDxGNC_SumXiMultipoles/xis_LDxGNC_L$L" * "_trap_noF_specific_ss.txt", comments=true)

          name = "calc_xis_LDxGNC_L$L" * "_trap_noF_specific_ss.txt"
          isfile(name) && rm(name)

          GaPSE.print_map_sum_ξ_LDxGNC_multipole(COSMO, name, ss_trap; 
               single=true, L=L, alg=:trap, kwargs...)
          calc_ss_trap, calc_sums_trap, calc_xis_trap =  GaPSE.readxyall(name, comments=true)

          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(ss_trap, calc_ss_trap)])
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_sums_trap, calc_sums_trap)])
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[1], calc_xis_trap[1])]) # auto_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[2], calc_xis_trap[2])]) # auto_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[3], calc_xis_trap[3])]) # auto_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[4], calc_xis_trap[4])]) # auto_integrated
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[5], calc_xis_trap[5])]) # lensing_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[6], calc_xis_trap[6])]) # doppler_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[7], calc_xis_trap[7])]) # doppler_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[8], calc_xis_trap[8])]) # localgp_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[9], calc_xis_trap[9])]) # doppler_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[10], calc_xis_trap[10])]) # integratedgp_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[11], calc_xis_trap[11])]) # lensing_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[12], calc_xis_trap[12])]) # localgp_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[13], calc_xis_trap[13])]) # lensing_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[14], calc_xis_trap[14])]) # integratedgp_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[15], calc_xis_trap[15])]) # localgp_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[16], calc_xis_trap[16])]) # integratedgp_localgp

          rm(name)

          #### lob ####
          ss_lob, res_sums_lob, res_xis_lob = GaPSE.readxyall(
               "datatest/LDxGNC_SumXiMultipoles/xis_LDxGNC_L$L" * "_lob_noF_specific_ss.txt", comments=true)

          name = "calc_xis_LDxGNC_L$L" * "_lob_noF_specific_ss.txt"
          isfile(name) && rm(name)

          GaPSE.print_map_sum_ξ_LDxGNC_multipole(COSMO, name, ss_lob; 
               single=true, L=L, alg=:lobatto, kwargs...)
          calc_ss_lob, calc_sums_lob, calc_xis_lob =  GaPSE.readxyall(name,  comments=true)

          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(ss_lob, calc_ss_lob)])
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_sums_lob, calc_sums_lob)])
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[1], calc_xis_lob[1])]) # auto_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[2], calc_xis_lob[2])]) # auto_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[3], calc_xis_lob[3])]) # auto_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[4], calc_xis_lob[4])]) # auto_integrated
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[5], calc_xis_lob[5])]) # lensing_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[6], calc_xis_lob[6])]) # doppler_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[7], calc_xis_lob[7])]) # doppler_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[8], calc_xis_lob[8])]) # localgp_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[9], calc_xis_lob[9])]) # doppler_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[10], calc_xis_lob[10])]) # integratedgp_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[11], calc_xis_lob[11])]) # lensing_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[12], calc_xis_lob[12])]) # localgp_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[13], calc_xis_lob[13])]) # lensing_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[14], calc_xis_lob[14])]) # integratedgp_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[15], calc_xis_lob[15])]) # localgp_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[16], calc_xis_lob[16])]) # integratedgp_localgp

          rm(name)

     end

     @testset "L = 4 no_window" begin
          L = 4

          #### trap ####
          ss_trap, res_sums_trap, res_xis_trap = GaPSE.readxyall(
               "datatest/LDxGNC_SumXiMultipoles/xis_LDxGNC_L$L" * "_trap_noF_specific_ss.txt", comments=true)

          name = "calc_xis_LDxGNC_L$L" * "_trap_noF_specific_ss.txt"
          isfile(name) && rm(name)

          GaPSE.print_map_sum_ξ_LDxGNC_multipole(COSMO, name, ss_trap; 
               single=true, L=L, alg=:trap, kwargs...)
          calc_ss_trap, calc_sums_trap, calc_xis_trap =  GaPSE.readxyall(name, comments=true)

          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(ss_trap, calc_ss_trap)])
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_sums_trap, calc_sums_trap)])
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[1], calc_xis_trap[1])]) # auto_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[2], calc_xis_trap[2])]) # auto_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[3], calc_xis_trap[3])]) # auto_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[4], calc_xis_trap[4])]) # auto_integrated
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[5], calc_xis_trap[5])]) # lensing_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[6], calc_xis_trap[6])]) # doppler_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[7], calc_xis_trap[7])]) # doppler_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[8], calc_xis_trap[8])]) # localgp_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[9], calc_xis_trap[9])]) # doppler_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[10], calc_xis_trap[10])]) # integratedgp_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[11], calc_xis_trap[11])]) # lensing_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[12], calc_xis_trap[12])]) # localgp_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[13], calc_xis_trap[13])]) # lensing_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[14], calc_xis_trap[14])]) # integratedgp_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[15], calc_xis_trap[15])]) # localgp_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[16], calc_xis_trap[16])]) # integratedgp_localgp

          rm(name)

          #### lob ####
          ss_lob, res_sums_lob, res_xis_lob = GaPSE.readxyall(
               "datatest/LDxGNC_SumXiMultipoles/xis_LDxGNC_L$L" * "_lob_noF_specific_ss.txt", comments=true)

          name = "calc_xis_LDxGNC_L$L" * "_lob_noF_specific_ss.txt"
          isfile(name) && rm(name)

          GaPSE.print_map_sum_ξ_LDxGNC_multipole(COSMO, name, ss_lob; 
               single=true, L=L, alg=:lobatto, kwargs...)
          calc_ss_lob, calc_sums_lob, calc_xis_lob =  GaPSE.readxyall(name,  comments=true)

          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(ss_lob, calc_ss_lob)])
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_sums_lob, calc_sums_lob)])
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[1], calc_xis_lob[1])]) # auto_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[2], calc_xis_lob[2])]) # auto_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[3], calc_xis_lob[3])]) # auto_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[4], calc_xis_lob[4])]) # auto_integrated
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[5], calc_xis_lob[5])]) # lensing_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[6], calc_xis_lob[6])]) # doppler_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[7], calc_xis_lob[7])]) # doppler_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[8], calc_xis_lob[8])]) # localgp_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[9], calc_xis_lob[9])]) # doppler_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[10], calc_xis_lob[10])]) # integratedgp_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[11], calc_xis_lob[11])]) # lensing_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[12], calc_xis_lob[12])]) # localgp_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[13], calc_xis_lob[13])]) # lensing_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[14], calc_xis_lob[14])]) # integratedgp_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[15], calc_xis_lob[15])]) # localgp_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[16], calc_xis_lob[16])]) # integratedgp_localgp

          rm(name)

     end
end

println("\nprocessing...")

@testset "test print_map_sum_ξ_LDxGNC_multipole with_window, first method" begin
     RTOL = 2e-2

     kwargs = Dict(
          :use_windows => true,
          :enhancer => 1e8, :s1 => nothing,
          :N_trap => 30, :N_lob => 30,
          :atol_quad => 0.0, :rtol_quad => 1e-2,
          :N_χs => 40, :N_χs_2 => 20,
          :pr => false,
     )

     @testset "test_zeros" begin
          @test_throws AssertionError GaPSE.GaPSE.print_map_sum_ξ_LDxGNC_multipole(COSMO, "NONEXISTINGDIR/file.txt", [1.0, 2.0, 3.0])
     end

     @testset "L = 0 with_window" begin
          L = 0

          ##### quad ####
          ss_quad, res_sums_quad, res_xis_quad = GaPSE.readxyall(
               "datatest/LDxGNC_SumXiMultipoles/xis_LDxGNC_L$L" * "_quad_withF_specific_ss.txt", comments=true)

          name = "calc_xis_LDxGNC_L$L" * "_quad_withF_specific_ss.txt"
          isfile(name) && rm(name)

          GaPSE.print_map_sum_ξ_LDxGNC_multipole(COSMO, name, ss_quad; 
               single=true, L=L, alg=:quad, kwargs...)
          calc_ss_quad, calc_sums_quad, calc_xis_quad =  GaPSE.readxyall(name, comments=true)

          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(ss_quad, calc_ss_quad)])
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_sums_quad, calc_sums_quad)])
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[1], calc_xis_quad[1])]) # auto_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[2], calc_xis_quad[2])]) # auto_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[3], calc_xis_quad[3])]) # auto_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[4], calc_xis_quad[4])]) # auto_integrated
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[5], calc_xis_quad[5])]) # lensing_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[6], calc_xis_quad[6])]) # doppler_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[7], calc_xis_quad[7])]) # doppler_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[8], calc_xis_quad[8])]) # localgp_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[9], calc_xis_quad[9])]) # doppler_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[10], calc_xis_quad[10])]) # integratedgp_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[11], calc_xis_quad[11])]) # lensing_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[12], calc_xis_quad[12])]) # localgp_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[13], calc_xis_quad[13])]) # lensing_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[14], calc_xis_quad[14])]) # integratedgp_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[15], calc_xis_quad[15])]) # localgp_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[16], calc_xis_quad[16])]) # integratedgp_localgp

          rm(name)

          #### trap ####
          ss_trap, res_sums_trap, res_xis_trap = GaPSE.readxyall(
               "datatest/LDxGNC_SumXiMultipoles/xis_LDxGNC_L$L" * "_trap_withF_specific_ss.txt", comments=true)

          name = "calc_xis_LDxGNC_L$L" * "_trap_withF_specific_ss.txt"
          isfile(name) && rm(name)

          GaPSE.print_map_sum_ξ_LDxGNC_multipole(COSMO, name, ss_trap; 
               single=true, L=L, alg=:trap, kwargs...)
          calc_ss_trap, calc_sums_trap, calc_xis_trap =  GaPSE.readxyall(name, comments=true)

          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(ss_trap, calc_ss_trap)])
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_sums_trap, calc_sums_trap)])
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[1], calc_xis_trap[1])]) # auto_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[2], calc_xis_trap[2])]) # auto_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[3], calc_xis_trap[3])]) # auto_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[4], calc_xis_trap[4])]) # auto_integrated
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[5], calc_xis_trap[5])]) # lensing_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[6], calc_xis_trap[6])]) # doppler_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[7], calc_xis_trap[7])]) # doppler_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[8], calc_xis_trap[8])]) # localgp_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[9], calc_xis_trap[9])]) # doppler_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[10], calc_xis_trap[10])]) # integratedgp_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[11], calc_xis_trap[11])]) # lensing_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[12], calc_xis_trap[12])]) # localgp_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[13], calc_xis_trap[13])]) # lensing_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[14], calc_xis_trap[14])]) # integratedgp_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[15], calc_xis_trap[15])]) # localgp_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[16], calc_xis_trap[16])]) # integratedgp_localgp

          rm(name)

          #### lob ####
          ss_lob, res_sums_lob, res_xis_lob = GaPSE.readxyall(
               "datatest/LDxGNC_SumXiMultipoles/xis_LDxGNC_L$L" * "_lob_withF_specific_ss.txt", comments=true)

          name = "calc_xis_LDxGNC_L$L" * "_lob_withF_specific_ss.txt"
          isfile(name) && rm(name)

          GaPSE.print_map_sum_ξ_LDxGNC_multipole(COSMO, name, ss_lob; 
               single=true, L=L, alg=:lobatto, kwargs...)
          calc_ss_lob, calc_sums_lob, calc_xis_lob =  GaPSE.readxyall(name,  comments=true)

          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(ss_lob, calc_ss_lob)])
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_sums_lob, calc_sums_lob)])
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[1], calc_xis_lob[1])]) # auto_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[2], calc_xis_lob[2])]) # auto_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[3], calc_xis_lob[3])]) # auto_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[4], calc_xis_lob[4])]) # auto_integrated
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[5], calc_xis_lob[5])]) # lensing_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[6], calc_xis_lob[6])]) # doppler_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[7], calc_xis_lob[7])]) # doppler_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[8], calc_xis_lob[8])]) # localgp_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[9], calc_xis_lob[9])]) # doppler_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[10], calc_xis_lob[10])]) # integratedgp_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[11], calc_xis_lob[11])]) # lensing_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[12], calc_xis_lob[12])]) # localgp_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[13], calc_xis_lob[13])]) # lensing_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[14], calc_xis_lob[14])]) # integratedgp_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[15], calc_xis_lob[15])]) # localgp_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[16], calc_xis_lob[16])]) # integratedgp_localgp

          rm(name)

     end

     @testset "L = 1 with_window" begin
          L = 1

          #### trap ####
          ss_trap, res_sums_trap, res_xis_trap = GaPSE.readxyall(
               "datatest/LDxGNC_SumXiMultipoles/xis_LDxGNC_L$L" * "_trap_withF_specific_ss.txt", comments=true)

          name = "calc_xis_LDxGNC_L$L" * "_trap_withF_specific_ss.txt"
          isfile(name) && rm(name)

          GaPSE.print_map_sum_ξ_LDxGNC_multipole(COSMO, name, ss_trap; 
               single=true, L=L, alg=:trap, kwargs...)
          calc_ss_trap, calc_sums_trap, calc_xis_trap =  GaPSE.readxyall(name, comments=true)

          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(ss_trap, calc_ss_trap)])
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_sums_trap, calc_sums_trap)])
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[1], calc_xis_trap[1])]) # auto_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[2], calc_xis_trap[2])]) # auto_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[3], calc_xis_trap[3])]) # auto_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[4], calc_xis_trap[4])]) # auto_integrated
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[5], calc_xis_trap[5])]) # lensing_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[6], calc_xis_trap[6])]) # doppler_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[7], calc_xis_trap[7])]) # doppler_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[8], calc_xis_trap[8])]) # localgp_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[9], calc_xis_trap[9])]) # doppler_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[10], calc_xis_trap[10])]) # integratedgp_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[11], calc_xis_trap[11])]) # lensing_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[12], calc_xis_trap[12])]) # localgp_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[13], calc_xis_trap[13])]) # lensing_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[14], calc_xis_trap[14])]) # integratedgp_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[15], calc_xis_trap[15])]) # localgp_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[16], calc_xis_trap[16])]) # integratedgp_localgp

          rm(name)

          #### lob ####
          ss_lob, res_sums_lob, res_xis_lob = GaPSE.readxyall(
               "datatest/LDxGNC_SumXiMultipoles/xis_LDxGNC_L$L" * "_lob_withF_specific_ss.txt", comments=true)

          name = "calc_xis_LDxGNC_L$L" * "_lob_withF_specific_ss.txt"
          isfile(name) && rm(name)

          GaPSE.print_map_sum_ξ_LDxGNC_multipole(COSMO, name, ss_lob; 
               single=true, L=L, alg=:lobatto, kwargs...)
          calc_ss_lob, calc_sums_lob, calc_xis_lob =  GaPSE.readxyall(name,  comments=true)

          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(ss_lob, calc_ss_lob)])
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_sums_lob, calc_sums_lob)])
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[1], calc_xis_lob[1])]) # auto_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[2], calc_xis_lob[2])]) # auto_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[3], calc_xis_lob[3])]) # auto_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[4], calc_xis_lob[4])]) # auto_integrated
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[5], calc_xis_lob[5])]) # lensing_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[6], calc_xis_lob[6])]) # doppler_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[7], calc_xis_lob[7])]) # doppler_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[8], calc_xis_lob[8])]) # localgp_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[9], calc_xis_lob[9])]) # doppler_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[10], calc_xis_lob[10])]) # integratedgp_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[11], calc_xis_lob[11])]) # lensing_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[12], calc_xis_lob[12])]) # localgp_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[13], calc_xis_lob[13])]) # lensing_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[14], calc_xis_lob[14])]) # integratedgp_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[15], calc_xis_lob[15])]) # localgp_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[16], calc_xis_lob[16])]) # integratedgp_localgp

          rm(name)

     end

     @testset "L = 2 with_window" begin
          L = 2

          ##### quad ####
          ss_quad, res_sums_quad, res_xis_quad = GaPSE.readxyall(
               "datatest/LDxGNC_SumXiMultipoles/xis_LDxGNC_L$L" * "_quad_withF_specific_ss.txt", comments=true)

          name = "calc_xis_LDxGNC_L$L" * "_quad_withF_specific_ss.txt"
          isfile(name) && rm(name)

          GaPSE.print_map_sum_ξ_LDxGNC_multipole(COSMO, name, ss_quad; 
               single=true, L=L, alg=:quad, kwargs...)
          calc_ss_quad, calc_sums_quad, calc_xis_quad =  GaPSE.readxyall(name, comments=true)

          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(ss_quad, calc_ss_quad)])
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_sums_quad, calc_sums_quad)])
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[1], calc_xis_quad[1])]) # auto_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[2], calc_xis_quad[2])]) # auto_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[3], calc_xis_quad[3])]) # auto_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[4], calc_xis_quad[4])]) # auto_integrated
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[5], calc_xis_quad[5])]) # lensing_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[6], calc_xis_quad[6])]) # doppler_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[7], calc_xis_quad[7])]) # doppler_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[8], calc_xis_quad[8])]) # localgp_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[9], calc_xis_quad[9])]) # doppler_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[10], calc_xis_quad[10])]) # integratedgp_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[11], calc_xis_quad[11])]) # lensing_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[12], calc_xis_quad[12])]) # localgp_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[13], calc_xis_quad[13])]) # lensing_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[14], calc_xis_quad[14])]) # integratedgp_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[15], calc_xis_quad[15])]) # localgp_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[16], calc_xis_quad[16])]) # integratedgp_localgp

          rm(name)

          #### trap ####
          ss_trap, res_sums_trap, res_xis_trap = GaPSE.readxyall(
               "datatest/LDxGNC_SumXiMultipoles/xis_LDxGNC_L$L" * "_trap_withF_specific_ss.txt", comments=true)

          name = "calc_xis_LDxGNC_L$L" * "_trap_withF_specific_ss.txt"
          isfile(name) && rm(name)

          GaPSE.print_map_sum_ξ_LDxGNC_multipole(COSMO, name, ss_trap; 
               single=true, L=L, alg=:trap, kwargs...)
          calc_ss_trap, calc_sums_trap, calc_xis_trap =  GaPSE.readxyall(name, comments=true)

          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(ss_trap, calc_ss_trap)])
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_sums_trap, calc_sums_trap)])
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[1], calc_xis_trap[1])]) # auto_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[2], calc_xis_trap[2])]) # auto_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[3], calc_xis_trap[3])]) # auto_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[4], calc_xis_trap[4])]) # auto_integrated
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[5], calc_xis_trap[5])]) # lensing_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[6], calc_xis_trap[6])]) # doppler_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[7], calc_xis_trap[7])]) # doppler_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[8], calc_xis_trap[8])]) # localgp_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[9], calc_xis_trap[9])]) # doppler_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[10], calc_xis_trap[10])]) # integratedgp_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[11], calc_xis_trap[11])]) # lensing_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[12], calc_xis_trap[12])]) # localgp_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[13], calc_xis_trap[13])]) # lensing_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[14], calc_xis_trap[14])]) # integratedgp_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[15], calc_xis_trap[15])]) # localgp_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[16], calc_xis_trap[16])]) # integratedgp_localgp

          rm(name)

          #### lob ####
          ss_lob, res_sums_lob, res_xis_lob = GaPSE.readxyall(
               "datatest/LDxGNC_SumXiMultipoles/xis_LDxGNC_L$L" * "_lob_withF_specific_ss.txt", comments=true)

          name = "calc_xis_LDxGNC_L$L" * "_lob_withF_specific_ss.txt"
          isfile(name) && rm(name)

          GaPSE.print_map_sum_ξ_LDxGNC_multipole(COSMO, name, ss_lob; 
               single=true, L=L, alg=:lobatto, kwargs...)
          calc_ss_lob, calc_sums_lob, calc_xis_lob =  GaPSE.readxyall(name,  comments=true)

          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(ss_lob, calc_ss_lob)])
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_sums_lob, calc_sums_lob)])
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[1], calc_xis_lob[1])]) # auto_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[2], calc_xis_lob[2])]) # auto_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[3], calc_xis_lob[3])]) # auto_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[4], calc_xis_lob[4])]) # auto_integrated
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[5], calc_xis_lob[5])]) # lensing_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[6], calc_xis_lob[6])]) # doppler_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[7], calc_xis_lob[7])]) # doppler_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[8], calc_xis_lob[8])]) # localgp_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[9], calc_xis_lob[9])]) # doppler_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[10], calc_xis_lob[10])]) # integratedgp_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[11], calc_xis_lob[11])]) # lensing_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[12], calc_xis_lob[12])]) # localgp_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[13], calc_xis_lob[13])]) # lensing_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[14], calc_xis_lob[14])]) # integratedgp_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[15], calc_xis_lob[15])]) # localgp_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[16], calc_xis_lob[16])]) # integratedgp_localgp

          rm(name)

     end

     @testset "L = 3 with_window" begin
          L = 3

          #### trap ####
          ss_trap, res_sums_trap, res_xis_trap = GaPSE.readxyall(
               "datatest/LDxGNC_SumXiMultipoles/xis_LDxGNC_L$L" * "_trap_withF_specific_ss.txt", comments=true)

          name = "calc_xis_LDxGNC_L$L" * "_trap_withF_specific_ss.txt"
          isfile(name) && rm(name)

          GaPSE.print_map_sum_ξ_LDxGNC_multipole(COSMO, name, ss_trap; 
               single=true, L=L, alg=:trap, kwargs...)
          calc_ss_trap, calc_sums_trap, calc_xis_trap =  GaPSE.readxyall(name, comments=true)

          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(ss_trap, calc_ss_trap)])
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_sums_trap, calc_sums_trap)])
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[1], calc_xis_trap[1])]) # auto_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[2], calc_xis_trap[2])]) # auto_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[3], calc_xis_trap[3])]) # auto_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[4], calc_xis_trap[4])]) # auto_integrated
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[5], calc_xis_trap[5])]) # lensing_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[6], calc_xis_trap[6])]) # doppler_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[7], calc_xis_trap[7])]) # doppler_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[8], calc_xis_trap[8])]) # localgp_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[9], calc_xis_trap[9])]) # doppler_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[10], calc_xis_trap[10])]) # integratedgp_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[11], calc_xis_trap[11])]) # lensing_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[12], calc_xis_trap[12])]) # localgp_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[13], calc_xis_trap[13])]) # lensing_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[14], calc_xis_trap[14])]) # integratedgp_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[15], calc_xis_trap[15])]) # localgp_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[16], calc_xis_trap[16])]) # integratedgp_localgp

          rm(name)

          #### lob ####
          ss_lob, res_sums_lob, res_xis_lob = GaPSE.readxyall(
               "datatest/LDxGNC_SumXiMultipoles/xis_LDxGNC_L$L" * "_lob_withF_specific_ss.txt", comments=true)

          name = "calc_xis_LDxGNC_L$L" * "_lob_withF_specific_ss.txt"
          isfile(name) && rm(name)

          GaPSE.print_map_sum_ξ_LDxGNC_multipole(COSMO, name, ss_lob; 
               single=true, L=L, alg=:lobatto, kwargs...)
          calc_ss_lob, calc_sums_lob, calc_xis_lob =  GaPSE.readxyall(name,  comments=true)

          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(ss_lob, calc_ss_lob)])
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_sums_lob, calc_sums_lob)])
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[1], calc_xis_lob[1])]) # auto_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[2], calc_xis_lob[2])]) # auto_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[3], calc_xis_lob[3])]) # auto_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[4], calc_xis_lob[4])]) # auto_integrated
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[5], calc_xis_lob[5])]) # lensing_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[6], calc_xis_lob[6])]) # doppler_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[7], calc_xis_lob[7])]) # doppler_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[8], calc_xis_lob[8])]) # localgp_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[9], calc_xis_lob[9])]) # doppler_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[10], calc_xis_lob[10])]) # integratedgp_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[11], calc_xis_lob[11])]) # lensing_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[12], calc_xis_lob[12])]) # localgp_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[13], calc_xis_lob[13])]) # lensing_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[14], calc_xis_lob[14])]) # integratedgp_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[15], calc_xis_lob[15])]) # localgp_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[16], calc_xis_lob[16])]) # integratedgp_localgp

          rm(name)

     end

     @testset "L = 4 with_window" begin
          L = 4

          #### trap ####
          ss_trap, res_sums_trap, res_xis_trap = GaPSE.readxyall(
               "datatest/LDxGNC_SumXiMultipoles/xis_LDxGNC_L$L" * "_trap_withF_specific_ss.txt", comments=true)

          name = "calc_xis_LDxGNC_L$L" * "_trap_withF_specific_ss.txt"
          isfile(name) && rm(name)

          GaPSE.print_map_sum_ξ_LDxGNC_multipole(COSMO, name, ss_trap; 
               single=true, L=L, alg=:trap, kwargs...)
          calc_ss_trap, calc_sums_trap, calc_xis_trap =  GaPSE.readxyall(name, comments=true)

          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(ss_trap, calc_ss_trap)])
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_sums_trap, calc_sums_trap)])
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[1], calc_xis_trap[1])]) # auto_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[2], calc_xis_trap[2])]) # auto_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[3], calc_xis_trap[3])]) # auto_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[4], calc_xis_trap[4])]) # auto_integrated
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[5], calc_xis_trap[5])]) # lensing_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[6], calc_xis_trap[6])]) # doppler_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[7], calc_xis_trap[7])]) # doppler_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[8], calc_xis_trap[8])]) # localgp_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[9], calc_xis_trap[9])]) # doppler_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[10], calc_xis_trap[10])]) # integratedgp_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[11], calc_xis_trap[11])]) # lensing_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[12], calc_xis_trap[12])]) # localgp_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[13], calc_xis_trap[13])]) # lensing_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[14], calc_xis_trap[14])]) # integratedgp_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[15], calc_xis_trap[15])]) # localgp_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[16], calc_xis_trap[16])]) # integratedgp_localgp

          rm(name)

          #### lob ####
          ss_lob, res_sums_lob, res_xis_lob = GaPSE.readxyall(
               "datatest/LDxGNC_SumXiMultipoles/xis_LDxGNC_L$L" * "_lob_withF_specific_ss.txt", comments=true)

          name = "calc_xis_LDxGNC_L$L" * "_lob_withF_specific_ss.txt"
          isfile(name) && rm(name)

          GaPSE.print_map_sum_ξ_LDxGNC_multipole(COSMO, name, ss_lob; 
               single=true, L=L, alg=:lobatto, kwargs...)
          calc_ss_lob, calc_sums_lob, calc_xis_lob =  GaPSE.readxyall(name,  comments=true)

          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(ss_lob, calc_ss_lob)])
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_sums_lob, calc_sums_lob)])
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[1], calc_xis_lob[1])]) # auto_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[2], calc_xis_lob[2])]) # auto_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[3], calc_xis_lob[3])]) # auto_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[4], calc_xis_lob[4])]) # auto_integrated
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[5], calc_xis_lob[5])]) # lensing_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[6], calc_xis_lob[6])]) # doppler_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[7], calc_xis_lob[7])]) # doppler_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[8], calc_xis_lob[8])]) # localgp_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[9], calc_xis_lob[9])]) # doppler_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[10], calc_xis_lob[10])]) # integratedgp_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[11], calc_xis_lob[11])]) # lensing_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[12], calc_xis_lob[12])]) # localgp_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[13], calc_xis_lob[13])]) # lensing_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[14], calc_xis_lob[14])]) # integratedgp_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[15], calc_xis_lob[15])]) # localgp_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[16], calc_xis_lob[16])]) # integratedgp_localgp

          rm(name)

     end
end

println("\nI am at half of this path...")

##########################################################################################92


@testset "test print_map_sum_ξ_LDxGNC_multipole no_window second method" begin
     RTOL = 2e-2

     kwargs = Dict(
          :use_windows => false,
          :enhancer => 1e8, :s1 => nothing,
          :N_trap => 30, :N_lob => 30,
          :atol_quad => 0.0, :rtol_quad => 1e-2,
          :N_χs => 40, :N_χs_2 => 20,
          :pr => false,
     )

     @testset "test_zeros" begin
          @test_throws AssertionError GaPSE.GaPSE.print_map_sum_ξ_LDxGNC_multipole(COSMO, "NONEXISTINGDIR/file.txt", [1.0, 2.0, 3.0])
     end

     @testset "L = 0 no_window" begin
          L = 0

          ##### quad ####
          ss_quad, res_sums_quad, res_xis_quad = GaPSE.readxyall(
               "datatest/LDxGNC_SumXiMultipoles/xis_LDxGNC_L$L" * "_quad_noF_specific_ss.txt", comments=true)

          name = "calc_xis_LDxGNC_L$L" * "_quad_noF_specific_ss.txt"
          isfile(name) && rm(name)

          GaPSE.print_map_sum_ξ_LDxGNC_multipole(COSMO, name, ss_quad; 
               single=false, L=L, alg=:quad, kwargs...)

          calc_ss_quad, calc_sums_quad =  GaPSE.readxy(name, comments=true)
          calc_xis_quad = [
               begin
                    a_name = "all_LDxGNC_standalones_CFs/xi_LDxGNC_" * effect * "_L$L" * ".txt"
                    _, xis = GaPSE.readxy(a_name, comments = true)
                    xis
               end for effect in GaPSE.GR_EFFECTS_LDxGNC
          ]

          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(ss_quad, calc_ss_quad)])
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_sums_quad, calc_sums_quad)])
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[1], calc_xis_quad[1])]) # auto_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[2], calc_xis_quad[2])]) # auto_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[3], calc_xis_quad[3])]) # auto_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[4], calc_xis_quad[4])]) # auto_integrated
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[5], calc_xis_quad[5])]) # lensing_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[6], calc_xis_quad[6])]) # doppler_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[7], calc_xis_quad[7])]) # doppler_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[8], calc_xis_quad[8])]) # localgp_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[9], calc_xis_quad[9])]) # doppler_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[10], calc_xis_quad[10])]) # integratedgp_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[11], calc_xis_quad[11])]) # lensing_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[12], calc_xis_quad[12])]) # localgp_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[13], calc_xis_quad[13])]) # lensing_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[14], calc_xis_quad[14])]) # integratedgp_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[15], calc_xis_quad[15])]) # localgp_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[16], calc_xis_quad[16])]) # integratedgp_localgp

        
          rm(name)
          for effect in GaPSE.GR_EFFECTS_LDxGNC
               rm("all_LDxGNC_standalones_CFs/xi_LDxGNC_" * effect * "_L$L" * ".txt")
          end
          rm("all_LDxGNC_standalones_CFs")

          #### trap ####
          ss_trap, res_sums_trap, res_xis_trap = GaPSE.readxyall(
               "datatest/LDxGNC_SumXiMultipoles/xis_LDxGNC_L$L" * "_trap_noF_specific_ss.txt", comments=true)

          name = "calc_xis_LDxGNC_L$L" * "_trap_noF_specific_ss.txt"
          isfile(name) && rm(name)

          GaPSE.print_map_sum_ξ_LDxGNC_multipole(COSMO, name, ss_trap; 
               single=false, L=L, alg=:trap, kwargs...)
          calc_ss_trap, calc_sums_trap =  GaPSE.readxy(name, comments=true)
          calc_xis_trap = [
               begin
                    a_name = "all_LDxGNC_standalones_CFs/xi_LDxGNC_" * effect * "_L$L" * ".txt"
                    _, xis = GaPSE.readxy(a_name, comments = true)
                    xis
               end for effect in GaPSE.GR_EFFECTS_LDxGNC
          ]

          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(ss_trap, calc_ss_trap)])
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_sums_trap, calc_sums_trap)])
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[1], calc_xis_trap[1])]) # auto_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[2], calc_xis_trap[2])]) # auto_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[3], calc_xis_trap[3])]) # auto_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[4], calc_xis_trap[4])]) # auto_integrated
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[5], calc_xis_trap[5])]) # lensing_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[6], calc_xis_trap[6])]) # doppler_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[7], calc_xis_trap[7])]) # doppler_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[8], calc_xis_trap[8])]) # localgp_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[9], calc_xis_trap[9])]) # doppler_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[10], calc_xis_trap[10])]) # integratedgp_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[11], calc_xis_trap[11])]) # lensing_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[12], calc_xis_trap[12])]) # localgp_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[13], calc_xis_trap[13])]) # lensing_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[14], calc_xis_trap[14])]) # integratedgp_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[15], calc_xis_trap[15])]) # localgp_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[16], calc_xis_trap[16])]) # integratedgp_localgp

          rm(name)
          for effect in GaPSE.GR_EFFECTS_LDxGNC
               rm("all_LDxGNC_standalones_CFs/xi_LDxGNC_" * effect * "_L$L" * ".txt")
          end
          rm("all_LDxGNC_standalones_CFs")

          #### lob ####
          ss_lob, res_sums_lob, res_xis_lob = GaPSE.readxyall(
               "datatest/LDxGNC_SumXiMultipoles/xis_LDxGNC_L$L" * "_lob_noF_specific_ss.txt", comments=true)

          name = "calc_xis_LDxGNC_L$L" * "_lob_noF_specific_ss.txt"
          isfile(name) && rm(name)

          GaPSE.print_map_sum_ξ_LDxGNC_multipole(COSMO, name, ss_lob; 
               single=false, L=L, alg=:lobatto, kwargs...)
          calc_ss_lob, calc_sums_lob, calc_xis_lob =  GaPSE.readxyall(name,  comments=true)
          calc_ss_lob, calc_sums_lob =  GaPSE.readxy(name, comments=true)
          calc_xis_lob = [
               begin
                    a_name = "all_LDxGNC_standalones_CFs/xi_LDxGNC_" * effect * "_L$L" * ".txt"
                    _, xis = GaPSE.readxy(a_name, comments = true)
                    xis
               end for effect in GaPSE.GR_EFFECTS_LDxGNC
          ]

          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(ss_lob, calc_ss_lob)])
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_sums_lob, calc_sums_lob)])
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[1], calc_xis_lob[1])]) # auto_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[2], calc_xis_lob[2])]) # auto_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[3], calc_xis_lob[3])]) # auto_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[4], calc_xis_lob[4])]) # auto_integrated
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[5], calc_xis_lob[5])]) # lensing_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[6], calc_xis_lob[6])]) # doppler_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[7], calc_xis_lob[7])]) # doppler_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[8], calc_xis_lob[8])]) # localgp_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[9], calc_xis_lob[9])]) # doppler_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[10], calc_xis_lob[10])]) # integratedgp_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[11], calc_xis_lob[11])]) # lensing_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[12], calc_xis_lob[12])]) # localgp_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[13], calc_xis_lob[13])]) # lensing_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[14], calc_xis_lob[14])]) # integratedgp_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[15], calc_xis_lob[15])]) # localgp_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[16], calc_xis_lob[16])]) # integratedgp_localgp

          rm(name)
          for effect in GaPSE.GR_EFFECTS_LDxGNC
               rm("all_LDxGNC_standalones_CFs/xi_LDxGNC_" * effect * "_L$L" * ".txt")
          end
          rm("all_LDxGNC_standalones_CFs")

     end

     @testset "L = 1 no_window" begin
          L = 1

          #### trap ####
          ss_trap, res_sums_trap, res_xis_trap = GaPSE.readxyall(
               "datatest/LDxGNC_SumXiMultipoles/xis_LDxGNC_L$L" * "_trap_noF_specific_ss.txt", comments=true)

          name = "calc_xis_LDxGNC_L$L" * "_trap_noF_specific_ss.txt"
          isfile(name) && rm(name)

          GaPSE.print_map_sum_ξ_LDxGNC_multipole(COSMO, name, ss_trap; 
               single=false, L=L, alg=:trap, kwargs...)
          calc_ss_trap, calc_sums_trap =  GaPSE.readxy(name, comments=true)
          calc_xis_trap = [
               begin
                    a_name = "all_LDxGNC_standalones_CFs/xi_LDxGNC_" * effect * "_L$L" * ".txt"
                    _, xis = GaPSE.readxy(a_name, comments = true)
                    xis
               end for effect in GaPSE.GR_EFFECTS_LDxGNC
          ]

          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(ss_trap, calc_ss_trap)])
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_sums_trap, calc_sums_trap)])
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[1], calc_xis_trap[1])]) # auto_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[2], calc_xis_trap[2])]) # auto_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[3], calc_xis_trap[3])]) # auto_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[4], calc_xis_trap[4])]) # auto_integrated
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[5], calc_xis_trap[5])]) # lensing_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[6], calc_xis_trap[6])]) # doppler_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[7], calc_xis_trap[7])]) # doppler_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[8], calc_xis_trap[8])]) # localgp_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[9], calc_xis_trap[9])]) # doppler_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[10], calc_xis_trap[10])]) # integratedgp_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[11], calc_xis_trap[11])]) # lensing_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[12], calc_xis_trap[12])]) # localgp_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[13], calc_xis_trap[13])]) # lensing_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[14], calc_xis_trap[14])]) # integratedgp_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[15], calc_xis_trap[15])]) # localgp_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[16], calc_xis_trap[16])]) # integratedgp_localgp

          rm(name)
          for effect in GaPSE.GR_EFFECTS_LDxGNC
               rm("all_LDxGNC_standalones_CFs/xi_LDxGNC_" * effect * "_L$L" * ".txt")
          end
          rm("all_LDxGNC_standalones_CFs")

          #### lob ####
          ss_lob, res_sums_lob, res_xis_lob = GaPSE.readxyall(
               "datatest/LDxGNC_SumXiMultipoles/xis_LDxGNC_L$L" * "_lob_noF_specific_ss.txt", comments=true)

          name = "calc_xis_LDxGNC_L$L" * "_lob_noF_specific_ss.txt"
          isfile(name) && rm(name)

          GaPSE.print_map_sum_ξ_LDxGNC_multipole(COSMO, name, ss_lob; 
               single=false, L=L, alg=:lobatto, kwargs...)
          calc_ss_lob, calc_sums_lob =  GaPSE.readxy(name, comments=true)
          calc_xis_lob = [
               begin
                    a_name = "all_LDxGNC_standalones_CFs/xi_LDxGNC_" * effect * "_L$L" * ".txt"
                    _, xis = GaPSE.readxy(a_name, comments = true)
                    xis
               end for effect in GaPSE.GR_EFFECTS_LDxGNC
          ]

          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(ss_lob, calc_ss_lob)])
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_sums_lob, calc_sums_lob)])
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[1], calc_xis_lob[1])]) # auto_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[2], calc_xis_lob[2])]) # auto_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[3], calc_xis_lob[3])]) # auto_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[4], calc_xis_lob[4])]) # auto_integrated
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[5], calc_xis_lob[5])]) # lensing_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[6], calc_xis_lob[6])]) # doppler_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[7], calc_xis_lob[7])]) # doppler_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[8], calc_xis_lob[8])]) # localgp_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[9], calc_xis_lob[9])]) # doppler_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[10], calc_xis_lob[10])]) # integratedgp_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[11], calc_xis_lob[11])]) # lensing_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[12], calc_xis_lob[12])]) # localgp_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[13], calc_xis_lob[13])]) # lensing_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[14], calc_xis_lob[14])]) # integratedgp_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[15], calc_xis_lob[15])]) # localgp_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[16], calc_xis_lob[16])]) # integratedgp_localgp

          rm(name)
          for effect in GaPSE.GR_EFFECTS_LDxGNC
               rm("all_LDxGNC_standalones_CFs/xi_LDxGNC_" * effect * "_L$L" * ".txt")
          end
          rm("all_LDxGNC_standalones_CFs")

     end

     @testset "L = 2 no_window" begin
          L = 2

          ##### quad ####
          ss_quad, res_sums_quad, res_xis_quad = GaPSE.readxyall(
               "datatest/LDxGNC_SumXiMultipoles/xis_LDxGNC_L$L" * "_quad_noF_specific_ss.txt", comments=true)

          name = "calc_xis_LDxGNC_L$L" * "_quad_noF_specific_ss.txt"
          isfile(name) && rm(name)

          GaPSE.print_map_sum_ξ_LDxGNC_multipole(COSMO, name, ss_quad; 
               single=false, L=L, alg=:quad, kwargs...)
          calc_ss_quad, calc_sums_quad =  GaPSE.readxy(name, comments=true)
          calc_xis_quad = [
               begin
                    a_name = "all_LDxGNC_standalones_CFs/xi_LDxGNC_" * effect * "_L$L" * ".txt"
                    _, xis = GaPSE.readxy(a_name, comments = true)
                    xis
               end for effect in GaPSE.GR_EFFECTS_LDxGNC
          ]

          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(ss_quad, calc_ss_quad)])
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_sums_quad, calc_sums_quad)])
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[1], calc_xis_quad[1])]) # auto_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[2], calc_xis_quad[2])]) # auto_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[3], calc_xis_quad[3])]) # auto_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[4], calc_xis_quad[4])]) # auto_integrated
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[5], calc_xis_quad[5])]) # lensing_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[6], calc_xis_quad[6])]) # doppler_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[7], calc_xis_quad[7])]) # doppler_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[8], calc_xis_quad[8])]) # localgp_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[9], calc_xis_quad[9])]) # doppler_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[10], calc_xis_quad[10])]) # integratedgp_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[11], calc_xis_quad[11])]) # lensing_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[12], calc_xis_quad[12])]) # localgp_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[13], calc_xis_quad[13])]) # lensing_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[14], calc_xis_quad[14])]) # integratedgp_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[15], calc_xis_quad[15])]) # localgp_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[16], calc_xis_quad[16])]) # integratedgp_localgp

          rm(name)
          for effect in GaPSE.GR_EFFECTS_LDxGNC
               rm("all_LDxGNC_standalones_CFs/xi_LDxGNC_" * effect * "_L$L" * ".txt")
          end
          rm("all_LDxGNC_standalones_CFs")

          #### trap ####
          ss_trap, res_sums_trap, res_xis_trap = GaPSE.readxyall(
               "datatest/LDxGNC_SumXiMultipoles/xis_LDxGNC_L$L" * "_trap_noF_specific_ss.txt", comments=true)

          name = "calc_xis_LDxGNC_L$L" * "_trap_noF_specific_ss.txt"
          isfile(name) && rm(name)

          GaPSE.print_map_sum_ξ_LDxGNC_multipole(COSMO, name, ss_trap; 
               single=false, L=L, alg=:trap, kwargs...)
          calc_ss_trap, calc_sums_trap =  GaPSE.readxy(name, comments=true)
          calc_xis_trap = [
               begin
                    a_name = "all_LDxGNC_standalones_CFs/xi_LDxGNC_" * effect * "_L$L" * ".txt"
                    _, xis = GaPSE.readxy(a_name, comments = true)
                    xis
               end for effect in GaPSE.GR_EFFECTS_LDxGNC
          ]

          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(ss_trap, calc_ss_trap)])
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_sums_trap, calc_sums_trap)])
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[1], calc_xis_trap[1])]) # auto_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[2], calc_xis_trap[2])]) # auto_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[3], calc_xis_trap[3])]) # auto_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[4], calc_xis_trap[4])]) # auto_integrated
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[5], calc_xis_trap[5])]) # lensing_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[6], calc_xis_trap[6])]) # doppler_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[7], calc_xis_trap[7])]) # doppler_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[8], calc_xis_trap[8])]) # localgp_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[9], calc_xis_trap[9])]) # doppler_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[10], calc_xis_trap[10])]) # integratedgp_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[11], calc_xis_trap[11])]) # lensing_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[12], calc_xis_trap[12])]) # localgp_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[13], calc_xis_trap[13])]) # lensing_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[14], calc_xis_trap[14])]) # integratedgp_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[15], calc_xis_trap[15])]) # localgp_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[16], calc_xis_trap[16])]) # integratedgp_localgp

          rm(name)
          for effect in GaPSE.GR_EFFECTS_LDxGNC
               rm("all_LDxGNC_standalones_CFs/xi_LDxGNC_" * effect * "_L$L" * ".txt")
          end
          rm("all_LDxGNC_standalones_CFs")

          #### lob ####
          ss_lob, res_sums_lob, res_xis_lob = GaPSE.readxyall(
               "datatest/LDxGNC_SumXiMultipoles/xis_LDxGNC_L$L" * "_lob_noF_specific_ss.txt", comments=true)

          name = "calc_xis_LDxGNC_L$L" * "_lob_noF_specific_ss.txt"
          isfile(name) && rm(name)

          GaPSE.print_map_sum_ξ_LDxGNC_multipole(COSMO, name, ss_lob; 
               single=false, L=L, alg=:lobatto, kwargs...)
          calc_ss_lob, calc_sums_lob =  GaPSE.readxy(name, comments=true)
          calc_xis_lob = [
               begin
                    a_name = "all_LDxGNC_standalones_CFs/xi_LDxGNC_" * effect * "_L$L" * ".txt"
                    _, xis = GaPSE.readxy(a_name, comments = true)
                    xis
               end for effect in GaPSE.GR_EFFECTS_LDxGNC
          ]

          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(ss_lob, calc_ss_lob)])
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_sums_lob, calc_sums_lob)])
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[1], calc_xis_lob[1])]) # auto_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[2], calc_xis_lob[2])]) # auto_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[3], calc_xis_lob[3])]) # auto_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[4], calc_xis_lob[4])]) # auto_integrated
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[5], calc_xis_lob[5])]) # lensing_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[6], calc_xis_lob[6])]) # doppler_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[7], calc_xis_lob[7])]) # doppler_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[8], calc_xis_lob[8])]) # localgp_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[9], calc_xis_lob[9])]) # doppler_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[10], calc_xis_lob[10])]) # integratedgp_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[11], calc_xis_lob[11])]) # lensing_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[12], calc_xis_lob[12])]) # localgp_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[13], calc_xis_lob[13])]) # lensing_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[14], calc_xis_lob[14])]) # integratedgp_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[15], calc_xis_lob[15])]) # localgp_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[16], calc_xis_lob[16])]) # integratedgp_localgp

          rm(name)
          for effect in GaPSE.GR_EFFECTS_LDxGNC
               rm("all_LDxGNC_standalones_CFs/xi_LDxGNC_" * effect * "_L$L" * ".txt")
          end
          rm("all_LDxGNC_standalones_CFs")

     end

     @testset "L = 3 no_window" begin
          L = 3

          #### trap ####
          ss_trap, res_sums_trap, res_xis_trap = GaPSE.readxyall(
               "datatest/LDxGNC_SumXiMultipoles/xis_LDxGNC_L$L" * "_trap_noF_specific_ss.txt", comments=true)

          name = "calc_xis_LDxGNC_L$L" * "_trap_noF_specific_ss.txt"
          isfile(name) && rm(name)

          GaPSE.print_map_sum_ξ_LDxGNC_multipole(COSMO, name, ss_trap; 
               single=false, L=L, alg=:trap, kwargs...)
          calc_ss_trap, calc_sums_trap =  GaPSE.readxy(name, comments=true)
          calc_xis_trap = [
               begin
                    a_name = "all_LDxGNC_standalones_CFs/xi_LDxGNC_" * effect * "_L$L" * ".txt"
                    _, xis = GaPSE.readxy(a_name, comments = true)
                    xis
               end for effect in GaPSE.GR_EFFECTS_LDxGNC
          ]

          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(ss_trap, calc_ss_trap)])
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_sums_trap, calc_sums_trap)])
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[1], calc_xis_trap[1])]) # auto_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[2], calc_xis_trap[2])]) # auto_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[3], calc_xis_trap[3])]) # auto_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[4], calc_xis_trap[4])]) # auto_integrated
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[5], calc_xis_trap[5])]) # lensing_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[6], calc_xis_trap[6])]) # doppler_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[7], calc_xis_trap[7])]) # doppler_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[8], calc_xis_trap[8])]) # localgp_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[9], calc_xis_trap[9])]) # doppler_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[10], calc_xis_trap[10])]) # integratedgp_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[11], calc_xis_trap[11])]) # lensing_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[12], calc_xis_trap[12])]) # localgp_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[13], calc_xis_trap[13])]) # lensing_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[14], calc_xis_trap[14])]) # integratedgp_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[15], calc_xis_trap[15])]) # localgp_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[16], calc_xis_trap[16])]) # integratedgp_localgp

          rm(name)
          for effect in GaPSE.GR_EFFECTS_LDxGNC
               rm("all_LDxGNC_standalones_CFs/xi_LDxGNC_" * effect * "_L$L" * ".txt")
          end
          rm("all_LDxGNC_standalones_CFs")

          #### lob ####
          ss_lob, res_sums_lob, res_xis_lob = GaPSE.readxyall(
               "datatest/LDxGNC_SumXiMultipoles/xis_LDxGNC_L$L" * "_lob_noF_specific_ss.txt", comments=true)

          name = "calc_xis_LDxGNC_L$L" * "_lob_noF_specific_ss.txt"
          isfile(name) && rm(name)

          GaPSE.print_map_sum_ξ_LDxGNC_multipole(COSMO, name, ss_lob; 
               single=false, L=L, alg=:lobatto, kwargs...)
          calc_ss_lob, calc_sums_lob =  GaPSE.readxy(name, comments=true)
          calc_xis_lob = [
               begin
                    a_name = "all_LDxGNC_standalones_CFs/xi_LDxGNC_" * effect * "_L$L" * ".txt"
                    _, xis = GaPSE.readxy(a_name, comments = true)
                    xis
               end for effect in GaPSE.GR_EFFECTS_LDxGNC
          ]

          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(ss_lob, calc_ss_lob)])
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_sums_lob, calc_sums_lob)])
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[1], calc_xis_lob[1])]) # auto_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[2], calc_xis_lob[2])]) # auto_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[3], calc_xis_lob[3])]) # auto_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[4], calc_xis_lob[4])]) # auto_integrated
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[5], calc_xis_lob[5])]) # lensing_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[6], calc_xis_lob[6])]) # doppler_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[7], calc_xis_lob[7])]) # doppler_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[8], calc_xis_lob[8])]) # localgp_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[9], calc_xis_lob[9])]) # doppler_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[10], calc_xis_lob[10])]) # integratedgp_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[11], calc_xis_lob[11])]) # lensing_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[12], calc_xis_lob[12])]) # localgp_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[13], calc_xis_lob[13])]) # lensing_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[14], calc_xis_lob[14])]) # integratedgp_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[15], calc_xis_lob[15])]) # localgp_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[16], calc_xis_lob[16])]) # integratedgp_localgp

          rm(name)
          for effect in GaPSE.GR_EFFECTS_LDxGNC
               rm("all_LDxGNC_standalones_CFs/xi_LDxGNC_" * effect * "_L$L" * ".txt")
          end
          rm("all_LDxGNC_standalones_CFs")

     end

     @testset "L = 4 no_window" begin
          L = 4

          #### trap ####
          ss_trap, res_sums_trap, res_xis_trap = GaPSE.readxyall(
               "datatest/LDxGNC_SumXiMultipoles/xis_LDxGNC_L$L" * "_trap_noF_specific_ss.txt", comments=true)

          name = "calc_xis_LDxGNC_L$L" * "_trap_noF_specific_ss.txt"
          isfile(name) && rm(name)

          GaPSE.print_map_sum_ξ_LDxGNC_multipole(COSMO, name, ss_trap; 
               single=false, L=L, alg=:trap, kwargs...)
          calc_ss_trap, calc_sums_trap =  GaPSE.readxy(name, comments=true)
          calc_xis_trap = [
               begin
                    a_name = "all_LDxGNC_standalones_CFs/xi_LDxGNC_" * effect * "_L$L" * ".txt"
                    _, xis = GaPSE.readxy(a_name, comments = true)
                    xis
               end for effect in GaPSE.GR_EFFECTS_LDxGNC
          ]

          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(ss_trap, calc_ss_trap)])
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_sums_trap, calc_sums_trap)])
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[1], calc_xis_trap[1])]) # auto_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[2], calc_xis_trap[2])]) # auto_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[3], calc_xis_trap[3])]) # auto_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[4], calc_xis_trap[4])]) # auto_integrated
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[5], calc_xis_trap[5])]) # lensing_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[6], calc_xis_trap[6])]) # doppler_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[7], calc_xis_trap[7])]) # doppler_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[8], calc_xis_trap[8])]) # localgp_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[9], calc_xis_trap[9])]) # doppler_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[10], calc_xis_trap[10])]) # integratedgp_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[11], calc_xis_trap[11])]) # lensing_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[12], calc_xis_trap[12])]) # localgp_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[13], calc_xis_trap[13])]) # lensing_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[14], calc_xis_trap[14])]) # integratedgp_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[15], calc_xis_trap[15])]) # localgp_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[16], calc_xis_trap[16])]) # integratedgp_localgp

          rm(name)
          for effect in GaPSE.GR_EFFECTS_LDxGNC
               rm("all_LDxGNC_standalones_CFs/xi_LDxGNC_" * effect * "_L$L" * ".txt")
          end
          rm("all_LDxGNC_standalones_CFs")

          #### lob ####
          ss_lob, res_sums_lob, res_xis_lob = GaPSE.readxyall(
               "datatest/LDxGNC_SumXiMultipoles/xis_LDxGNC_L$L" * "_lob_noF_specific_ss.txt", comments=true)

          name = "calc_xis_LDxGNC_L$L" * "_lob_noF_specific_ss.txt"
          isfile(name) && rm(name)

          GaPSE.print_map_sum_ξ_LDxGNC_multipole(COSMO, name, ss_lob; 
               single=false, L=L, alg=:lobatto, kwargs...)
          calc_ss_lob, calc_sums_lob =  GaPSE.readxy(name, comments=true)
          calc_xis_lob = [
               begin
                    a_name = "all_LDxGNC_standalones_CFs/xi_LDxGNC_" * effect * "_L$L" * ".txt"
                    _, xis = GaPSE.readxy(a_name, comments = true)
                    xis
               end for effect in GaPSE.GR_EFFECTS_LDxGNC
          ]

          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(ss_lob, calc_ss_lob)])
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_sums_lob, calc_sums_lob)])
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[1], calc_xis_lob[1])]) # auto_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[2], calc_xis_lob[2])]) # auto_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[3], calc_xis_lob[3])]) # auto_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[4], calc_xis_lob[4])]) # auto_integrated
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[5], calc_xis_lob[5])]) # lensing_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[6], calc_xis_lob[6])]) # doppler_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[7], calc_xis_lob[7])]) # doppler_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[8], calc_xis_lob[8])]) # localgp_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[9], calc_xis_lob[9])]) # doppler_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[10], calc_xis_lob[10])]) # integratedgp_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[11], calc_xis_lob[11])]) # lensing_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[12], calc_xis_lob[12])]) # localgp_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[13], calc_xis_lob[13])]) # lensing_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[14], calc_xis_lob[14])]) # integratedgp_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[15], calc_xis_lob[15])]) # localgp_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[16], calc_xis_lob[16])]) # integratedgp_localgp

          rm(name)
          for effect in GaPSE.GR_EFFECTS_LDxGNC
               rm("all_LDxGNC_standalones_CFs/xi_LDxGNC_" * effect * "_L$L" * ".txt")
          end
          rm("all_LDxGNC_standalones_CFs")

     end
end

println("\nAlmost finished...")


@testset "test print_map_sum_ξ_LDxGNC_multipole with_window second method" begin
     RTOL = 2e-2

     kwargs = Dict(
          :use_windows => true,
          :enhancer => 1e8, :s1 => nothing,
          :N_trap => 30, :N_lob => 30,
          :atol_quad => 0.0, :rtol_quad => 1e-2,
          :N_χs => 40, :N_χs_2 => 20,
          :pr => false,
     )

     @testset "test_zeros" begin
          @test_throws AssertionError GaPSE.GaPSE.print_map_sum_ξ_LDxGNC_multipole(COSMO, "NONEXISTINGDIR/file.txt", [1.0, 2.0, 3.0])
     end

     @testset "L = 0 with_window" begin
          L = 0

          ##### quad ####
          ss_quad, res_sums_quad, res_xis_quad = GaPSE.readxyall(
               "datatest/LDxGNC_SumXiMultipoles/xis_LDxGNC_L$L" * "_quad_withF_specific_ss.txt", comments=true)

          name = "calc_xis_LDxGNC_L$L" * "_quad_withF_specific_ss.txt"
          isfile(name) && rm(name)

          GaPSE.print_map_sum_ξ_LDxGNC_multipole(COSMO, name, ss_quad; 
               single=false, L=L, alg=:quad, kwargs...)

          calc_ss_quad, calc_sums_quad =  GaPSE.readxy(name, comments=true)
          calc_xis_quad = [
               begin
                    a_name = "all_LDxGNC_standalones_CFs/xi_LDxGNC_" * effect * "_L$L" * ".txt"
                    _, xis = GaPSE.readxy(a_name, comments = true)
                    xis
               end for effect in GaPSE.GR_EFFECTS_LDxGNC
          ]

          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(ss_quad, calc_ss_quad)])
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_sums_quad, calc_sums_quad)])
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[1], calc_xis_quad[1])]) # auto_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[2], calc_xis_quad[2])]) # auto_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[3], calc_xis_quad[3])]) # auto_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[4], calc_xis_quad[4])]) # auto_integrated
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[5], calc_xis_quad[5])]) # lensing_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[6], calc_xis_quad[6])]) # doppler_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[7], calc_xis_quad[7])]) # doppler_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[8], calc_xis_quad[8])]) # localgp_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[9], calc_xis_quad[9])]) # doppler_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[10], calc_xis_quad[10])]) # integratedgp_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[11], calc_xis_quad[11])]) # lensing_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[12], calc_xis_quad[12])]) # localgp_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[13], calc_xis_quad[13])]) # lensing_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[14], calc_xis_quad[14])]) # integratedgp_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[15], calc_xis_quad[15])]) # localgp_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[16], calc_xis_quad[16])]) # integratedgp_localgp

          rm(name)
          for effect in GaPSE.GR_EFFECTS_LDxGNC
               rm("all_LDxGNC_standalones_CFs/xi_LDxGNC_" * effect * "_L$L" * ".txt")
          end
          rm("all_LDxGNC_standalones_CFs")

          #### trap ####
          ss_trap, res_sums_trap, res_xis_trap = GaPSE.readxyall(
               "datatest/LDxGNC_SumXiMultipoles/xis_LDxGNC_L$L" * "_trap_withF_specific_ss.txt", comments=true)

          name = "calc_xis_LDxGNC_L$L" * "_trap_withF_specific_ss.txt"
          isfile(name) && rm(name)

          GaPSE.print_map_sum_ξ_LDxGNC_multipole(COSMO, name, ss_trap; 
               single=false, L=L, alg=:trap, kwargs...)
          calc_ss_trap, calc_sums_trap =  GaPSE.readxy(name, comments=true)
          calc_xis_trap = [
               begin
                    a_name = "all_LDxGNC_standalones_CFs/xi_LDxGNC_" * effect * "_L$L" * ".txt"
                    _, xis = GaPSE.readxy(a_name, comments = true)
                    xis
               end for effect in GaPSE.GR_EFFECTS_LDxGNC
          ]

          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(ss_trap, calc_ss_trap)])
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_sums_trap, calc_sums_trap)])
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[1], calc_xis_trap[1])]) # auto_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[2], calc_xis_trap[2])]) # auto_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[3], calc_xis_trap[3])]) # auto_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[4], calc_xis_trap[4])]) # auto_integrated
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[5], calc_xis_trap[5])]) # lensing_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[6], calc_xis_trap[6])]) # doppler_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[7], calc_xis_trap[7])]) # doppler_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[8], calc_xis_trap[8])]) # localgp_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[9], calc_xis_trap[9])]) # doppler_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[10], calc_xis_trap[10])]) # integratedgp_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[11], calc_xis_trap[11])]) # lensing_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[12], calc_xis_trap[12])]) # localgp_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[13], calc_xis_trap[13])]) # lensing_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[14], calc_xis_trap[14])]) # integratedgp_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[15], calc_xis_trap[15])]) # localgp_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[16], calc_xis_trap[16])]) # integratedgp_localgp

          rm(name)
          for effect in GaPSE.GR_EFFECTS_LDxGNC
               rm("all_LDxGNC_standalones_CFs/xi_LDxGNC_" * effect * "_L$L" * ".txt")
          end
          rm("all_LDxGNC_standalones_CFs")

          #### lob ####
          ss_lob, res_sums_lob, res_xis_lob = GaPSE.readxyall(
               "datatest/LDxGNC_SumXiMultipoles/xis_LDxGNC_L$L" * "_lob_withF_specific_ss.txt", comments=true)

          name = "calc_xis_LDxGNC_L$L" * "_lob_withF_specific_ss.txt"
          isfile(name) && rm(name)

          GaPSE.print_map_sum_ξ_LDxGNC_multipole(COSMO, name, ss_lob; 
               single=false, L=L, alg=:lobatto, kwargs...)
          calc_ss_lob, calc_sums_lob, calc_xis_lob =  GaPSE.readxyall(name,  comments=true)
          calc_ss_lob, calc_sums_lob =  GaPSE.readxy(name, comments=true)
          calc_xis_lob = [
               begin
                    a_name = "all_LDxGNC_standalones_CFs/xi_LDxGNC_" * effect * "_L$L" * ".txt"
                    _, xis = GaPSE.readxy(a_name, comments = true)
                    xis
               end for effect in GaPSE.GR_EFFECTS_LDxGNC
          ]

          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(ss_lob, calc_ss_lob)])
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_sums_lob, calc_sums_lob)])
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[1], calc_xis_lob[1])]) # auto_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[2], calc_xis_lob[2])]) # auto_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[3], calc_xis_lob[3])]) # auto_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[4], calc_xis_lob[4])]) # auto_integrated
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[5], calc_xis_lob[5])]) # lensing_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[6], calc_xis_lob[6])]) # doppler_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[7], calc_xis_lob[7])]) # doppler_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[8], calc_xis_lob[8])]) # localgp_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[9], calc_xis_lob[9])]) # doppler_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[10], calc_xis_lob[10])]) # integratedgp_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[11], calc_xis_lob[11])]) # lensing_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[12], calc_xis_lob[12])]) # localgp_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[13], calc_xis_lob[13])]) # lensing_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[14], calc_xis_lob[14])]) # integratedgp_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[15], calc_xis_lob[15])]) # localgp_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[16], calc_xis_lob[16])]) # integratedgp_localgp

          rm(name)
          for effect in GaPSE.GR_EFFECTS_LDxGNC
               rm("all_LDxGNC_standalones_CFs/xi_LDxGNC_" * effect * "_L$L" * ".txt")
          end
          rm("all_LDxGNC_standalones_CFs")

     end

     @testset "L = 1 with_window" begin
          L = 1

          #### trap ####
          ss_trap, res_sums_trap, res_xis_trap = GaPSE.readxyall(
               "datatest/LDxGNC_SumXiMultipoles/xis_LDxGNC_L$L" * "_trap_withF_specific_ss.txt", comments=true)

          name = "calc_xis_LDxGNC_L$L" * "_trap_withF_specific_ss.txt"
          isfile(name) && rm(name)

          GaPSE.print_map_sum_ξ_LDxGNC_multipole(COSMO, name, ss_trap; 
               single=false, L=L, alg=:trap, kwargs...)
          calc_ss_trap, calc_sums_trap =  GaPSE.readxy(name, comments=true)
          calc_xis_trap = [
               begin
                    a_name = "all_LDxGNC_standalones_CFs/xi_LDxGNC_" * effect * "_L$L" * ".txt"
                    _, xis = GaPSE.readxy(a_name, comments = true)
                    xis
               end for effect in GaPSE.GR_EFFECTS_LDxGNC
          ]

          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(ss_trap, calc_ss_trap)])
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_sums_trap, calc_sums_trap)])
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[1], calc_xis_trap[1])]) # auto_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[2], calc_xis_trap[2])]) # auto_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[3], calc_xis_trap[3])]) # auto_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[4], calc_xis_trap[4])]) # auto_integrated
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[5], calc_xis_trap[5])]) # lensing_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[6], calc_xis_trap[6])]) # doppler_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[7], calc_xis_trap[7])]) # doppler_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[8], calc_xis_trap[8])]) # localgp_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[9], calc_xis_trap[9])]) # doppler_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[10], calc_xis_trap[10])]) # integratedgp_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[11], calc_xis_trap[11])]) # lensing_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[12], calc_xis_trap[12])]) # localgp_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[13], calc_xis_trap[13])]) # lensing_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[14], calc_xis_trap[14])]) # integratedgp_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[15], calc_xis_trap[15])]) # localgp_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[16], calc_xis_trap[16])]) # integratedgp_localgp

          rm(name)
          for effect in GaPSE.GR_EFFECTS_LDxGNC
               rm("all_LDxGNC_standalones_CFs/xi_LDxGNC_" * effect * "_L$L" * ".txt")
          end
          rm("all_LDxGNC_standalones_CFs")

          #### lob ####
          ss_lob, res_sums_lob, res_xis_lob = GaPSE.readxyall(
               "datatest/LDxGNC_SumXiMultipoles/xis_LDxGNC_L$L" * "_lob_withF_specific_ss.txt", comments=true)

          name = "calc_xis_LDxGNC_L$L" * "_lob_withF_specific_ss.txt"
          isfile(name) && rm(name)

          GaPSE.print_map_sum_ξ_LDxGNC_multipole(COSMO, name, ss_lob; 
               single=false, L=L, alg=:lobatto, kwargs...)
          calc_ss_lob, calc_sums_lob =  GaPSE.readxy(name, comments=true)
          calc_xis_lob = [
               begin
                    a_name = "all_LDxGNC_standalones_CFs/xi_LDxGNC_" * effect * "_L$L" * ".txt"
                    _, xis = GaPSE.readxy(a_name, comments = true)
                    xis
               end for effect in GaPSE.GR_EFFECTS_LDxGNC
          ]

          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(ss_lob, calc_ss_lob)])
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_sums_lob, calc_sums_lob)])
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[1], calc_xis_lob[1])]) # auto_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[2], calc_xis_lob[2])]) # auto_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[3], calc_xis_lob[3])]) # auto_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[4], calc_xis_lob[4])]) # auto_integrated
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[5], calc_xis_lob[5])]) # lensing_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[6], calc_xis_lob[6])]) # doppler_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[7], calc_xis_lob[7])]) # doppler_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[8], calc_xis_lob[8])]) # localgp_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[9], calc_xis_lob[9])]) # doppler_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[10], calc_xis_lob[10])]) # integratedgp_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[11], calc_xis_lob[11])]) # lensing_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[12], calc_xis_lob[12])]) # localgp_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[13], calc_xis_lob[13])]) # lensing_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[14], calc_xis_lob[14])]) # integratedgp_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[15], calc_xis_lob[15])]) # localgp_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[16], calc_xis_lob[16])]) # integratedgp_localgp

          rm(name)
          for effect in GaPSE.GR_EFFECTS_LDxGNC
               rm("all_LDxGNC_standalones_CFs/xi_LDxGNC_" * effect * "_L$L" * ".txt")
          end
          rm("all_LDxGNC_standalones_CFs")

     end

     @testset "L = 2 with_window" begin
          L = 2

          ##### quad ####
          ss_quad, res_sums_quad, res_xis_quad = GaPSE.readxyall(
               "datatest/LDxGNC_SumXiMultipoles/xis_LDxGNC_L$L" * "_quad_withF_specific_ss.txt", comments=true)

          name = "calc_xis_LDxGNC_L$L" * "_quad_withF_specific_ss.txt"
          isfile(name) && rm(name)

          GaPSE.print_map_sum_ξ_LDxGNC_multipole(COSMO, name, ss_quad; 
               single=false, L=L, alg=:quad, kwargs...)
          calc_ss_quad, calc_sums_quad =  GaPSE.readxy(name, comments=true)
          calc_xis_quad = [
               begin
                    a_name = "all_LDxGNC_standalones_CFs/xi_LDxGNC_" * effect * "_L$L" * ".txt"
                    _, xis = GaPSE.readxy(a_name, comments = true)
                    xis
               end for effect in GaPSE.GR_EFFECTS_LDxGNC
          ]

          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(ss_quad, calc_ss_quad)])
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_sums_quad, calc_sums_quad)])
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[1], calc_xis_quad[1])]) # auto_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[2], calc_xis_quad[2])]) # auto_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[3], calc_xis_quad[3])]) # auto_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[4], calc_xis_quad[4])]) # auto_integrated
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[5], calc_xis_quad[5])]) # lensing_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[6], calc_xis_quad[6])]) # doppler_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[7], calc_xis_quad[7])]) # doppler_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[8], calc_xis_quad[8])]) # localgp_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[9], calc_xis_quad[9])]) # doppler_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[10], calc_xis_quad[10])]) # integratedgp_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[11], calc_xis_quad[11])]) # lensing_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[12], calc_xis_quad[12])]) # localgp_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[13], calc_xis_quad[13])]) # lensing_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[14], calc_xis_quad[14])]) # integratedgp_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[15], calc_xis_quad[15])]) # localgp_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_quad[16], calc_xis_quad[16])]) # integratedgp_localgp

          rm(name)
          for effect in GaPSE.GR_EFFECTS_LDxGNC
               rm("all_LDxGNC_standalones_CFs/xi_LDxGNC_" * effect * "_L$L" * ".txt")
          end
          rm("all_LDxGNC_standalones_CFs")

          #### trap ####
          ss_trap, res_sums_trap, res_xis_trap = GaPSE.readxyall(
               "datatest/LDxGNC_SumXiMultipoles/xis_LDxGNC_L$L" * "_trap_withF_specific_ss.txt", comments=true)

          name = "calc_xis_LDxGNC_L$L" * "_trap_withF_specific_ss.txt"
          isfile(name) && rm(name)

          GaPSE.print_map_sum_ξ_LDxGNC_multipole(COSMO, name, ss_trap; 
               single=false, L=L, alg=:trap, kwargs...)
          calc_ss_trap, calc_sums_trap =  GaPSE.readxy(name, comments=true)
          calc_xis_trap = [
               begin
                    a_name = "all_LDxGNC_standalones_CFs/xi_LDxGNC_" * effect * "_L$L" * ".txt"
                    _, xis = GaPSE.readxy(a_name, comments = true)
                    xis
               end for effect in GaPSE.GR_EFFECTS_LDxGNC
          ]

          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(ss_trap, calc_ss_trap)])
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_sums_trap, calc_sums_trap)])
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[1], calc_xis_trap[1])]) # auto_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[2], calc_xis_trap[2])]) # auto_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[3], calc_xis_trap[3])]) # auto_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[4], calc_xis_trap[4])]) # auto_integrated
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[5], calc_xis_trap[5])]) # lensing_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[6], calc_xis_trap[6])]) # doppler_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[7], calc_xis_trap[7])]) # doppler_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[8], calc_xis_trap[8])]) # localgp_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[9], calc_xis_trap[9])]) # doppler_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[10], calc_xis_trap[10])]) # integratedgp_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[11], calc_xis_trap[11])]) # lensing_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[12], calc_xis_trap[12])]) # localgp_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[13], calc_xis_trap[13])]) # lensing_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[14], calc_xis_trap[14])]) # integratedgp_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[15], calc_xis_trap[15])]) # localgp_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[16], calc_xis_trap[16])]) # integratedgp_localgp

          rm(name)
          for effect in GaPSE.GR_EFFECTS_LDxGNC
               rm("all_LDxGNC_standalones_CFs/xi_LDxGNC_" * effect * "_L$L" * ".txt")
          end
          rm("all_LDxGNC_standalones_CFs")

          #### lob ####
          ss_lob, res_sums_lob, res_xis_lob = GaPSE.readxyall(
               "datatest/LDxGNC_SumXiMultipoles/xis_LDxGNC_L$L" * "_lob_withF_specific_ss.txt", comments=true)

          name = "calc_xis_LDxGNC_L$L" * "_lob_withF_specific_ss.txt"
          isfile(name) && rm(name)

          GaPSE.print_map_sum_ξ_LDxGNC_multipole(COSMO, name, ss_lob; 
               single=false, L=L, alg=:lobatto, kwargs...)
          calc_ss_lob, calc_sums_lob =  GaPSE.readxy(name, comments=true)
          calc_xis_lob = [
               begin
                    a_name = "all_LDxGNC_standalones_CFs/xi_LDxGNC_" * effect * "_L$L" * ".txt"
                    _, xis = GaPSE.readxy(a_name, comments = true)
                    xis
               end for effect in GaPSE.GR_EFFECTS_LDxGNC
          ]

          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(ss_lob, calc_ss_lob)])
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_sums_lob, calc_sums_lob)])
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[1], calc_xis_lob[1])]) # auto_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[2], calc_xis_lob[2])]) # auto_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[3], calc_xis_lob[3])]) # auto_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[4], calc_xis_lob[4])]) # auto_integrated
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[5], calc_xis_lob[5])]) # lensing_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[6], calc_xis_lob[6])]) # doppler_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[7], calc_xis_lob[7])]) # doppler_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[8], calc_xis_lob[8])]) # localgp_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[9], calc_xis_lob[9])]) # doppler_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[10], calc_xis_lob[10])]) # integratedgp_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[11], calc_xis_lob[11])]) # lensing_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[12], calc_xis_lob[12])]) # localgp_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[13], calc_xis_lob[13])]) # lensing_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[14], calc_xis_lob[14])]) # integratedgp_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[15], calc_xis_lob[15])]) # localgp_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[16], calc_xis_lob[16])]) # integratedgp_localgp

          rm(name)
          for effect in GaPSE.GR_EFFECTS_LDxGNC
               rm("all_LDxGNC_standalones_CFs/xi_LDxGNC_" * effect * "_L$L" * ".txt")
          end
          rm("all_LDxGNC_standalones_CFs")

     end

     @testset "L = 3 with_window" begin
          L = 3

          #### trap ####
          ss_trap, res_sums_trap, res_xis_trap = GaPSE.readxyall(
               "datatest/LDxGNC_SumXiMultipoles/xis_LDxGNC_L$L" * "_trap_withF_specific_ss.txt", comments=true)

          name = "calc_xis_LDxGNC_L$L" * "_trap_withF_specific_ss.txt"
          isfile(name) && rm(name)

          GaPSE.print_map_sum_ξ_LDxGNC_multipole(COSMO, name, ss_trap; 
               single=false, L=L, alg=:trap, kwargs...)
          calc_ss_trap, calc_sums_trap =  GaPSE.readxy(name, comments=true)
          calc_xis_trap = [
               begin
                    a_name = "all_LDxGNC_standalones_CFs/xi_LDxGNC_" * effect * "_L$L" * ".txt"
                    _, xis = GaPSE.readxy(a_name, comments = true)
                    xis
               end for effect in GaPSE.GR_EFFECTS_LDxGNC
          ]

          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(ss_trap, calc_ss_trap)])
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_sums_trap, calc_sums_trap)])
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[1], calc_xis_trap[1])]) # auto_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[2], calc_xis_trap[2])]) # auto_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[3], calc_xis_trap[3])]) # auto_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[4], calc_xis_trap[4])]) # auto_integrated
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[5], calc_xis_trap[5])]) # lensing_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[6], calc_xis_trap[6])]) # doppler_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[7], calc_xis_trap[7])]) # doppler_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[8], calc_xis_trap[8])]) # localgp_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[9], calc_xis_trap[9])]) # doppler_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[10], calc_xis_trap[10])]) # integratedgp_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[11], calc_xis_trap[11])]) # lensing_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[12], calc_xis_trap[12])]) # localgp_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[13], calc_xis_trap[13])]) # lensing_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[14], calc_xis_trap[14])]) # integratedgp_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[15], calc_xis_trap[15])]) # localgp_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[16], calc_xis_trap[16])]) # integratedgp_localgp

          rm(name)
          for effect in GaPSE.GR_EFFECTS_LDxGNC
               rm("all_LDxGNC_standalones_CFs/xi_LDxGNC_" * effect * "_L$L" * ".txt")
          end
          rm("all_LDxGNC_standalones_CFs")

          #### lob ####
          ss_lob, res_sums_lob, res_xis_lob = GaPSE.readxyall(
               "datatest/LDxGNC_SumXiMultipoles/xis_LDxGNC_L$L" * "_lob_withF_specific_ss.txt", comments=true)

          name = "calc_xis_LDxGNC_L$L" * "_lob_withF_specific_ss.txt"
          isfile(name) && rm(name)

          GaPSE.print_map_sum_ξ_LDxGNC_multipole(COSMO, name, ss_lob; 
               single=false, L=L, alg=:lobatto, kwargs...)
          calc_ss_lob, calc_sums_lob =  GaPSE.readxy(name, comments=true)
          calc_xis_lob = [
               begin
                    a_name = "all_LDxGNC_standalones_CFs/xi_LDxGNC_" * effect * "_L$L" * ".txt"
                    _, xis = GaPSE.readxy(a_name, comments = true)
                    xis
               end for effect in GaPSE.GR_EFFECTS_LDxGNC
          ]

          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(ss_lob, calc_ss_lob)])
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_sums_lob, calc_sums_lob)])
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[1], calc_xis_lob[1])]) # auto_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[2], calc_xis_lob[2])]) # auto_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[3], calc_xis_lob[3])]) # auto_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[4], calc_xis_lob[4])]) # auto_integrated
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[5], calc_xis_lob[5])]) # lensing_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[6], calc_xis_lob[6])]) # doppler_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[7], calc_xis_lob[7])]) # doppler_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[8], calc_xis_lob[8])]) # localgp_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[9], calc_xis_lob[9])]) # doppler_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[10], calc_xis_lob[10])]) # integratedgp_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[11], calc_xis_lob[11])]) # lensing_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[12], calc_xis_lob[12])]) # localgp_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[13], calc_xis_lob[13])]) # lensing_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[14], calc_xis_lob[14])]) # integratedgp_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[15], calc_xis_lob[15])]) # localgp_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[16], calc_xis_lob[16])]) # integratedgp_localgp

          rm(name)
          for effect in GaPSE.GR_EFFECTS_LDxGNC
               rm("all_LDxGNC_standalones_CFs/xi_LDxGNC_" * effect * "_L$L" * ".txt")
          end
          rm("all_LDxGNC_standalones_CFs")

     end

     @testset "L = 4 with_window" begin
          L = 4

          #### trap ####
          ss_trap, res_sums_trap, res_xis_trap = GaPSE.readxyall(
               "datatest/LDxGNC_SumXiMultipoles/xis_LDxGNC_L$L" * "_trap_withF_specific_ss.txt", comments=true)

          name = "calc_xis_LDxGNC_L$L" * "_trap_withF_specific_ss.txt"
          isfile(name) && rm(name)

          GaPSE.print_map_sum_ξ_LDxGNC_multipole(COSMO, name, ss_trap; 
               single=false, L=L, alg=:trap, kwargs...)
          calc_ss_trap, calc_sums_trap =  GaPSE.readxy(name, comments=true)
          calc_xis_trap = [
               begin
                    a_name = "all_LDxGNC_standalones_CFs/xi_LDxGNC_" * effect * "_L$L" * ".txt"
                    _, xis = GaPSE.readxy(a_name, comments = true)
                    xis
               end for effect in GaPSE.GR_EFFECTS_LDxGNC
          ]

          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(ss_trap, calc_ss_trap)])
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_sums_trap, calc_sums_trap)])
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[1], calc_xis_trap[1])]) # auto_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[2], calc_xis_trap[2])]) # auto_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[3], calc_xis_trap[3])]) # auto_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[4], calc_xis_trap[4])]) # auto_integrated
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[5], calc_xis_trap[5])]) # lensing_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[6], calc_xis_trap[6])]) # doppler_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[7], calc_xis_trap[7])]) # doppler_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[8], calc_xis_trap[8])]) # localgp_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[9], calc_xis_trap[9])]) # doppler_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[10], calc_xis_trap[10])]) # integratedgp_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[11], calc_xis_trap[11])]) # lensing_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[12], calc_xis_trap[12])]) # localgp_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[13], calc_xis_trap[13])]) # lensing_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[14], calc_xis_trap[14])]) # integratedgp_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[15], calc_xis_trap[15])]) # localgp_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_trap[16], calc_xis_trap[16])]) # integratedgp_localgp

          rm(name)
          for effect in GaPSE.GR_EFFECTS_LDxGNC
               rm("all_LDxGNC_standalones_CFs/xi_LDxGNC_" * effect * "_L$L" * ".txt")
          end
          rm("all_LDxGNC_standalones_CFs")

          #### lob ####
          ss_lob, res_sums_lob, res_xis_lob = GaPSE.readxyall(
               "datatest/LDxGNC_SumXiMultipoles/xis_LDxGNC_L$L" * "_lob_withF_specific_ss.txt", comments=true)

          name = "calc_xis_LDxGNC_L$L" * "_lob_withF_specific_ss.txt"
          isfile(name) && rm(name)

          GaPSE.print_map_sum_ξ_LDxGNC_multipole(COSMO, name, ss_lob; 
               single=false, L=L, alg=:lobatto, kwargs...)
          calc_ss_lob, calc_sums_lob =  GaPSE.readxy(name, comments=true)
          calc_xis_lob = [
               begin
                    a_name = "all_LDxGNC_standalones_CFs/xi_LDxGNC_" * effect * "_L$L" * ".txt"
                    _, xis = GaPSE.readxy(a_name, comments = true)
                    xis
               end for effect in GaPSE.GR_EFFECTS_LDxGNC
          ]

          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(ss_lob, calc_ss_lob)])
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_sums_lob, calc_sums_lob)])
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[1], calc_xis_lob[1])]) # auto_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[2], calc_xis_lob[2])]) # auto_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[3], calc_xis_lob[3])]) # auto_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[4], calc_xis_lob[4])]) # auto_integrated
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[5], calc_xis_lob[5])]) # lensing_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[6], calc_xis_lob[6])]) # doppler_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[7], calc_xis_lob[7])]) # doppler_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[8], calc_xis_lob[8])]) # localgp_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[9], calc_xis_lob[9])]) # doppler_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[10], calc_xis_lob[10])]) # integratedgp_doppler
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[11], calc_xis_lob[11])]) # lensing_localgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[12], calc_xis_lob[12])]) # localgp_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[13], calc_xis_lob[13])]) # lensing_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[14], calc_xis_lob[14])]) # integratedgp_lensing
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[15], calc_xis_lob[15])]) # localgp_integratedgp
          @test all([isapprox(a, r, rtol=RTOL) for (a, r) in zip(res_xis_lob[16], calc_xis_lob[16])]) # integratedgp_localgp

          rm(name)
          for effect in GaPSE.GR_EFFECTS_LDxGNC
               rm("all_LDxGNC_standalones_CFs/xi_LDxGNC_" * effect * "_L$L" * ".txt")
          end
          rm("all_LDxGNC_standalones_CFs")

     end
end


println("\nEnded tests on print_map_sum_ξ_LDxGNC_multipole!\n")
