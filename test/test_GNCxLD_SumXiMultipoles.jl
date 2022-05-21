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


@testset "test sum_ξ_GNCxLD_multipole no_window" begin
     table = readdlm("datatest/GNCxLD_SumXiMultipoles/map_sum_xi_GNCxLD_L0_noF_specific_ss.txt", comments = true)
     ss = convert(Vector{Float64}, table[:, 1])
     all_res_sums = convert(Vector{Float64}, table[:, 2])
     all_res_xis = [convert(Vector{Float64}, table[:, i]) for i in 3:22]
     RTOL = 1e-2

     kwargs = Dict(
          :L => 0, :use_windows => false,
          :enhancer => 1e8, :N_μs => 30,
          :μ_atol => 0.0, :μ_rtol => 1e-2,
          :N_χs => 100, :N_χs_2 => 50,
          :pr => true,
     )

     #=
     @testset "zeros" begin
          a_func(x) = x
          @test_throws AssertionError GaPSE.sum_ξ_GNCxLD_multipole(COSMO.s_eff, 10, COSMO;
               L = 0, use_windows = false, SPLINE = true, kwargs...)
          @test_throws AssertionError GaPSE.integral_on_mu(COSMO.s_eff, 10, a_func, COSMO;
               L = 0, use_windows = false, SPLINE = true, kwargs...)
     end
     =#

     @testset "s = 10, L = 0, no_window" begin
          s = 10

          ind = findfirst(x -> x≈s, ss)
          res_sum =  all_res_sums[ind]
          res_xis = [vec[ind] for vec in all_res_xis] 

          calc_res_sum, calc_res_xis = GaPSE.sum_ξ_GNCxLD_multipole(COSMO.s_eff, s, COSMO; kwargs...)

          @test isapprox(res_sum, calc_res_sum; rtol = RTOL)
          @test all([isapprox(a, r; rtol = RTOL) for (a, r) in zip(calc_res_xis, res_xis)])
     end

     @testset "s = 100, L = 0, no_window" begin
          s = 100

          ind = findfirst(x -> x≈s, ss)
          res_sum =  all_res_sums[ind]
          res_xis = [vec[ind] for vec in all_res_xis]

          calc_res_sum, calc_res_xis = GaPSE.sum_ξ_GNCxLD_multipole(COSMO.s_eff, s, COSMO; kwargs...)

          @test isapprox(res_sum, calc_res_sum; rtol = RTOL)
          @test all([isapprox(a, r; rtol = RTOL) for (a, r) in zip(calc_res_xis, res_xis)])
     end

     @testset "s = 500, L = 0, no_window" begin
          s = 500

          ind = findfirst(x -> x≈s, ss)
          res_sum =  all_res_sums[ind]
          res_xis = [vec[ind] for vec in all_res_xis]

          calc_res_sum, calc_res_xis = GaPSE.sum_ξ_GNCxLD_multipole(COSMO.s_eff, s, COSMO; kwargs...)

          @test isapprox(res_sum, calc_res_sum; rtol = RTOL)
          @test all([isapprox(a, r; rtol = RTOL) for (a, r) in zip(calc_res_xis, res_xis)])
     end


     @testset "s = 700, L = 0, no_window" begin
          s = 700

          ind = findfirst(x -> x≈s, ss)
          res_sum =  all_res_sums[ind]
          res_xis = [vec[ind] for vec in all_res_xis]

          calc_res_sum, calc_res_xis = GaPSE.sum_ξ_GNCxLD_multipole(COSMO.s_eff, s, COSMO; kwargs...)

          @test isapprox(res_sum, calc_res_sum; rtol = RTOL)
          @test all([isapprox(a, r; rtol = RTOL) for (a, r) in zip(calc_res_xis, res_xis)])
     end

     @testset "s = 1000, L = 0, no_window" begin
          s = 1000

          ind = findfirst(x -> x≈s, ss)
          res_sum =  all_res_sums[ind]
          res_xis = [vec[ind] for vec in all_res_xis]

          calc_res_sum, calc_res_xis = GaPSE.sum_ξ_GNCxLD_multipole(COSMO.s_eff, s, COSMO; kwargs...)

          @test isapprox(res_sum, calc_res_sum; rtol = RTOL)
          @test all([isapprox(a, r; rtol = RTOL) for (a, r) in zip(calc_res_xis, res_xis)])
     end

end

@testset "test sum_ξ_GNCxLD_multipole with_window" begin
     table = readdlm("datatest/GNCxLD_SumXiMultipoles/map_sum_xi_GNCxLD_L0_withF_specific_ss.txt", comments = true)
     ss = convert(Vector{Float64}, table[:, 1])
     all_res_sums = convert(Vector{Float64}, table[:, 2])
     all_res_xis = [convert(Vector{Float64}, table[:, i]) for i in 3:22]
     RTOL = 1e-2

     kwargs = Dict(
          :L => 0, :use_windows => true,
          :enhancer => 1e8, :N_μs => 30,
          :μ_atol => 0.0, :μ_rtol => 1e-2,
          :N_χs => 100, :N_χs_2 => 50,
          :pr => true,
     )

     #=
     @testset "zeros" begin
          a_func(x) = x
          @test_throws AssertionError GaPSE.sum_ξ_GNCxLD_multipole(COSMO.s_eff, 10, COSMO;
               L = 0, use_windows = false, SPLINE = true, kwargs...)
          @test_throws AssertionError GaPSE.integral_on_mu(COSMO.s_eff, 10, a_func, COSMO;
               L = 0, use_windows = false, SPLINE = true, kwargs...)
     end
     =#

     @testset "s = 10, L = 0, with_window" begin
          s = 10

          ind = findfirst(x -> x≈s, ss)
          res_sum =  all_res_sums[ind]
          res_xis = [vec[ind] for vec in all_res_xis] 

          calc_res_sum, calc_res_xis = GaPSE.sum_ξ_GNCxLD_multipole(COSMO.s_eff, s, COSMO; kwargs...)

          @test isapprox(res_sum, calc_res_sum; rtol = RTOL)
          @test all([isapprox(a, r; rtol = RTOL) for (a, r) in zip(calc_res_xis, res_xis)])
     end

     @testset "s = 100, L = 0, with_window" begin
          s = 100

          ind = findfirst(x -> x≈s, ss)
          res_sum =  all_res_sums[ind]
          res_xis = [vec[ind] for vec in all_res_xis]

          calc_res_sum, calc_res_xis = GaPSE.sum_ξ_GNCxLD_multipole(COSMO.s_eff, s, COSMO; kwargs...)

          @test isapprox(res_sum, calc_res_sum; rtol = RTOL)
          @test all([isapprox(a, r; rtol = RTOL) for (a, r) in zip(calc_res_xis, res_xis)])
     end

     @testset "s = 500, L = 0, with_window" begin
          s = 500

          ind = findfirst(x -> x≈s, ss)
          res_sum =  all_res_sums[ind]
          res_xis = [vec[ind] for vec in all_res_xis]

          calc_res_sum, calc_res_xis = GaPSE.sum_ξ_GNCxLD_multipole(COSMO.s_eff, s, COSMO; kwargs...)

          @test isapprox(res_sum, calc_res_sum; rtol = RTOL)
          @test all([isapprox(a, r; rtol = RTOL) for (a, r) in zip(calc_res_xis, res_xis)])
     end


     @testset "s = 700, L = 0, with_window" begin
          s = 700

          ind = findfirst(x -> x≈s, ss)
          res_sum =  all_res_sums[ind]
          res_xis = [vec[ind] for vec in all_res_xis]

          calc_res_sum, calc_res_xis = GaPSE.sum_ξ_GNCxLD_multipole(COSMO.s_eff, s, COSMO; kwargs...)

          @test isapprox(res_sum, calc_res_sum; rtol = RTOL)
          @test all([isapprox(a, r; rtol = RTOL) for (a, r) in zip(calc_res_xis, res_xis)])
     end

     @testset "s = 1000, L = 0, with_window" begin
          s = 1000

          ind = findfirst(x -> x≈s, ss)
          res_sum =  all_res_sums[ind]
          res_xis = [vec[ind] for vec in all_res_xis]

          calc_res_sum, calc_res_xis = GaPSE.sum_ξ_GNCxLD_multipole(COSMO.s_eff, s, COSMO; kwargs...)

          @test isapprox(res_sum, calc_res_sum; rtol = RTOL)
          @test all([isapprox(a, r; rtol = RTOL) for (a, r) in zip(calc_res_xis, res_xis)])
     end

end


@testset "test map_sum_ξ_GNCxLD_multipole no_window" begin
     table = readdlm("datatest/GNCxLD_SumXiMultipoles/map_sum_xi_GNCxLD_L0_noF_specific_ss.txt", comments = true)
     ss = convert(Vector{Float64}, table[:, 1])
     res_sums = convert(Vector{Float64}, table[:, 2])
     res_xis = [convert(Vector{Float64}, table[:, i]) for i in 3:18]

     RTOL = 1e-2

     kwargs = Dict(
          :s1 => nothing,
          :N_log => 3, :use_windows => false,
          :enhancer => 1e8, :N_μs => 30,
          :μ_atol => 0.0, :μ_rtol => 1e-2,
          :pr => true,
     )


     calc_ss, calc_sums, calc_xis = GaPSE.map_sum_ξ_GNCxLD_multipole(COSMO, ss; kwargs...)

     @test all([isapprox(a, r, rtol = RTOL) for (a, r) in zip(ss, calc_ss)])
     @test all([isapprox(a, r, rtol = RTOL) for (a, r) in zip(res_sums, calc_sums)])
     @test all([isapprox(a, r, rtol = RTOL) for (a, r) in zip(res_xis[1], calc_xis[1])]) # auto_doppler
     @test all([isapprox(a, r, rtol = RTOL) for (a, r) in zip(res_xis[2], calc_xis[2])]) # auto_lensing
     @test all([isapprox(a, r, rtol = RTOL) for (a, r) in zip(res_xis[3], calc_xis[3])]) # auto_localgp
     @test all([isapprox(a, r, rtol = RTOL) for (a, r) in zip(res_xis[4], calc_xis[4])]) # auto_integrated
     @test all([isapprox(a, r, rtol = RTOL) for (a, r) in zip(res_xis[5], calc_xis[5])]) # lensing_doppler
     @test all([isapprox(a, r, rtol = RTOL) for (a, r) in zip(res_xis[6], calc_xis[6])]) # doppler_lensing
     @test all([isapprox(a, r, rtol = RTOL) for (a, r) in zip(res_xis[7], calc_xis[7])]) # doppler_localgp
     @test all([isapprox(a, r, rtol = RTOL) for (a, r) in zip(res_xis[8], calc_xis[8])]) # localgp_doppler
     @test all([isapprox(a, r, rtol = RTOL) for (a, r) in zip(res_xis[9], calc_xis[9])]) # doppler_integratedgp
     @test all([isapprox(a, r, rtol = RTOL) for (a, r) in zip(res_xis[10], calc_xis[10])]) # integratedgp_doppler
     @test all([isapprox(a, r, rtol = RTOL) for (a, r) in zip(res_xis[11], calc_xis[11])]) # lensing_localgp
     @test all([isapprox(a, r, rtol = RTOL) for (a, r) in zip(res_xis[12], calc_xis[12])]) # localgp_lensing
     @test all([isapprox(a, r, rtol = RTOL) for (a, r) in zip(res_xis[13], calc_xis[13])]) # lensing_integratedgp
     @test all([isapprox(a, r, rtol = RTOL) for (a, r) in zip(res_xis[14], calc_xis[14])]) # integratedgp_lensing
     @test all([isapprox(a, r, rtol = RTOL) for (a, r) in zip(res_xis[15], calc_xis[15])]) # localgp_integatedgp
     @test all([isapprox(a, r, rtol = RTOL) for (a, r) in zip(res_xis[16], calc_xis[16])]) # integratedgp_localgp
end



@testset "test print_map_sum_ξ_GNCxLD_multipole no_window" begin
     table = readdlm("datatest/GNCxLD_SumXiMultipoles/map_sum_xi_GNCxLD_L0_noF_specific_ss.txt", comments = true)
     ss = convert(Vector{Float64}, table[:, 1])
     res_sums = convert(Vector{Float64}, table[:, 2])
     res_xis = [convert(Vector{Float64}, table[:, i]) for i in 3:18]

     RTOL = 1e-2
     kwargs = Dict(
          :s1 => nothing,
          :N_log => 3, :use_windows => false,
          :enhancer => 1e8, :N_μs => 30,
          :μ_atol => 0.0, :μ_rtol => 1e-2,
          :pr => true,
     )

     @testset "first" begin
          name = "first_map_sum_xi_GNCxLD_L0_noF_specific_ss.txt"
          isfile(name) && rm(name)

          GaPSE.print_map_sum_ξ_GNCxLD_multipole(COSMO, name, ss;
               single = true, kwargs...)

          calc_table = readdlm(name, comments = true)
          calc_ss = convert(Vector{Float64}, calc_table[:, 1])
          calc_sums = convert(Vector{Float64}, calc_table[:, 2])
          calc_xis = [convert(Vector{Float64}, calc_table[:, i]) for i in 3:18]

          @test all([isapprox(a, r, rtol = RTOL) for (a, r) in zip(ss, calc_ss)])
          @test all([isapprox(a, r, rtol = RTOL) for (a, r) in zip(res_sums, calc_sums)])
          @test all([isapprox(a, r, rtol = RTOL) for (a, r) in zip(res_xis[1], calc_xis[1])]) # auto_doppler
          @test all([isapprox(a, r, rtol = RTOL) for (a, r) in zip(res_xis[2], calc_xis[2])]) # auto_lensing
          @test all([isapprox(a, r, rtol = RTOL) for (a, r) in zip(res_xis[3], calc_xis[3])]) # auto_localgp
          @test all([isapprox(a, r, rtol = RTOL) for (a, r) in zip(res_xis[4], calc_xis[4])]) # auto_integrated
          @test all([isapprox(a, r, rtol = RTOL) for (a, r) in zip(res_xis[5], calc_xis[5])]) # lensing_doppler
          @test all([isapprox(a, r, rtol = RTOL) for (a, r) in zip(res_xis[6], calc_xis[6])]) # doppler_lensing
          @test all([isapprox(a, r, rtol = RTOL) for (a, r) in zip(res_xis[7], calc_xis[7])]) # doppler_localgp
          @test all([isapprox(a, r, rtol = RTOL) for (a, r) in zip(res_xis[8], calc_xis[8])]) # localgp_doppler
          @test all([isapprox(a, r, rtol = RTOL) for (a, r) in zip(res_xis[9], calc_xis[9])]) # doppler_integratedgp
          @test all([isapprox(a, r, rtol = RTOL) for (a, r) in zip(res_xis[10], calc_xis[10])]) # integratedgp_doppler
          @test all([isapprox(a, r, rtol = RTOL) for (a, r) in zip(res_xis[11], calc_xis[11])]) # lensing_localgp
          @test all([isapprox(a, r, rtol = RTOL) for (a, r) in zip(res_xis[12], calc_xis[12])]) # localgp_lensing
          @test all([isapprox(a, r, rtol = RTOL) for (a, r) in zip(res_xis[13], calc_xis[13])]) # lensing_integratedgp
          @test all([isapprox(a, r, rtol = RTOL) for (a, r) in zip(res_xis[14], calc_xis[14])]) # integratedgp_lensing
          @test all([isapprox(a, r, rtol = RTOL) for (a, r) in zip(res_xis[15], calc_xis[15])]) # localgp_integatedgp
          @test all([isapprox(a, r, rtol = RTOL) for (a, r) in zip(res_xis[16], calc_xis[16])]) # integratedgp_localgp

          rm(name)
     end

     @testset "second" begin
          name = "second_map_sum_xi_GNCxLD_L0_noF_specific_ss.txt"
          isfile(name) && rm(name)

          GaPSE.print_map_sum_ξ_GNCxLD_multipole(COSMO, name, ss;
               single = false, kwargs...)

          calc_table_1 = readdlm(name, comments = true)
          calc_ss = convert(Vector{Float64}, calc_table_1[:, 1])
          calc_sums = convert(Vector{Float64}, calc_table_1[:, 2])


          calc_xis = [
               begin
                    a_name = "all_GNCxLD_standalones_CFs/xi_GNCxLD_" * effect * "_L0" * ".txt"
                    calc_table = readdlm(a_name, comments = true)
                    #calc_ss = convert(Vector{Float64}, calc_table[:, 1])
                    xis = convert(Vector{Float64}, calc_table[:, 2])
                    xis
               end for effect in GaPSE.GR_EFFECTS_GNCxLD
          ]

          @test all([isapprox(a, r, rtol = RTOL) for (a, r) in zip(ss, calc_ss)])
          @test all([isapprox(a, r, rtol = RTOL) for (a, r) in zip(res_sums, calc_sums)])
          @test all([isapprox(a, r, rtol = RTOL) for (a, r) in zip(res_xis[1], calc_xis[1])]) # auto_doppler
          @test all([isapprox(a, r, rtol = RTOL) for (a, r) in zip(res_xis[2], calc_xis[2])]) # auto_lensing
          @test all([isapprox(a, r, rtol = RTOL) for (a, r) in zip(res_xis[3], calc_xis[3])]) # auto_localgp
          @test all([isapprox(a, r, rtol = RTOL) for (a, r) in zip(res_xis[4], calc_xis[4])]) # auto_integrated
          @test all([isapprox(a, r, rtol = RTOL) for (a, r) in zip(res_xis[5], calc_xis[5])]) # lensing_doppler
          @test all([isapprox(a, r, rtol = RTOL) for (a, r) in zip(res_xis[6], calc_xis[6])]) # doppler_lensing
          @test all([isapprox(a, r, rtol = RTOL) for (a, r) in zip(res_xis[7], calc_xis[7])]) # doppler_localgp
          @test all([isapprox(a, r, rtol = RTOL) for (a, r) in zip(res_xis[8], calc_xis[8])]) # localgp_doppler
          @test all([isapprox(a, r, rtol = RTOL) for (a, r) in zip(res_xis[9], calc_xis[9])]) # doppler_integratedgp
          @test all([isapprox(a, r, rtol = RTOL) for (a, r) in zip(res_xis[10], calc_xis[10])]) # integratedgp_doppler
          @test all([isapprox(a, r, rtol = RTOL) for (a, r) in zip(res_xis[11], calc_xis[11])]) # lensing_localgp
          @test all([isapprox(a, r, rtol = RTOL) for (a, r) in zip(res_xis[12], calc_xis[12])]) # localgp_lensing
          @test all([isapprox(a, r, rtol = RTOL) for (a, r) in zip(res_xis[13], calc_xis[13])]) # lensing_integratedgp
          @test all([isapprox(a, r, rtol = RTOL) for (a, r) in zip(res_xis[14], calc_xis[14])]) # integratedgp_lensing
          @test all([isapprox(a, r, rtol = RTOL) for (a, r) in zip(res_xis[15], calc_xis[15])]) # localgp_integatedgp
          @test all([isapprox(a, r, rtol = RTOL) for (a, r) in zip(res_xis[16], calc_xis[16])]) # integratedgp_localgp

          rm(name)
          for effect in GaPSE.GR_EFFECTS_GNCxLD
               rm("all_GNCxLD_standalones_CFs/xi_GNCxLD_" * effect * "_L0" * ".txt")
          end
          rm("all_GNCxLD_standalones_CFs")
     end
end

