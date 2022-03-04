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

@testset "test sum_ξ_multipole" begin
     RTOL = 1e-2

     kwargs = Dict(
          :L => 0, :use_windows => false,
          :enhancer => 1e8, :N_μs => 30,
          :μ_atol => 0.0, :μ_rtol => 1e-2,
          #:pr => false,
     )

     #=
     @testset "zeros" begin
          a_func(x) = x
          @test_throws AssertionError GaPSE.sum_ξ_multipole(COSMO.s_eff, 10, COSMO;
               L = 0, use_windows = false, SPLINE = true, kwargs...)
          @test_throws AssertionError GaPSE.integral_on_mu(COSMO.s_eff, 10, a_func, COSMO;
               L = 0, use_windows = false, SPLINE = true, kwargs...)
     end
     =#

     @testset "s = 10, L = 0, no_window" begin
          s, L, use_windows = 10, 0, false
          res_sum = 3.100057865775196e-5
          res_xis = [3.1083647323130395e-5, 1.617540202369208e-7, -1.21383701230076e-10,
               -3.49295089623204e-10, -1.6900809977615647e-7, -1.6924376436917377e-7,
               8.264905880080302e-10, 6.044013181911759e-10, 3.906461398291596e-8,
               3.906284620812885e-8, 7.272492754980249e-9, 7.27786631713832e-9,
               2.1861059508231394e-9, 2.187002470696239e-9, -2.2909254896728633e-9,
               -2.2910367803886207e-9]

          calc_res_sum, calc_res_xis = GaPSE.sum_ξ_multipole(COSMO.s_eff, s, COSMO;
               all = true, kwargs...)

          @test isapprox(res_sum, calc_res_sum; rtol = RTOL)
          @test all([isapprox(a, r; rtol = RTOL) for (a, r) in zip(calc_res_xis, res_xis)])
     end

     @testset "s = 100, L = 0, no_window" begin
          s, L, use_windows = 100, 0, false
          res_sum = 6.006178411066864e-6
          res_xis = [6.024865181025425e-6, 1.9780419647453728e-8, -4.9633712604631205e-9,
               -3.870682784239748e-10, -7.051019541134061e-8, -7.711583800345777e-8,
               2.292720951467199e-8, 1.698413720953249e-8, 3.2354497106195044e-8,
               3.224582344522693e-8, 5.74351175585059e-9, 6.104777136983635e-9,
               1.537461723506664e-9, 1.595320758808456e-9, -2.4870626291800128e-9,
               -2.496392673925445e-9]

          calc_res_sum, calc_res_xis = GaPSE.sum_ξ_multipole(COSMO.s_eff, s, COSMO;
               all = true, kwargs...)

          @test isapprox(res_sum, calc_res_sum; rtol = RTOL)
          @test all([isapprox(a, r; rtol = RTOL) for (a, r) in zip(calc_res_xis, res_xis)])
     end

     @testset "s = 700, L = 0, no_window" begin
          s, L, use_windows = 700, 0, false
          res_sum = 1.9540157917742675e-7
          res_xis = [1.4055178762472735e-7, 1.65683447635497e-9, -1.387147283648418e-8,
               -9.952329450745797e-10, 5.971511328265444e-10, -5.111667985927547e-9,
               1.3513756515145347e-8, 4.035832620379043e-8, 1.3147582840969097e-8,
               1.0360300419959947e-8, -5.233533247831251e-11, 2.498063130210084e-9,
               2.358619426410098e-10, 4.856778020909242e-10, -5.07087713847432e-9,
               -2.902176672849967e-9]

          calc_res_sum, calc_res_xis = GaPSE.sum_ξ_multipole(COSMO.s_eff, s, COSMO;
               all = true, kwargs...)

          @test isapprox(res_sum, calc_res_sum; rtol = RTOL)
          @test all([isapprox(a, r; rtol = RTOL) for (a, r) in zip(calc_res_xis, res_xis)])
     end
end

@testset "test map_sum_ξ_multipole" begin
     table = readdlm("datatest/map_sum_xi_L0.txt", comments = true)
     ss = convert(Vector{Float64}, table[:, 1])
     res_sums = convert(Vector{Float64}, table[:, 2])
     res_xis = [convert(Vector{Float64}, table[:, i]) for i in 3:18]

     RTOL = 1e-2

     kwargs = Dict(
          :s_1 => nothing,
          :N_log => 3, :use_windows => false,
          :enhancer => 1e8, :N_μs => 30,
          :μ_atol => 0.0, :μ_rtol => 1e-2,
          :pr => false,
     )

     #=
     res_sums = [3.100057865775196e-5, 6.006178411066864e-6, 9.834804209819796e-7, 1.9540157917742675e-7, 6.161330318538807e-8]
     res_xis = [
          [3.1083647323130395e-5, 6.024865181025425e-6, 8.992670696422657e-7, 1.4055178762472735e-7, 3.70429288148483e-8],
          [1.617540202369208e-7, 1.9780419647453728e-8, 4.390158923408416e-9, 1.65683447635497e-9, 1.148203496583022e-9],
          [-1.21383701230076e-10, -4.9633712604631205e-9, -1.526949939249473e-8, -1.387147283648418e-8, -8.701923426808867e-9],
          [-3.49295089623204e-10, -3.870682784239748e-10, -5.650486208049296e-10, -9.952329450745797e-10, -1.3189966155445926e-9],
          [-1.6900809977615647e-7, -7.051019541134061e-8, -7.554056283454021e-9, 5.971511328265444e-10, 2.437012122147995e-10],
          [-1.6924376436917377e-7, -7.711583800345777e-8, -2.230891563010909e-8, -5.111667985927547e-9, -2.1595251409863386e-9],
          [8.264905880080302e-10, 2.292720951467199e-8, 4.447934241479703e-8, 1.3513756515145347e-8, 5.172174705786513e-9],
          [6.044013181911759e-10, 1.698413720953249e-8, 3.7317043925737376e-8, 4.035832620379043e-8, 2.048989062874215e-8],
          [3.906461398291596e-8, 3.2354497106195044e-8, 2.1126957902449834e-8, 1.3147582840969097e-8, 1.0752716768250802e-8],
          [3.906284620812885e-8, 3.224582344522693e-8, 2.1113063300348826e-8, 1.0360300419959947e-8, 4.8564103242756225e-9],
          [7.272492754980249e-9, 5.74351175585059e-9, 2.3652633246115675e-9, -5.233533247831251e-11, -5.031090591479135e-11],
          [7.27786631713832e-9, 6.104777136983635e-9, 4.094282218195058e-9, 2.498063130210084e-9, 1.9502602679213973e-9],
          [2.1861059508231394e-9, 1.537461723506664e-9, 6.949249973729436e-10, 2.358619426410098e-10, 1.4896878162932332e-10],
          [2.187002470696239e-9, 1.595320758808456e-9, 9.061931506712828e-10, 4.856778020909242e-10, 3.647613731455436e-10],
          [-2.2909254896728633e-9, -2.4870626291800128e-9, -3.2701624458780454e-9, -5.07087713847432e-9, -6.476931202798078e-9],
          [-2.2910367803886207e-9, -2.496392673925445e-9, -3.3061964451371894e-9, -2.902176672849967e-9, -1.8490258959567266e-9]
     ]
     =#

     calc_ss, calc_sums, calc_xis = GaPSE.map_sum_ξ_multipole(COSMO, ss;
          all = true, kwargs...)

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



@testset "test print_map_sum_ξ_multipole" begin
     table = readdlm("datatest/map_sum_xi_L0.txt", comments = true)
     ss = convert(Vector{Float64}, table[:, 1])
     res_sums = convert(Vector{Float64}, table[:, 2])
     res_xis = [convert(Vector{Float64}, table[:, i]) for i in 3:18]

     RTOL = 1e-2
     kwargs = Dict(
          :s_1 => nothing,
          :N_log => 3, :use_windows => false,
          :enhancer => 1e8, :N_μs => 30,
          :μ_atol => 0.0, :μ_rtol => 1e-2,
          :pr => false,
     )

     @testset "first" begin
          name = "first_map_sum_xi_L0.txt"
          isfile(name) && rm(name)

          GaPSE.print_map_sum_ξ_multipole(COSMO, name, ss;
               all = true, single = true, kwargs...)

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
          name = "second_map_sum_xi_L0.txt"
          isfile(name) && rm(name)

          GaPSE.print_map_sum_ξ_multipole(COSMO, name, ss;
               all = true, single = false, kwargs...)

          calc_table_1 = readdlm(name, comments = true)
          calc_ss = convert(Vector{Float64}, calc_table_1[:, 1])
          calc_sums = convert(Vector{Float64}, calc_table_1[:, 2])


          calc_xis = [
               begin
                    a_name = "all_standalones_CFs/xi_" * effect * "_L0" * ".txt"
                    calc_table = readdlm(a_name, comments = true)
                    #calc_ss = convert(Vector{Float64}, calc_table[:, 1])
                    xis = convert(Vector{Float64}, calc_table[:, 2])
                    xis
               end for effect in GaPSE.IMPLEMENTED_GR_EFFECTS
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
          for effect in GaPSE.IMPLEMENTED_GR_EFFECTS
               rm("all_standalones_CFs/xi_" * effect * "_L0" * ".txt")
          end
          rm("all_standalones_CFs")
     end
end

