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

println("\n\nNow I start the tests of sum_ξ_LDxGNC_multipole and map_sum_ξ_LDxGNC_multipole")
println("It will take a while, but do not worry: I'm working.")

@testset "test sum_ξ_LDxGNC_multipole no_window" begin
     RTOL = 1.2e-2
     kwargs = Dict(
          :use_windows => false,
          :enhancer => 1e8,
          :N_trap => 30, :N_lob => 30,
          :atol_quad => 0.0, :rtol_quad => 1e-2,
          :N_χs => 40, :N_χs_2 => 20,
          #:pr => false,
     )

     #=
     @testset "zeros" begin
          a_func(x) = x
          @test_throws AssertionError GaPSE.sum_ξ_LDxGNC_multipole(COSMO.s_eff, 10, COSMO;
               L = 0, use_windows = false, SPLINE = true, kwargs...)
          @test_throws AssertionError GaPSE.integral_on_mu(COSMO.s_eff, 10, a_func, COSMO;
               L = 0, use_windows = false, SPLINE = true, kwargs...)
     end
     =#

     @testset "L = 0 no_window" begin
          L = 0
          ss_quad, res_sums_quad, res_xis_quad = GaPSE.readxyall(
               "datatest/LDxGNC_SumXiMultipoles/xis_LDxGNC_L$L" * "_quad_noF_specific_ss.txt", comments=true)
          ss_trap, res_sums_trap, res_xis_trap = GaPSE.readxyall(
               "datatest/LDxGNC_SumXiMultipoles/xis_LDxGNC_L$L" * "_trap_noF_specific_ss.txt", comments=true)
          ss_lob, res_sums_lob, res_xis_lob = GaPSE.readxyall(
               "datatest/LDxGNC_SumXiMultipoles/xis_LDxGNC_L$L" * "_lob_noF_specific_ss.txt", comments=true)


          @testset "s = 10, L = 0, no_window" begin
               s = 10

               ### quad ###
               ind_quad = findfirst(x -> x ≈ s, ss_quad)
               res_sum_spec_ss_quad = res_sums_quad[ind_quad]
               res_xis_spec_ss_quad = [vec[ind_quad] for vec in res_xis_quad]

               calc_res_sum_spec_ss_quad, calc_res_xis_spec_ss_quad =
                    GaPSE.sum_ξ_LDxGNC_multipole(COSMO.s_eff, s, COSMO; L=L, alg=:quad, kwargs...)

               @test isapprox(res_sum_spec_ss_quad, calc_res_sum_spec_ss_quad; rtol=RTOL)
               @test all([isapprox(a, r; rtol=RTOL) for (a, r) in zip(calc_res_xis_spec_ss_quad, res_xis_spec_ss_quad)])

               ### trap ###
               ind_trap = findfirst(x -> x ≈ s, ss_trap)
               res_sum_spec_ss_trap = res_sums_trap[ind_trap]
               res_xis_spec_ss_trap = [vec[ind_trap] for vec in res_xis_trap]

               calc_res_sum_spec_ss_trap, calc_res_xis_spec_ss_trap =
                    GaPSE.sum_ξ_LDxGNC_multipole(COSMO.s_eff, s, COSMO; L=L, alg=:trap, kwargs...)

               @test isapprox(res_sum_spec_ss_trap, calc_res_sum_spec_ss_trap; rtol=RTOL)
               @test all([isapprox(a, r; rtol=RTOL) for (a, r) in zip(calc_res_xis_spec_ss_trap, res_xis_spec_ss_trap)])

               #println("calc_res_xis_spec_ss_trap = $calc_res_xis_spec_ss_trap ;")
               #println("res_xis_spec_ss_trap = $res_xis_spec_ss_trap ;")


               ### lob ###
               ind_lob = findfirst(x -> x ≈ s, ss_lob)
               res_sum_spec_ss_lob = res_sums_lob[ind_lob]
               res_xis_spec_ss_lob = [vec[ind_lob] for vec in res_xis_lob]

               calc_res_sum_spec_ss_lob, calc_res_xis_spec_ss_lob =
                    GaPSE.sum_ξ_LDxGNC_multipole(COSMO.s_eff, s, COSMO; L=L, alg=:lobatto, kwargs...)

               @test isapprox(res_sum_spec_ss_lob, calc_res_sum_spec_ss_lob; rtol=RTOL)
               @test all([isapprox(a, r; rtol=RTOL) for (a, r) in zip(calc_res_xis_spec_ss_lob, res_xis_spec_ss_lob)])

          end

          @testset "s = 500, L = 0, no_window" begin
               s = 500

               ### quad ###
               ind_quad = findfirst(x -> x ≈ s, ss_quad)
               res_sum_spec_ss_quad = res_sums_quad[ind_quad]
               res_xis_spec_ss_quad = [vec[ind_quad] for vec in res_xis_quad]

               calc_res_sum_spec_ss_quad, calc_res_xis_spec_ss_quad =
                    GaPSE.sum_ξ_LDxGNC_multipole(COSMO.s_eff, s, COSMO; L=L, alg=:quad, kwargs...)

               @test isapprox(res_sum_spec_ss_quad, calc_res_sum_spec_ss_quad; rtol=RTOL)
               @test all([isapprox(a, r; rtol=RTOL) for (a, r) in zip(calc_res_xis_spec_ss_quad, res_xis_spec_ss_quad)])

               ### trap ###
               ind_trap = findfirst(x -> x ≈ s, ss_trap)
               res_sum_spec_ss_trap = res_sums_trap[ind_trap]
               res_xis_spec_ss_trap = [vec[ind_trap] for vec in res_xis_trap]

               calc_res_sum_spec_ss_trap, calc_res_xis_spec_ss_trap =
                    GaPSE.sum_ξ_LDxGNC_multipole(COSMO.s_eff, s, COSMO; L=L, alg=:trap, kwargs...)

               @test isapprox(res_sum_spec_ss_trap, calc_res_sum_spec_ss_trap; rtol=RTOL)
               @test all([isapprox(a, r; rtol=RTOL) for (a, r) in zip(calc_res_xis_spec_ss_trap, res_xis_spec_ss_trap)])

               #println("calc_res_xis_spec_ss_trap = $calc_res_xis_spec_ss_trap ;")
               #println("res_xis_spec_ss_trap = $res_xis_spec_ss_trap ;")


               ### lob ###
               ind_lob = findfirst(x -> x ≈ s, ss_lob)
               res_sum_spec_ss_lob = res_sums_lob[ind_lob]
               res_xis_spec_ss_lob = [vec[ind_lob] for vec in res_xis_lob]

               calc_res_sum_spec_ss_lob, calc_res_xis_spec_ss_lob =
                    GaPSE.sum_ξ_LDxGNC_multipole(COSMO.s_eff, s, COSMO; L=L, alg=:lobatto, kwargs...)

               @test isapprox(res_sum_spec_ss_lob, calc_res_sum_spec_ss_lob; rtol=RTOL)
               @test all([isapprox(a, r; rtol=RTOL) for (a, r) in zip(calc_res_xis_spec_ss_lob, res_xis_spec_ss_lob)])

          end

          @testset "s = 1000, L = 0, no_window" begin
               s = 1000

               ### quad ###
               ind_quad = findfirst(x -> x ≈ s, ss_quad)
               res_sum_spec_ss_quad = res_sums_quad[ind_quad]
               res_xis_spec_ss_quad = [vec[ind_quad] for vec in res_xis_quad]

               calc_res_sum_spec_ss_quad, calc_res_xis_spec_ss_quad =
                    GaPSE.sum_ξ_LDxGNC_multipole(COSMO.s_eff, s, COSMO; L=L, alg=:quad, kwargs...)

               @test isapprox(res_sum_spec_ss_quad, calc_res_sum_spec_ss_quad; rtol=RTOL)
               @test all([isapprox(a, r; rtol=RTOL) for (a, r) in zip(calc_res_xis_spec_ss_quad, res_xis_spec_ss_quad)])

               ### trap ###
               ind_trap = findfirst(x -> x ≈ s, ss_trap)
               res_sum_spec_ss_trap = res_sums_trap[ind_trap]
               res_xis_spec_ss_trap = [vec[ind_trap] for vec in res_xis_trap]

               calc_res_sum_spec_ss_trap, calc_res_xis_spec_ss_trap =
                    GaPSE.sum_ξ_LDxGNC_multipole(COSMO.s_eff, s, COSMO; L=L, alg=:trap, kwargs...)

               @test isapprox(res_sum_spec_ss_trap, calc_res_sum_spec_ss_trap; rtol=RTOL)
               @test all([isapprox(a, r; rtol=RTOL) for (a, r) in zip(calc_res_xis_spec_ss_trap, res_xis_spec_ss_trap)])

               #println("calc_res_xis_spec_ss_trap = $calc_res_xis_spec_ss_trap ;")
               #println("res_xis_spec_ss_trap = $res_xis_spec_ss_trap ;")


               ### lob ###
               ind_lob = findfirst(x -> x ≈ s, ss_lob)
               res_sum_spec_ss_lob = res_sums_lob[ind_lob]
               res_xis_spec_ss_lob = [vec[ind_lob] for vec in res_xis_lob]

               calc_res_sum_spec_ss_lob, calc_res_xis_spec_ss_lob =
                    GaPSE.sum_ξ_LDxGNC_multipole(COSMO.s_eff, s, COSMO; L=L, alg=:lobatto, kwargs...)

               @test isapprox(res_sum_spec_ss_lob, calc_res_sum_spec_ss_lob; rtol=RTOL)
               @test all([isapprox(a, r; rtol=RTOL) for (a, r) in zip(calc_res_xis_spec_ss_lob, res_xis_spec_ss_lob)])

          end
     end

     @testset "L = 1 no_window" begin
          L = 1

          ss_trap, res_sums_trap, res_xis_trap = GaPSE.readxyall(
               "datatest/LDxGNC_SumXiMultipoles/xis_LDxGNC_L$L" * "_trap_noF_specific_ss.txt", comments=true)
          ss_lob, res_sums_lob, res_xis_lob = GaPSE.readxyall(
               "datatest/LDxGNC_SumXiMultipoles/xis_LDxGNC_L$L" * "_lob_noF_specific_ss.txt", comments=true)


          @testset "s = 10, L = 1, no_window" begin
               s = 10

               ### trap ###
               ind_trap = findfirst(x -> x ≈ s, ss_trap)
               res_sum_spec_ss_trap = res_sums_trap[ind_trap]
               res_xis_spec_ss_trap = [vec[ind_trap] for vec in res_xis_trap]

               calc_res_sum_spec_ss_trap, calc_res_xis_spec_ss_trap =
                    GaPSE.sum_ξ_LDxGNC_multipole(COSMO.s_eff, s, COSMO; L=L, alg=:trap, kwargs...)

               @test isapprox(res_sum_spec_ss_trap, calc_res_sum_spec_ss_trap; rtol=RTOL)
               @test all([isapprox(a, r; rtol=RTOL) for (a, r) in zip(calc_res_xis_spec_ss_trap, res_xis_spec_ss_trap)])

               #println("calc_res_xis_spec_ss_trap = $calc_res_xis_spec_ss_trap ;")
               #println("res_xis_spec_ss_trap = $res_xis_spec_ss_trap ;")


               ### lob ###
               ind_lob = findfirst(x -> x ≈ s, ss_lob)
               res_sum_spec_ss_lob = res_sums_lob[ind_lob]
               res_xis_spec_ss_lob = [vec[ind_lob] for vec in res_xis_lob]

               calc_res_sum_spec_ss_lob, calc_res_xis_spec_ss_lob =
                    GaPSE.sum_ξ_LDxGNC_multipole(COSMO.s_eff, s, COSMO; L=L, alg=:lobatto, kwargs...)

               @test isapprox(res_sum_spec_ss_lob, calc_res_sum_spec_ss_lob; rtol=RTOL)
               @test all([isapprox(a, r; rtol=RTOL) for (a, r) in zip(calc_res_xis_spec_ss_lob, res_xis_spec_ss_lob)])

          end

          @testset "s = 500, L = 1, no_window" begin
               s = 500


               ### trap ###
               ind_trap = findfirst(x -> x ≈ s, ss_trap)
               res_sum_spec_ss_trap = res_sums_trap[ind_trap]
               res_xis_spec_ss_trap = [vec[ind_trap] for vec in res_xis_trap]

               calc_res_sum_spec_ss_trap, calc_res_xis_spec_ss_trap =
                    GaPSE.sum_ξ_LDxGNC_multipole(COSMO.s_eff, s, COSMO; L=L, alg=:trap, kwargs...)

               @test isapprox(res_sum_spec_ss_trap, calc_res_sum_spec_ss_trap; rtol=RTOL)
               @test all([isapprox(a, r; rtol=RTOL) for (a, r) in zip(calc_res_xis_spec_ss_trap, res_xis_spec_ss_trap)])

               #println("calc_res_xis_spec_ss_trap = $calc_res_xis_spec_ss_trap ;")
               #println("res_xis_spec_ss_trap = $res_xis_spec_ss_trap ;")


               ### lob ###
               ind_lob = findfirst(x -> x ≈ s, ss_lob)
               res_sum_spec_ss_lob = res_sums_lob[ind_lob]
               res_xis_spec_ss_lob = [vec[ind_lob] for vec in res_xis_lob]

               calc_res_sum_spec_ss_lob, calc_res_xis_spec_ss_lob =
                    GaPSE.sum_ξ_LDxGNC_multipole(COSMO.s_eff, s, COSMO; L=L, alg=:lobatto, kwargs...)

               @test isapprox(res_sum_spec_ss_lob, calc_res_sum_spec_ss_lob; rtol=RTOL)
               @test all([isapprox(a, r; rtol=RTOL) for (a, r) in zip(calc_res_xis_spec_ss_lob, res_xis_spec_ss_lob)])

          end

          @testset "s = 1000, L = 1, no_window" begin
               s = 1000

               ### trap ###
               ind_trap = findfirst(x -> x ≈ s, ss_trap)
               res_sum_spec_ss_trap = res_sums_trap[ind_trap]
               res_xis_spec_ss_trap = [vec[ind_trap] for vec in res_xis_trap]

               calc_res_sum_spec_ss_trap, calc_res_xis_spec_ss_trap =
                    GaPSE.sum_ξ_LDxGNC_multipole(COSMO.s_eff, s, COSMO; L=L, alg=:trap, kwargs...)

               @test isapprox(res_sum_spec_ss_trap, calc_res_sum_spec_ss_trap; rtol=RTOL)
               @test all([isapprox(a, r; rtol=RTOL) for (a, r) in zip(calc_res_xis_spec_ss_trap, res_xis_spec_ss_trap)])

               #println("calc_res_xis_spec_ss_trap = $calc_res_xis_spec_ss_trap ;")
               #println("res_xis_spec_ss_trap = $res_xis_spec_ss_trap ;")


               ### lob ###
               ind_lob = findfirst(x -> x ≈ s, ss_lob)
               res_sum_spec_ss_lob = res_sums_lob[ind_lob]
               res_xis_spec_ss_lob = [vec[ind_lob] for vec in res_xis_lob]

               calc_res_sum_spec_ss_lob, calc_res_xis_spec_ss_lob =
                    GaPSE.sum_ξ_LDxGNC_multipole(COSMO.s_eff, s, COSMO; L=L, alg=:lobatto, kwargs...)

               @test isapprox(res_sum_spec_ss_lob, calc_res_sum_spec_ss_lob; rtol=RTOL)
               @test all([isapprox(a, r; rtol=RTOL) for (a, r) in zip(calc_res_xis_spec_ss_lob, res_xis_spec_ss_lob)])

          end
     end

     @testset "L = 2 no_window" begin
          L = 2
          ss_quad, res_sums_quad, res_xis_quad = GaPSE.readxyall(
               "datatest/LDxGNC_SumXiMultipoles/xis_LDxGNC_L$L" * "_quad_noF_specific_ss.txt", comments=true)
          ss_trap, res_sums_trap, res_xis_trap = GaPSE.readxyall(
               "datatest/LDxGNC_SumXiMultipoles/xis_LDxGNC_L$L" * "_trap_noF_specific_ss.txt", comments=true)
          ss_lob, res_sums_lob, res_xis_lob = GaPSE.readxyall(
               "datatest/LDxGNC_SumXiMultipoles/xis_LDxGNC_L$L" * "_lob_noF_specific_ss.txt", comments=true)


          @testset "s = 10, L = 2, no_window" begin
               s = 10

               ### quad ###
               ind_quad = findfirst(x -> x ≈ s, ss_quad)
               res_sum_spec_ss_quad = res_sums_quad[ind_quad]
               res_xis_spec_ss_quad = [vec[ind_quad] for vec in res_xis_quad]

               calc_res_sum_spec_ss_quad, calc_res_xis_spec_ss_quad =
                    GaPSE.sum_ξ_LDxGNC_multipole(COSMO.s_eff, s, COSMO; L=L, alg=:quad, kwargs...)

               @test isapprox(res_sum_spec_ss_quad, calc_res_sum_spec_ss_quad; rtol=RTOL)
               @test all([isapprox(a, r; rtol=RTOL) for (a, r) in zip(calc_res_xis_spec_ss_quad, res_xis_spec_ss_quad)])

               ### trap ###
               ind_trap = findfirst(x -> x ≈ s, ss_trap)
               res_sum_spec_ss_trap = res_sums_trap[ind_trap]
               res_xis_spec_ss_trap = [vec[ind_trap] for vec in res_xis_trap]

               calc_res_sum_spec_ss_trap, calc_res_xis_spec_ss_trap =
                    GaPSE.sum_ξ_LDxGNC_multipole(COSMO.s_eff, s, COSMO; L=L, alg=:trap, kwargs...)

               @test isapprox(res_sum_spec_ss_trap, calc_res_sum_spec_ss_trap; rtol=RTOL)
               @test all([isapprox(a, r; rtol=RTOL) for (a, r) in zip(calc_res_xis_spec_ss_trap, res_xis_spec_ss_trap)])

               #println("calc_res_xis_spec_ss_trap = $calc_res_xis_spec_ss_trap ;")
               #println("res_xis_spec_ss_trap = $res_xis_spec_ss_trap ;")


               ### lob ###
               ind_lob = findfirst(x -> x ≈ s, ss_lob)
               res_sum_spec_ss_lob = res_sums_lob[ind_lob]
               res_xis_spec_ss_lob = [vec[ind_lob] for vec in res_xis_lob]

               calc_res_sum_spec_ss_lob, calc_res_xis_spec_ss_lob =
                    GaPSE.sum_ξ_LDxGNC_multipole(COSMO.s_eff, s, COSMO; L=L, alg=:lobatto, kwargs...)

               @test isapprox(res_sum_spec_ss_lob, calc_res_sum_spec_ss_lob; rtol=RTOL)
               @test all([isapprox(a, r; rtol=RTOL) for (a, r) in zip(calc_res_xis_spec_ss_lob, res_xis_spec_ss_lob)])

          end

          @testset "s = 500, L = 2, no_window" begin
               s = 500

               ### quad ###
               ind_quad = findfirst(x -> x ≈ s, ss_quad)
               res_sum_spec_ss_quad = res_sums_quad[ind_quad]
               res_xis_spec_ss_quad = [vec[ind_quad] for vec in res_xis_quad]

               calc_res_sum_spec_ss_quad, calc_res_xis_spec_ss_quad =
                    GaPSE.sum_ξ_LDxGNC_multipole(COSMO.s_eff, s, COSMO; L=L, alg=:quad, kwargs...)

               @test isapprox(res_sum_spec_ss_quad, calc_res_sum_spec_ss_quad; rtol=RTOL)
               @test all([isapprox(a, r; rtol=RTOL) for (a, r) in zip(calc_res_xis_spec_ss_quad, res_xis_spec_ss_quad)])

               ### trap ###
               ind_trap = findfirst(x -> x ≈ s, ss_trap)
               res_sum_spec_ss_trap = res_sums_trap[ind_trap]
               res_xis_spec_ss_trap = [vec[ind_trap] for vec in res_xis_trap]

               calc_res_sum_spec_ss_trap, calc_res_xis_spec_ss_trap =
                    GaPSE.sum_ξ_LDxGNC_multipole(COSMO.s_eff, s, COSMO; L=L, alg=:trap, kwargs...)

               @test isapprox(res_sum_spec_ss_trap, calc_res_sum_spec_ss_trap; rtol=RTOL)
               @test all([isapprox(a, r; rtol=RTOL) for (a, r) in zip(calc_res_xis_spec_ss_trap, res_xis_spec_ss_trap)])

               #println("calc_res_xis_spec_ss_trap = $calc_res_xis_spec_ss_trap ;")
               #println("res_xis_spec_ss_trap = $res_xis_spec_ss_trap ;")


               ### lob ###
               ind_lob = findfirst(x -> x ≈ s, ss_lob)
               res_sum_spec_ss_lob = res_sums_lob[ind_lob]
               res_xis_spec_ss_lob = [vec[ind_lob] for vec in res_xis_lob]

               calc_res_sum_spec_ss_lob, calc_res_xis_spec_ss_lob =
                    GaPSE.sum_ξ_LDxGNC_multipole(COSMO.s_eff, s, COSMO; L=L, alg=:lobatto, kwargs...)

               @test isapprox(res_sum_spec_ss_lob, calc_res_sum_spec_ss_lob; rtol=RTOL)
               @test all([isapprox(a, r; rtol=RTOL) for (a, r) in zip(calc_res_xis_spec_ss_lob, res_xis_spec_ss_lob)])

          end

          @testset "s = 1000, L = 2, no_window" begin
               s = 1000

               ### quad ###
               ind_quad = findfirst(x -> x ≈ s, ss_quad)
               res_sum_spec_ss_quad = res_sums_quad[ind_quad]
               res_xis_spec_ss_quad = [vec[ind_quad] for vec in res_xis_quad]

               calc_res_sum_spec_ss_quad, calc_res_xis_spec_ss_quad =
                    GaPSE.sum_ξ_LDxGNC_multipole(COSMO.s_eff, s, COSMO; L=L, alg=:quad, kwargs...)

               @test isapprox(res_sum_spec_ss_quad, calc_res_sum_spec_ss_quad; rtol=RTOL)
               @test all([isapprox(a, r; rtol=RTOL) for (a, r) in zip(calc_res_xis_spec_ss_quad, res_xis_spec_ss_quad)])

               ### trap ###
               ind_trap = findfirst(x -> x ≈ s, ss_trap)
               res_sum_spec_ss_trap = res_sums_trap[ind_trap]
               res_xis_spec_ss_trap = [vec[ind_trap] for vec in res_xis_trap]

               calc_res_sum_spec_ss_trap, calc_res_xis_spec_ss_trap =
                    GaPSE.sum_ξ_LDxGNC_multipole(COSMO.s_eff, s, COSMO; L=L, alg=:trap, kwargs...)

               @test isapprox(res_sum_spec_ss_trap, calc_res_sum_spec_ss_trap; rtol=RTOL)
               @test all([isapprox(a, r; rtol=RTOL) for (a, r) in zip(calc_res_xis_spec_ss_trap, res_xis_spec_ss_trap)])

               #println("calc_res_xis_spec_ss_trap = $calc_res_xis_spec_ss_trap ;")
               #println("res_xis_spec_ss_trap = $res_xis_spec_ss_trap ;")


               ### lob ###
               ind_lob = findfirst(x -> x ≈ s, ss_lob)
               res_sum_spec_ss_lob = res_sums_lob[ind_lob]
               res_xis_spec_ss_lob = [vec[ind_lob] for vec in res_xis_lob]

               calc_res_sum_spec_ss_lob, calc_res_xis_spec_ss_lob =
                    GaPSE.sum_ξ_LDxGNC_multipole(COSMO.s_eff, s, COSMO; L=L, alg=:lobatto, kwargs...)

               @test isapprox(res_sum_spec_ss_lob, calc_res_sum_spec_ss_lob; rtol=RTOL)
               @test all([isapprox(a, r; rtol=RTOL) for (a, r) in zip(calc_res_xis_spec_ss_lob, res_xis_spec_ss_lob)])

          end
     end

     @testset "L = 3 no_window" begin
          L = 3

          ss_trap, res_sums_trap, res_xis_trap = GaPSE.readxyall(
               "datatest/LDxGNC_SumXiMultipoles/xis_LDxGNC_L$L" * "_trap_noF_specific_ss.txt", comments=true)
          ss_lob, res_sums_lob, res_xis_lob = GaPSE.readxyall(
               "datatest/LDxGNC_SumXiMultipoles/xis_LDxGNC_L$L" * "_lob_noF_specific_ss.txt", comments=true)


          @testset "s = 10, L = 1, no_window" begin
               s = 10

               ### trap ###
               ind_trap = findfirst(x -> x ≈ s, ss_trap)
               res_sum_spec_ss_trap = res_sums_trap[ind_trap]
               res_xis_spec_ss_trap = [vec[ind_trap] for vec in res_xis_trap]

               calc_res_sum_spec_ss_trap, calc_res_xis_spec_ss_trap =
                    GaPSE.sum_ξ_LDxGNC_multipole(COSMO.s_eff, s, COSMO; L=L, alg=:trap, kwargs...)

               @test isapprox(res_sum_spec_ss_trap, calc_res_sum_spec_ss_trap; rtol=RTOL)
               @test all([isapprox(a, r; rtol=RTOL) for (a, r) in zip(calc_res_xis_spec_ss_trap, res_xis_spec_ss_trap)])

               #println("calc_res_xis_spec_ss_trap = $calc_res_xis_spec_ss_trap ;")
               #println("res_xis_spec_ss_trap = $res_xis_spec_ss_trap ;")


               ### lob ###
               ind_lob = findfirst(x -> x ≈ s, ss_lob)
               res_sum_spec_ss_lob = res_sums_lob[ind_lob]
               res_xis_spec_ss_lob = [vec[ind_lob] for vec in res_xis_lob]

               calc_res_sum_spec_ss_lob, calc_res_xis_spec_ss_lob =
                    GaPSE.sum_ξ_LDxGNC_multipole(COSMO.s_eff, s, COSMO; L=L, alg=:lobatto, kwargs...)

               @test isapprox(res_sum_spec_ss_lob, calc_res_sum_spec_ss_lob; rtol=RTOL)
               @test all([isapprox(a, r; rtol=RTOL) for (a, r) in zip(calc_res_xis_spec_ss_lob, res_xis_spec_ss_lob)])

          end

          @testset "s = 500, L = 1, no_window" begin
               s = 500


               ### trap ###
               ind_trap = findfirst(x -> x ≈ s, ss_trap)
               res_sum_spec_ss_trap = res_sums_trap[ind_trap]
               res_xis_spec_ss_trap = [vec[ind_trap] for vec in res_xis_trap]

               calc_res_sum_spec_ss_trap, calc_res_xis_spec_ss_trap =
                    GaPSE.sum_ξ_LDxGNC_multipole(COSMO.s_eff, s, COSMO; L=L, alg=:trap, kwargs...)

               @test isapprox(res_sum_spec_ss_trap, calc_res_sum_spec_ss_trap; rtol=RTOL)
               @test all([isapprox(a, r; rtol=RTOL) for (a, r) in zip(calc_res_xis_spec_ss_trap, res_xis_spec_ss_trap)])

               #println("calc_res_xis_spec_ss_trap = $calc_res_xis_spec_ss_trap ;")
               #println("res_xis_spec_ss_trap = $res_xis_spec_ss_trap ;")


               ### lob ###
               ind_lob = findfirst(x -> x ≈ s, ss_lob)
               res_sum_spec_ss_lob = res_sums_lob[ind_lob]
               res_xis_spec_ss_lob = [vec[ind_lob] for vec in res_xis_lob]

               calc_res_sum_spec_ss_lob, calc_res_xis_spec_ss_lob =
                    GaPSE.sum_ξ_LDxGNC_multipole(COSMO.s_eff, s, COSMO; L=L, alg=:lobatto, kwargs...)

               @test isapprox(res_sum_spec_ss_lob, calc_res_sum_spec_ss_lob; rtol=RTOL)
               @test all([isapprox(a, r; rtol=RTOL) for (a, r) in zip(calc_res_xis_spec_ss_lob, res_xis_spec_ss_lob)])

          end

          @testset "s = 1000, L = 1, no_window" begin
               s = 1000

               ### trap ###
               ind_trap = findfirst(x -> x ≈ s, ss_trap)
               res_sum_spec_ss_trap = res_sums_trap[ind_trap]
               res_xis_spec_ss_trap = [vec[ind_trap] for vec in res_xis_trap]

               calc_res_sum_spec_ss_trap, calc_res_xis_spec_ss_trap =
                    GaPSE.sum_ξ_LDxGNC_multipole(COSMO.s_eff, s, COSMO; L=L, alg=:trap, kwargs...)

               @test isapprox(res_sum_spec_ss_trap, calc_res_sum_spec_ss_trap; rtol=RTOL)
               @test all([isapprox(a, r; rtol=RTOL) for (a, r) in zip(calc_res_xis_spec_ss_trap, res_xis_spec_ss_trap)])

               #println("calc_res_xis_spec_ss_trap = $calc_res_xis_spec_ss_trap ;")
               #println("res_xis_spec_ss_trap = $res_xis_spec_ss_trap ;")


               ### lob ###
               ind_lob = findfirst(x -> x ≈ s, ss_lob)
               res_sum_spec_ss_lob = res_sums_lob[ind_lob]
               res_xis_spec_ss_lob = [vec[ind_lob] for vec in res_xis_lob]

               calc_res_sum_spec_ss_lob, calc_res_xis_spec_ss_lob =
                    GaPSE.sum_ξ_LDxGNC_multipole(COSMO.s_eff, s, COSMO; L=L, alg=:lobatto, kwargs...)

               @test isapprox(res_sum_spec_ss_lob, calc_res_sum_spec_ss_lob; rtol=RTOL)
               @test all([isapprox(a, r; rtol=RTOL) for (a, r) in zip(calc_res_xis_spec_ss_lob, res_xis_spec_ss_lob)])

          end
     end

     @testset "L = 4 no_window" begin
          L = 4

          ss_trap, res_sums_trap, res_xis_trap = GaPSE.readxyall(
               "datatest/LDxGNC_SumXiMultipoles/xis_LDxGNC_L$L" * "_trap_noF_specific_ss.txt", comments=true)
          ss_lob, res_sums_lob, res_xis_lob = GaPSE.readxyall(
               "datatest/LDxGNC_SumXiMultipoles/xis_LDxGNC_L$L" * "_lob_noF_specific_ss.txt", comments=true)


          @testset "s = 10, L = 1, no_window" begin
               s = 10

               ### trap ###
               ind_trap = findfirst(x -> x ≈ s, ss_trap)
               res_sum_spec_ss_trap = res_sums_trap[ind_trap]
               res_xis_spec_ss_trap = [vec[ind_trap] for vec in res_xis_trap]

               calc_res_sum_spec_ss_trap, calc_res_xis_spec_ss_trap =
                    GaPSE.sum_ξ_LDxGNC_multipole(COSMO.s_eff, s, COSMO; L=L, alg=:trap, kwargs...)

               @test isapprox(res_sum_spec_ss_trap, calc_res_sum_spec_ss_trap; rtol=RTOL)
               @test all([isapprox(a, r; rtol=RTOL) for (a, r) in zip(calc_res_xis_spec_ss_trap, res_xis_spec_ss_trap)])

               #println("calc_res_xis_spec_ss_trap = $calc_res_xis_spec_ss_trap ;")
               #println("res_xis_spec_ss_trap = $res_xis_spec_ss_trap ;")


               ### lob ###
               ind_lob = findfirst(x -> x ≈ s, ss_lob)
               res_sum_spec_ss_lob = res_sums_lob[ind_lob]
               res_xis_spec_ss_lob = [vec[ind_lob] for vec in res_xis_lob]

               calc_res_sum_spec_ss_lob, calc_res_xis_spec_ss_lob =
                    GaPSE.sum_ξ_LDxGNC_multipole(COSMO.s_eff, s, COSMO; L=L, alg=:lobatto, kwargs...)

               @test isapprox(res_sum_spec_ss_lob, calc_res_sum_spec_ss_lob; rtol=RTOL)
               @test all([isapprox(a, r; rtol=RTOL) for (a, r) in zip(calc_res_xis_spec_ss_lob, res_xis_spec_ss_lob)])

          end

          @testset "s = 500, L = 1, no_window" begin
               s = 500


               ### trap ###
               ind_trap = findfirst(x -> x ≈ s, ss_trap)
               res_sum_spec_ss_trap = res_sums_trap[ind_trap]
               res_xis_spec_ss_trap = [vec[ind_trap] for vec in res_xis_trap]

               calc_res_sum_spec_ss_trap, calc_res_xis_spec_ss_trap =
                    GaPSE.sum_ξ_LDxGNC_multipole(COSMO.s_eff, s, COSMO; L=L, alg=:trap, kwargs...)

               @test isapprox(res_sum_spec_ss_trap, calc_res_sum_spec_ss_trap; rtol=RTOL)
               @test all([isapprox(a, r; rtol=RTOL) for (a, r) in zip(calc_res_xis_spec_ss_trap, res_xis_spec_ss_trap)])

               #println("calc_res_xis_spec_ss_trap = $calc_res_xis_spec_ss_trap ;")
               #println("res_xis_spec_ss_trap = $res_xis_spec_ss_trap ;")


               ### lob ###
               ind_lob = findfirst(x -> x ≈ s, ss_lob)
               res_sum_spec_ss_lob = res_sums_lob[ind_lob]
               res_xis_spec_ss_lob = [vec[ind_lob] for vec in res_xis_lob]

               calc_res_sum_spec_ss_lob, calc_res_xis_spec_ss_lob =
                    GaPSE.sum_ξ_LDxGNC_multipole(COSMO.s_eff, s, COSMO; L=L, alg=:lobatto, kwargs...)

               @test isapprox(res_sum_spec_ss_lob, calc_res_sum_spec_ss_lob; rtol=RTOL)
               @test all([isapprox(a, r; rtol=RTOL) for (a, r) in zip(calc_res_xis_spec_ss_lob, res_xis_spec_ss_lob)])

          end

          @testset "s = 1000, L = 1, no_window" begin
               s = 1000

               ### trap ###
               ind_trap = findfirst(x -> x ≈ s, ss_trap)
               res_sum_spec_ss_trap = res_sums_trap[ind_trap]
               res_xis_spec_ss_trap = [vec[ind_trap] for vec in res_xis_trap]

               calc_res_sum_spec_ss_trap, calc_res_xis_spec_ss_trap =
                    GaPSE.sum_ξ_LDxGNC_multipole(COSMO.s_eff, s, COSMO; L=L, alg=:trap, kwargs...)

               @test isapprox(res_sum_spec_ss_trap, calc_res_sum_spec_ss_trap; rtol=RTOL)
               @test all([isapprox(a, r; rtol=RTOL) for (a, r) in zip(calc_res_xis_spec_ss_trap, res_xis_spec_ss_trap)])

               #println("calc_res_xis_spec_ss_trap = $calc_res_xis_spec_ss_trap ;")
               #println("res_xis_spec_ss_trap = $res_xis_spec_ss_trap ;")


               ### lob ###
               ind_lob = findfirst(x -> x ≈ s, ss_lob)
               res_sum_spec_ss_lob = res_sums_lob[ind_lob]
               res_xis_spec_ss_lob = [vec[ind_lob] for vec in res_xis_lob]

               calc_res_sum_spec_ss_lob, calc_res_xis_spec_ss_lob =
                    GaPSE.sum_ξ_LDxGNC_multipole(COSMO.s_eff, s, COSMO; L=L, alg=:lobatto, kwargs...)

               @test isapprox(res_sum_spec_ss_lob, calc_res_sum_spec_ss_lob; rtol=RTOL)
               @test all([isapprox(a, r; rtol=RTOL) for (a, r) in zip(calc_res_xis_spec_ss_lob, res_xis_spec_ss_lob)])

          end
     end
end

println("\nDon't worry, I am on it...")

@testset "test sum_ξ_LDxGNC_multipole with_window" begin
     RTOL = 1.2e-2
     kwargs = Dict(
          :use_windows => true,
          :enhancer => 1e8,
          :N_trap => 30, :N_lob => 30,
          :atol_quad => 0.0, :rtol_quad => 1e-2,
          :N_χs => 40, :N_χs_2 => 20,
          #:pr => false,
     )

     #=
     @testset "zeros" begin
          a_func(x) = x
          @test_throws AssertionError GaPSE.sum_ξ_LDxGNC_multipole(COSMO.s_eff, 10, COSMO;
               L = 0, use_windows = false, SPLINE = true, kwargs...)
          @test_throws AssertionError GaPSE.integral_on_mu(COSMO.s_eff, 10, a_func, COSMO;
               L = 0, use_windows = false, SPLINE = true, kwargs...)
     end
     =#

     @testset "L = 0 with_window" begin
          L = 0
          ss_quad, res_sums_quad, res_xis_quad = GaPSE.readxyall(
               "datatest/LDxGNC_SumXiMultipoles/xis_LDxGNC_L$L" * "_quad_withF_specific_ss.txt", comments=true)
          ss_trap, res_sums_trap, res_xis_trap = GaPSE.readxyall(
               "datatest/LDxGNC_SumXiMultipoles/xis_LDxGNC_L$L" * "_trap_withF_specific_ss.txt", comments=true)
          ss_lob, res_sums_lob, res_xis_lob = GaPSE.readxyall(
               "datatest/LDxGNC_SumXiMultipoles/xis_LDxGNC_L$L" * "_lob_withF_specific_ss.txt", comments=true)


          @testset "s = 10, L = 0, with_window" begin
               s = 10

               ### quad ###
               ind_quad = findfirst(x -> x ≈ s, ss_quad)
               res_sum_spec_ss_quad = res_sums_quad[ind_quad]
               res_xis_spec_ss_quad = [vec[ind_quad] for vec in res_xis_quad]

               calc_res_sum_spec_ss_quad, calc_res_xis_spec_ss_quad =
                    GaPSE.sum_ξ_LDxGNC_multipole(COSMO.s_eff, s, COSMO; L=L, alg=:quad, kwargs...)

               @test isapprox(res_sum_spec_ss_quad, calc_res_sum_spec_ss_quad; rtol=RTOL)
               @test all([isapprox(a, r; rtol=RTOL) for (a, r) in zip(calc_res_xis_spec_ss_quad, res_xis_spec_ss_quad)])

               ### trap ###
               ind_trap = findfirst(x -> x ≈ s, ss_trap)
               res_sum_spec_ss_trap = res_sums_trap[ind_trap]
               res_xis_spec_ss_trap = [vec[ind_trap] for vec in res_xis_trap]

               calc_res_sum_spec_ss_trap, calc_res_xis_spec_ss_trap =
                    GaPSE.sum_ξ_LDxGNC_multipole(COSMO.s_eff, s, COSMO; L=L, alg=:trap, kwargs...)

               @test isapprox(res_sum_spec_ss_trap, calc_res_sum_spec_ss_trap; rtol=RTOL)
               @test all([isapprox(a, r; rtol=RTOL) for (a, r) in zip(calc_res_xis_spec_ss_trap, res_xis_spec_ss_trap)])

               #println("calc_res_xis_spec_ss_trap = $calc_res_xis_spec_ss_trap ;")
               #println("res_xis_spec_ss_trap = $res_xis_spec_ss_trap ;")


               ### lob ###
               ind_lob = findfirst(x -> x ≈ s, ss_lob)
               res_sum_spec_ss_lob = res_sums_lob[ind_lob]
               res_xis_spec_ss_lob = [vec[ind_lob] for vec in res_xis_lob]

               calc_res_sum_spec_ss_lob, calc_res_xis_spec_ss_lob =
                    GaPSE.sum_ξ_LDxGNC_multipole(COSMO.s_eff, s, COSMO; L=L, alg=:lobatto, kwargs...)

               @test isapprox(res_sum_spec_ss_lob, calc_res_sum_spec_ss_lob; rtol=RTOL)
               @test all([isapprox(a, r; rtol=RTOL) for (a, r) in zip(calc_res_xis_spec_ss_lob, res_xis_spec_ss_lob)])

          end

          @testset "s = 500, L = 0, with_window" begin
               s = 500

               ### quad ###
               ind_quad = findfirst(x -> x ≈ s, ss_quad)
               res_sum_spec_ss_quad = res_sums_quad[ind_quad]
               res_xis_spec_ss_quad = [vec[ind_quad] for vec in res_xis_quad]

               calc_res_sum_spec_ss_quad, calc_res_xis_spec_ss_quad =
                    GaPSE.sum_ξ_LDxGNC_multipole(COSMO.s_eff, s, COSMO; L=L, alg=:quad, kwargs...)

               @test isapprox(res_sum_spec_ss_quad, calc_res_sum_spec_ss_quad; rtol=RTOL)
               @test all([isapprox(a, r; rtol=RTOL) for (a, r) in zip(calc_res_xis_spec_ss_quad, res_xis_spec_ss_quad)])

               ### trap ###
               ind_trap = findfirst(x -> x ≈ s, ss_trap)
               res_sum_spec_ss_trap = res_sums_trap[ind_trap]
               res_xis_spec_ss_trap = [vec[ind_trap] for vec in res_xis_trap]

               calc_res_sum_spec_ss_trap, calc_res_xis_spec_ss_trap =
                    GaPSE.sum_ξ_LDxGNC_multipole(COSMO.s_eff, s, COSMO; L=L, alg=:trap, kwargs...)

               @test isapprox(res_sum_spec_ss_trap, calc_res_sum_spec_ss_trap; rtol=RTOL)
               @test all([isapprox(a, r; rtol=RTOL) for (a, r) in zip(calc_res_xis_spec_ss_trap, res_xis_spec_ss_trap)])

               #println("calc_res_xis_spec_ss_trap = $calc_res_xis_spec_ss_trap ;")
               #println("res_xis_spec_ss_trap = $res_xis_spec_ss_trap ;")


               ### lob ###
               ind_lob = findfirst(x -> x ≈ s, ss_lob)
               res_sum_spec_ss_lob = res_sums_lob[ind_lob]
               res_xis_spec_ss_lob = [vec[ind_lob] for vec in res_xis_lob]

               calc_res_sum_spec_ss_lob, calc_res_xis_spec_ss_lob =
                    GaPSE.sum_ξ_LDxGNC_multipole(COSMO.s_eff, s, COSMO; L=L, alg=:lobatto, kwargs...)

               @test isapprox(res_sum_spec_ss_lob, calc_res_sum_spec_ss_lob; rtol=RTOL)
               @test all([isapprox(a, r; rtol=RTOL) for (a, r) in zip(calc_res_xis_spec_ss_lob, res_xis_spec_ss_lob)])

          end

          @testset "s = 1000, L = 0, with_window" begin
               s = 1000

               ### quad ###
               ind_quad = findfirst(x -> x ≈ s, ss_quad)
               res_sum_spec_ss_quad = res_sums_quad[ind_quad]
               res_xis_spec_ss_quad = [vec[ind_quad] for vec in res_xis_quad]

               calc_res_sum_spec_ss_quad, calc_res_xis_spec_ss_quad =
                    GaPSE.sum_ξ_LDxGNC_multipole(COSMO.s_eff, s, COSMO; L=L, alg=:quad, kwargs...)

               @test isapprox(res_sum_spec_ss_quad, calc_res_sum_spec_ss_quad; rtol=RTOL)
               @test all([isapprox(a, r; rtol=RTOL) for (a, r) in zip(calc_res_xis_spec_ss_quad, res_xis_spec_ss_quad)])

               ### trap ###
               ind_trap = findfirst(x -> x ≈ s, ss_trap)
               res_sum_spec_ss_trap = res_sums_trap[ind_trap]
               res_xis_spec_ss_trap = [vec[ind_trap] for vec in res_xis_trap]

               calc_res_sum_spec_ss_trap, calc_res_xis_spec_ss_trap =
                    GaPSE.sum_ξ_LDxGNC_multipole(COSMO.s_eff, s, COSMO; L=L, alg=:trap, kwargs...)

               @test isapprox(res_sum_spec_ss_trap, calc_res_sum_spec_ss_trap; rtol=RTOL)
               @test all([isapprox(a, r; rtol=RTOL) for (a, r) in zip(calc_res_xis_spec_ss_trap, res_xis_spec_ss_trap)])

               #println("calc_res_xis_spec_ss_trap = $calc_res_xis_spec_ss_trap ;")
               #println("res_xis_spec_ss_trap = $res_xis_spec_ss_trap ;")


               ### lob ###
               ind_lob = findfirst(x -> x ≈ s, ss_lob)
               res_sum_spec_ss_lob = res_sums_lob[ind_lob]
               res_xis_spec_ss_lob = [vec[ind_lob] for vec in res_xis_lob]

               calc_res_sum_spec_ss_lob, calc_res_xis_spec_ss_lob =
                    GaPSE.sum_ξ_LDxGNC_multipole(COSMO.s_eff, s, COSMO; L=L, alg=:lobatto, kwargs...)

               @test isapprox(res_sum_spec_ss_lob, calc_res_sum_spec_ss_lob; rtol=RTOL)
               @test all([isapprox(a, r; rtol=RTOL) for (a, r) in zip(calc_res_xis_spec_ss_lob, res_xis_spec_ss_lob)])

          end
     end

     @testset "L = 1 with_window" begin
          L = 1

          ss_trap, res_sums_trap, res_xis_trap = GaPSE.readxyall(
               "datatest/LDxGNC_SumXiMultipoles/xis_LDxGNC_L$L" * "_trap_withF_specific_ss.txt", comments=true)
          ss_lob, res_sums_lob, res_xis_lob = GaPSE.readxyall(
               "datatest/LDxGNC_SumXiMultipoles/xis_LDxGNC_L$L" * "_lob_withF_specific_ss.txt", comments=true)


          @testset "s = 10, L = 1, with_window" begin
               s = 10

               ### trap ###
               ind_trap = findfirst(x -> x ≈ s, ss_trap)
               res_sum_spec_ss_trap = res_sums_trap[ind_trap]
               res_xis_spec_ss_trap = [vec[ind_trap] for vec in res_xis_trap]

               calc_res_sum_spec_ss_trap, calc_res_xis_spec_ss_trap =
                    GaPSE.sum_ξ_LDxGNC_multipole(COSMO.s_eff, s, COSMO; L=L, alg=:trap, kwargs...)

               @test isapprox(res_sum_spec_ss_trap, calc_res_sum_spec_ss_trap; rtol=RTOL)
               @test all([isapprox(a, r; rtol=RTOL) for (a, r) in zip(calc_res_xis_spec_ss_trap, res_xis_spec_ss_trap)])

               #println("calc_res_xis_spec_ss_trap = $calc_res_xis_spec_ss_trap ;")
               #println("res_xis_spec_ss_trap = $res_xis_spec_ss_trap ;")


               ### lob ###
               ind_lob = findfirst(x -> x ≈ s, ss_lob)
               res_sum_spec_ss_lob = res_sums_lob[ind_lob]
               res_xis_spec_ss_lob = [vec[ind_lob] for vec in res_xis_lob]

               calc_res_sum_spec_ss_lob, calc_res_xis_spec_ss_lob =
                    GaPSE.sum_ξ_LDxGNC_multipole(COSMO.s_eff, s, COSMO; L=L, alg=:lobatto, kwargs...)

               @test isapprox(res_sum_spec_ss_lob, calc_res_sum_spec_ss_lob; rtol=RTOL)
               @test all([isapprox(a, r; rtol=RTOL) for (a, r) in zip(calc_res_xis_spec_ss_lob, res_xis_spec_ss_lob)])

          end

          @testset "s = 500, L = 1, with_window" begin
               s = 500

               ### trap ###
               ind_trap = findfirst(x -> x ≈ s, ss_trap)
               res_sum_spec_ss_trap = res_sums_trap[ind_trap]
               res_xis_spec_ss_trap = [vec[ind_trap] for vec in res_xis_trap]

               calc_res_sum_spec_ss_trap, calc_res_xis_spec_ss_trap =
                    GaPSE.sum_ξ_LDxGNC_multipole(COSMO.s_eff, s, COSMO; L=L, alg=:trap, kwargs...)

               @test isapprox(res_sum_spec_ss_trap, calc_res_sum_spec_ss_trap; rtol=RTOL)
               @test all([isapprox(a, r; rtol=RTOL) for (a, r) in zip(calc_res_xis_spec_ss_trap, res_xis_spec_ss_trap)])

               #println("calc_res_xis_spec_ss_trap = $calc_res_xis_spec_ss_trap ;")
               #println("res_xis_spec_ss_trap = $res_xis_spec_ss_trap ;")


               ### lob ###
               ind_lob = findfirst(x -> x ≈ s, ss_lob)
               res_sum_spec_ss_lob = res_sums_lob[ind_lob]
               res_xis_spec_ss_lob = [vec[ind_lob] for vec in res_xis_lob]

               calc_res_sum_spec_ss_lob, calc_res_xis_spec_ss_lob =
                    GaPSE.sum_ξ_LDxGNC_multipole(COSMO.s_eff, s, COSMO; L=L, alg=:lobatto, kwargs...)

               @test isapprox(res_sum_spec_ss_lob, calc_res_sum_spec_ss_lob; rtol=RTOL)
               @test all([isapprox(a, r; rtol=RTOL) for (a, r) in zip(calc_res_xis_spec_ss_lob, res_xis_spec_ss_lob)])

          end

          @testset "s = 1000, L = 1, with_window" begin
               s = 1000

               ### trap ###
               ind_trap = findfirst(x -> x ≈ s, ss_trap)
               res_sum_spec_ss_trap = res_sums_trap[ind_trap]
               res_xis_spec_ss_trap = [vec[ind_trap] for vec in res_xis_trap]

               calc_res_sum_spec_ss_trap, calc_res_xis_spec_ss_trap =
                    GaPSE.sum_ξ_LDxGNC_multipole(COSMO.s_eff, s, COSMO; L=L, alg=:trap, kwargs...)

               @test isapprox(res_sum_spec_ss_trap, calc_res_sum_spec_ss_trap; rtol=RTOL)
               @test all([isapprox(a, r; rtol=RTOL) for (a, r) in zip(calc_res_xis_spec_ss_trap, res_xis_spec_ss_trap)])

               #println("calc_res_xis_spec_ss_trap = $calc_res_xis_spec_ss_trap ;")
               #println("res_xis_spec_ss_trap = $res_xis_spec_ss_trap ;")


               ### lob ###
               ind_lob = findfirst(x -> x ≈ s, ss_lob)
               res_sum_spec_ss_lob = res_sums_lob[ind_lob]
               res_xis_spec_ss_lob = [vec[ind_lob] for vec in res_xis_lob]

               calc_res_sum_spec_ss_lob, calc_res_xis_spec_ss_lob =
                    GaPSE.sum_ξ_LDxGNC_multipole(COSMO.s_eff, s, COSMO; L=L, alg=:lobatto, kwargs...)

               @test isapprox(res_sum_spec_ss_lob, calc_res_sum_spec_ss_lob; rtol=RTOL)
               @test all([isapprox(a, r; rtol=RTOL) for (a, r) in zip(calc_res_xis_spec_ss_lob, res_xis_spec_ss_lob)])

          end
     end

     @testset "L = 2 with_window" begin
          L = 2
          ss_quad, res_sums_quad, res_xis_quad = GaPSE.readxyall(
               "datatest/LDxGNC_SumXiMultipoles/xis_LDxGNC_L$L" * "_quad_withF_specific_ss.txt", comments=true)
          ss_trap, res_sums_trap, res_xis_trap = GaPSE.readxyall(
               "datatest/LDxGNC_SumXiMultipoles/xis_LDxGNC_L$L" * "_trap_withF_specific_ss.txt", comments=true)
          ss_lob, res_sums_lob, res_xis_lob = GaPSE.readxyall(
               "datatest/LDxGNC_SumXiMultipoles/xis_LDxGNC_L$L" * "_lob_withF_specific_ss.txt", comments=true)


          @testset "s = 10, L = 2, with_window" begin
               s = 10

               ### quad ###
               ind_quad = findfirst(x -> x ≈ s, ss_quad)
               res_sum_spec_ss_quad = res_sums_quad[ind_quad]
               res_xis_spec_ss_quad = [vec[ind_quad] for vec in res_xis_quad]

               calc_res_sum_spec_ss_quad, calc_res_xis_spec_ss_quad =
                    GaPSE.sum_ξ_LDxGNC_multipole(COSMO.s_eff, s, COSMO; L=L, alg=:quad, kwargs...)

               @test isapprox(res_sum_spec_ss_quad, calc_res_sum_spec_ss_quad; rtol=RTOL)
               @test all([isapprox(a, r; rtol=RTOL) for (a, r) in zip(calc_res_xis_spec_ss_quad, res_xis_spec_ss_quad)])

               ### trap ###
               ind_trap = findfirst(x -> x ≈ s, ss_trap)
               res_sum_spec_ss_trap = res_sums_trap[ind_trap]
               res_xis_spec_ss_trap = [vec[ind_trap] for vec in res_xis_trap]

               calc_res_sum_spec_ss_trap, calc_res_xis_spec_ss_trap =
                    GaPSE.sum_ξ_LDxGNC_multipole(COSMO.s_eff, s, COSMO; L=L, alg=:trap, kwargs...)

               @test isapprox(res_sum_spec_ss_trap, calc_res_sum_spec_ss_trap; rtol=RTOL)
               @test all([isapprox(a, r; rtol=RTOL) for (a, r) in zip(calc_res_xis_spec_ss_trap, res_xis_spec_ss_trap)])

               #println("calc_res_xis_spec_ss_trap = $calc_res_xis_spec_ss_trap ;")
               #println("res_xis_spec_ss_trap = $res_xis_spec_ss_trap ;")


               ### lob ###
               ind_lob = findfirst(x -> x ≈ s, ss_lob)
               res_sum_spec_ss_lob = res_sums_lob[ind_lob]
               res_xis_spec_ss_lob = [vec[ind_lob] for vec in res_xis_lob]

               calc_res_sum_spec_ss_lob, calc_res_xis_spec_ss_lob =
                    GaPSE.sum_ξ_LDxGNC_multipole(COSMO.s_eff, s, COSMO; L=L, alg=:lobatto, kwargs...)

               @test isapprox(res_sum_spec_ss_lob, calc_res_sum_spec_ss_lob; rtol=RTOL)
               @test all([isapprox(a, r; rtol=RTOL) for (a, r) in zip(calc_res_xis_spec_ss_lob, res_xis_spec_ss_lob)])

          end

          @testset "s = 500, L = 2, with_window" begin
               s = 500

               ### quad ###
               ind_quad = findfirst(x -> x ≈ s, ss_quad)
               res_sum_spec_ss_quad = res_sums_quad[ind_quad]
               res_xis_spec_ss_quad = [vec[ind_quad] for vec in res_xis_quad]

               calc_res_sum_spec_ss_quad, calc_res_xis_spec_ss_quad =
                    GaPSE.sum_ξ_LDxGNC_multipole(COSMO.s_eff, s, COSMO; L=L, alg=:quad, kwargs...)

               @test isapprox(res_sum_spec_ss_quad, calc_res_sum_spec_ss_quad; rtol=RTOL)
               @test all([isapprox(a, r; rtol=RTOL) for (a, r) in zip(calc_res_xis_spec_ss_quad, res_xis_spec_ss_quad)])

               ### trap ###
               ind_trap = findfirst(x -> x ≈ s, ss_trap)
               res_sum_spec_ss_trap = res_sums_trap[ind_trap]
               res_xis_spec_ss_trap = [vec[ind_trap] for vec in res_xis_trap]

               calc_res_sum_spec_ss_trap, calc_res_xis_spec_ss_trap =
                    GaPSE.sum_ξ_LDxGNC_multipole(COSMO.s_eff, s, COSMO; L=L, alg=:trap, kwargs...)

               @test isapprox(res_sum_spec_ss_trap, calc_res_sum_spec_ss_trap; rtol=RTOL)
               @test all([isapprox(a, r; rtol=RTOL) for (a, r) in zip(calc_res_xis_spec_ss_trap, res_xis_spec_ss_trap)])

               #println("calc_res_xis_spec_ss_trap = $calc_res_xis_spec_ss_trap ;")
               #println("res_xis_spec_ss_trap = $res_xis_spec_ss_trap ;")


               ### lob ###
               ind_lob = findfirst(x -> x ≈ s, ss_lob)
               res_sum_spec_ss_lob = res_sums_lob[ind_lob]
               res_xis_spec_ss_lob = [vec[ind_lob] for vec in res_xis_lob]

               calc_res_sum_spec_ss_lob, calc_res_xis_spec_ss_lob =
                    GaPSE.sum_ξ_LDxGNC_multipole(COSMO.s_eff, s, COSMO; L=L, alg=:lobatto, kwargs...)

               @test isapprox(res_sum_spec_ss_lob, calc_res_sum_spec_ss_lob; rtol=RTOL)
               @test all([isapprox(a, r; rtol=RTOL) for (a, r) in zip(calc_res_xis_spec_ss_lob, res_xis_spec_ss_lob)])

          end

          @testset "s = 1000, L = 2, with_window" begin
               s = 1000

               ### quad ###
               ind_quad = findfirst(x -> x ≈ s, ss_quad)
               res_sum_spec_ss_quad = res_sums_quad[ind_quad]
               res_xis_spec_ss_quad = [vec[ind_quad] for vec in res_xis_quad]

               calc_res_sum_spec_ss_quad, calc_res_xis_spec_ss_quad =
                    GaPSE.sum_ξ_LDxGNC_multipole(COSMO.s_eff, s, COSMO; L=L, alg=:quad, kwargs...)

               @test isapprox(res_sum_spec_ss_quad, calc_res_sum_spec_ss_quad; rtol=RTOL)
               @test all([isapprox(a, r; rtol=RTOL) for (a, r) in zip(calc_res_xis_spec_ss_quad, res_xis_spec_ss_quad)])

               ### trap ###
               ind_trap = findfirst(x -> x ≈ s, ss_trap)
               res_sum_spec_ss_trap = res_sums_trap[ind_trap]
               res_xis_spec_ss_trap = [vec[ind_trap] for vec in res_xis_trap]

               calc_res_sum_spec_ss_trap, calc_res_xis_spec_ss_trap =
                    GaPSE.sum_ξ_LDxGNC_multipole(COSMO.s_eff, s, COSMO; L=L, alg=:trap, kwargs...)

               @test isapprox(res_sum_spec_ss_trap, calc_res_sum_spec_ss_trap; rtol=RTOL)
               @test all([isapprox(a, r; rtol=RTOL) for (a, r) in zip(calc_res_xis_spec_ss_trap, res_xis_spec_ss_trap)])

               #println("calc_res_xis_spec_ss_trap = $calc_res_xis_spec_ss_trap ;")
               #println("res_xis_spec_ss_trap = $res_xis_spec_ss_trap ;")


               ### lob ###
               ind_lob = findfirst(x -> x ≈ s, ss_lob)
               res_sum_spec_ss_lob = res_sums_lob[ind_lob]
               res_xis_spec_ss_lob = [vec[ind_lob] for vec in res_xis_lob]

               calc_res_sum_spec_ss_lob, calc_res_xis_spec_ss_lob =
                    GaPSE.sum_ξ_LDxGNC_multipole(COSMO.s_eff, s, COSMO; L=L, alg=:lobatto, kwargs...)

               @test isapprox(res_sum_spec_ss_lob, calc_res_sum_spec_ss_lob; rtol=RTOL)
               @test all([isapprox(a, r; rtol=RTOL) for (a, r) in zip(calc_res_xis_spec_ss_lob, res_xis_spec_ss_lob)])

          end
     end

     @testset "L = 3 with_window" begin
          L = 3

          ss_trap, res_sums_trap, res_xis_trap = GaPSE.readxyall(
               "datatest/LDxGNC_SumXiMultipoles/xis_LDxGNC_L$L" * "_trap_withF_specific_ss.txt", comments=true)
          ss_lob, res_sums_lob, res_xis_lob = GaPSE.readxyall(
               "datatest/LDxGNC_SumXiMultipoles/xis_LDxGNC_L$L" * "_lob_withF_specific_ss.txt", comments=true)


          @testset "s = 10, L = 1, with_window" begin
               s = 10

               ### trap ###
               ind_trap = findfirst(x -> x ≈ s, ss_trap)
               res_sum_spec_ss_trap = res_sums_trap[ind_trap]
               res_xis_spec_ss_trap = [vec[ind_trap] for vec in res_xis_trap]

               calc_res_sum_spec_ss_trap, calc_res_xis_spec_ss_trap =
                    GaPSE.sum_ξ_LDxGNC_multipole(COSMO.s_eff, s, COSMO; L=L, alg=:trap, kwargs...)

               @test isapprox(res_sum_spec_ss_trap, calc_res_sum_spec_ss_trap; rtol=RTOL)
               @test all([isapprox(a, r; rtol=RTOL) for (a, r) in zip(calc_res_xis_spec_ss_trap, res_xis_spec_ss_trap)])

               #println("calc_res_xis_spec_ss_trap = $calc_res_xis_spec_ss_trap ;")
               #println("res_xis_spec_ss_trap = $res_xis_spec_ss_trap ;")


               ### lob ###
               ind_lob = findfirst(x -> x ≈ s, ss_lob)
               res_sum_spec_ss_lob = res_sums_lob[ind_lob]
               res_xis_spec_ss_lob = [vec[ind_lob] for vec in res_xis_lob]

               calc_res_sum_spec_ss_lob, calc_res_xis_spec_ss_lob =
                    GaPSE.sum_ξ_LDxGNC_multipole(COSMO.s_eff, s, COSMO; L=L, alg=:lobatto, kwargs...)

               @test isapprox(res_sum_spec_ss_lob, calc_res_sum_spec_ss_lob; rtol=RTOL)
               @test all([isapprox(a, r; rtol=RTOL) for (a, r) in zip(calc_res_xis_spec_ss_lob, res_xis_spec_ss_lob)])

          end

          @testset "s = 500, L = 1, with_window" begin
               s = 500

               ### trap ###
               ind_trap = findfirst(x -> x ≈ s, ss_trap)
               res_sum_spec_ss_trap = res_sums_trap[ind_trap]
               res_xis_spec_ss_trap = [vec[ind_trap] for vec in res_xis_trap]

               calc_res_sum_spec_ss_trap, calc_res_xis_spec_ss_trap =
                    GaPSE.sum_ξ_LDxGNC_multipole(COSMO.s_eff, s, COSMO; L=L, alg=:trap, kwargs...)

               @test isapprox(res_sum_spec_ss_trap, calc_res_sum_spec_ss_trap; rtol=RTOL)
               @test all([isapprox(a, r; rtol=RTOL) for (a, r) in zip(calc_res_xis_spec_ss_trap, res_xis_spec_ss_trap)])

               #println("calc_res_xis_spec_ss_trap = $calc_res_xis_spec_ss_trap ;")
               #println("res_xis_spec_ss_trap = $res_xis_spec_ss_trap ;")


               ### lob ###
               ind_lob = findfirst(x -> x ≈ s, ss_lob)
               res_sum_spec_ss_lob = res_sums_lob[ind_lob]
               res_xis_spec_ss_lob = [vec[ind_lob] for vec in res_xis_lob]

               calc_res_sum_spec_ss_lob, calc_res_xis_spec_ss_lob =
                    GaPSE.sum_ξ_LDxGNC_multipole(COSMO.s_eff, s, COSMO; L=L, alg=:lobatto, kwargs...)

               @test isapprox(res_sum_spec_ss_lob, calc_res_sum_spec_ss_lob; rtol=RTOL)
               @test all([isapprox(a, r; rtol=RTOL) for (a, r) in zip(calc_res_xis_spec_ss_lob, res_xis_spec_ss_lob)])

          end

          @testset "s = 1000, L = 1, with_window" begin
               s = 1000

               ### trap ###
               ind_trap = findfirst(x -> x ≈ s, ss_trap)
               res_sum_spec_ss_trap = res_sums_trap[ind_trap]
               res_xis_spec_ss_trap = [vec[ind_trap] for vec in res_xis_trap]

               calc_res_sum_spec_ss_trap, calc_res_xis_spec_ss_trap =
                    GaPSE.sum_ξ_LDxGNC_multipole(COSMO.s_eff, s, COSMO; L=L, alg=:trap, kwargs...)

               @test isapprox(res_sum_spec_ss_trap, calc_res_sum_spec_ss_trap; rtol=RTOL)
               @test all([isapprox(a, r; rtol=RTOL) for (a, r) in zip(calc_res_xis_spec_ss_trap, res_xis_spec_ss_trap)])

               #println("calc_res_xis_spec_ss_trap = $calc_res_xis_spec_ss_trap ;")
               #println("res_xis_spec_ss_trap = $res_xis_spec_ss_trap ;")


               ### lob ###
               ind_lob = findfirst(x -> x ≈ s, ss_lob)
               res_sum_spec_ss_lob = res_sums_lob[ind_lob]
               res_xis_spec_ss_lob = [vec[ind_lob] for vec in res_xis_lob]

               calc_res_sum_spec_ss_lob, calc_res_xis_spec_ss_lob =
                    GaPSE.sum_ξ_LDxGNC_multipole(COSMO.s_eff, s, COSMO; L=L, alg=:lobatto, kwargs...)

               @test isapprox(res_sum_spec_ss_lob, calc_res_sum_spec_ss_lob; rtol=RTOL)
               @test all([isapprox(a, r; rtol=RTOL) for (a, r) in zip(calc_res_xis_spec_ss_lob, res_xis_spec_ss_lob)])

          end
     end

     @testset "L = 4 with_window" begin
          L = 4

          ss_trap, res_sums_trap, res_xis_trap = GaPSE.readxyall(
               "datatest/LDxGNC_SumXiMultipoles/xis_LDxGNC_L$L" * "_trap_withF_specific_ss.txt", comments=true)
          ss_lob, res_sums_lob, res_xis_lob = GaPSE.readxyall(
               "datatest/LDxGNC_SumXiMultipoles/xis_LDxGNC_L$L" * "_lob_withF_specific_ss.txt", comments=true)


          @testset "s = 10, L = 1, with_window" begin
               s = 10

               ### trap ###
               ind_trap = findfirst(x -> x ≈ s, ss_trap)
               res_sum_spec_ss_trap = res_sums_trap[ind_trap]
               res_xis_spec_ss_trap = [vec[ind_trap] for vec in res_xis_trap]

               calc_res_sum_spec_ss_trap, calc_res_xis_spec_ss_trap =
                    GaPSE.sum_ξ_LDxGNC_multipole(COSMO.s_eff, s, COSMO; L=L, alg=:trap, kwargs...)

               @test isapprox(res_sum_spec_ss_trap, calc_res_sum_spec_ss_trap; rtol=RTOL)
               @test all([isapprox(a, r; rtol=RTOL) for (a, r) in zip(calc_res_xis_spec_ss_trap, res_xis_spec_ss_trap)])

               #println("calc_res_xis_spec_ss_trap = $calc_res_xis_spec_ss_trap ;")
               #println("res_xis_spec_ss_trap = $res_xis_spec_ss_trap ;")


               ### lob ###
               ind_lob = findfirst(x -> x ≈ s, ss_lob)
               res_sum_spec_ss_lob = res_sums_lob[ind_lob]
               res_xis_spec_ss_lob = [vec[ind_lob] for vec in res_xis_lob]

               calc_res_sum_spec_ss_lob, calc_res_xis_spec_ss_lob =
                    GaPSE.sum_ξ_LDxGNC_multipole(COSMO.s_eff, s, COSMO; L=L, alg=:lobatto, kwargs...)

               @test isapprox(res_sum_spec_ss_lob, calc_res_sum_spec_ss_lob; rtol=RTOL)
               @test all([isapprox(a, r; rtol=RTOL) for (a, r) in zip(calc_res_xis_spec_ss_lob, res_xis_spec_ss_lob)])

          end

          @testset "s = 500, L = 1, with_window" begin
               s = 500

               ### trap ###
               ind_trap = findfirst(x -> x ≈ s, ss_trap)
               res_sum_spec_ss_trap = res_sums_trap[ind_trap]
               res_xis_spec_ss_trap = [vec[ind_trap] for vec in res_xis_trap]

               calc_res_sum_spec_ss_trap, calc_res_xis_spec_ss_trap =
                    GaPSE.sum_ξ_LDxGNC_multipole(COSMO.s_eff, s, COSMO; L=L, alg=:trap, kwargs...)

               @test isapprox(res_sum_spec_ss_trap, calc_res_sum_spec_ss_trap; rtol=RTOL)
               @test all([isapprox(a, r; rtol=RTOL) for (a, r) in zip(calc_res_xis_spec_ss_trap, res_xis_spec_ss_trap)])

               #println("calc_res_xis_spec_ss_trap = $calc_res_xis_spec_ss_trap ;")
               #println("res_xis_spec_ss_trap = $res_xis_spec_ss_trap ;")


               ### lob ###
               ind_lob = findfirst(x -> x ≈ s, ss_lob)
               res_sum_spec_ss_lob = res_sums_lob[ind_lob]
               res_xis_spec_ss_lob = [vec[ind_lob] for vec in res_xis_lob]

               calc_res_sum_spec_ss_lob, calc_res_xis_spec_ss_lob =
                    GaPSE.sum_ξ_LDxGNC_multipole(COSMO.s_eff, s, COSMO; L=L, alg=:lobatto, kwargs...)

               @test isapprox(res_sum_spec_ss_lob, calc_res_sum_spec_ss_lob; rtol=RTOL)
               @test all([isapprox(a, r; rtol=RTOL) for (a, r) in zip(calc_res_xis_spec_ss_lob, res_xis_spec_ss_lob)])

          end

          @testset "s = 1000, L = 1, with_window" begin
               s = 1000

               ### trap ###
               ind_trap = findfirst(x -> x ≈ s, ss_trap)
               res_sum_spec_ss_trap = res_sums_trap[ind_trap]
               res_xis_spec_ss_trap = [vec[ind_trap] for vec in res_xis_trap]

               calc_res_sum_spec_ss_trap, calc_res_xis_spec_ss_trap =
                    GaPSE.sum_ξ_LDxGNC_multipole(COSMO.s_eff, s, COSMO; L=L, alg=:trap, kwargs...)

               @test isapprox(res_sum_spec_ss_trap, calc_res_sum_spec_ss_trap; rtol=RTOL)
               @test all([isapprox(a, r; rtol=RTOL) for (a, r) in zip(calc_res_xis_spec_ss_trap, res_xis_spec_ss_trap)])

               #println("calc_res_xis_spec_ss_trap = $calc_res_xis_spec_ss_trap ;")
               #println("res_xis_spec_ss_trap = $res_xis_spec_ss_trap ;")


               ### lob ###
               ind_lob = findfirst(x -> x ≈ s, ss_lob)
               res_sum_spec_ss_lob = res_sums_lob[ind_lob]
               res_xis_spec_ss_lob = [vec[ind_lob] for vec in res_xis_lob]

               calc_res_sum_spec_ss_lob, calc_res_xis_spec_ss_lob =
                    GaPSE.sum_ξ_LDxGNC_multipole(COSMO.s_eff, s, COSMO; L=L, alg=:lobatto, kwargs...)

               @test isapprox(res_sum_spec_ss_lob, calc_res_sum_spec_ss_lob; rtol=RTOL)
               @test all([isapprox(a, r; rtol=RTOL) for (a, r) in zip(calc_res_xis_spec_ss_lob, res_xis_spec_ss_lob)])

          end
     end
end

println("\nJust finished the tests on sum_ξ_LDxGNC_multipole!")


##########################################################################################92

println("Now I work on map_sum_ξ_LDxGNC_multipole...")

@testset "test map_sum_ξ_LDxGNC_multipole no_window" begin
     RTOL = 1.2e-2

     kwargs = Dict(
          :use_windows => false,
          :enhancer => 1e8,
          :N_trap => 30, :N_lob => 30,
          :atol_quad => 0.0, :rtol_quad => 1e-2,
          :N_χs => 40, :N_χs_2 => 20,
          :pr => false,
     )

     @testset "L = 0 no_window" begin
          L = 0

          ##### quad ####
          ss_quad, res_sums_quad, res_xis_quad = GaPSE.readxyall(
               "datatest/LDxGNC_SumXiMultipoles/xis_LDxGNC_L$L" * "_quad_noF_specific_ss.txt", comments=true)

          calc_ss_quad, calc_sums_quad, calc_xis_quad = GaPSE.map_sum_ξ_LDxGNC_multipole(COSMO, ss_quad;
               L=L, alg=:quad, kwargs...)

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


          #### trap ####
          ss_trap, res_sums_trap, res_xis_trap = GaPSE.readxyall(
               "datatest/LDxGNC_SumXiMultipoles/xis_LDxGNC_L$L" * "_trap_noF_specific_ss.txt", comments=true)

          calc_ss_trap, calc_sums_trap, calc_xis_trap = GaPSE.map_sum_ξ_LDxGNC_multipole(COSMO, ss_trap;
               L=L, alg=:trap, kwargs...)

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


          #### lob ####
          ss_lob, res_sums_lob, res_xis_lob = GaPSE.readxyall(
               "datatest/LDxGNC_SumXiMultipoles/xis_LDxGNC_L$L" * "_lob_noF_specific_ss.txt", comments=true)

          calc_ss_lob, calc_sums_lob, calc_xis_lob = GaPSE.map_sum_ξ_LDxGNC_multipole(COSMO, ss_lob;
               L=L, alg=:lobatto, kwargs...)

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

     end

     @testset "L = 1 no_window" begin
          L = 1

          #### trap ####
          ss_trap, res_sums_trap, res_xis_trap = GaPSE.readxyall(
               "datatest/LDxGNC_SumXiMultipoles/xis_LDxGNC_L$L" * "_trap_noF_specific_ss.txt", comments=true)

          calc_ss_trap, calc_sums_trap, calc_xis_trap = GaPSE.map_sum_ξ_LDxGNC_multipole(COSMO, ss_trap;
               L=L, alg=:trap, kwargs...)

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


          #### lob ####
          ss_lob, res_sums_lob, res_xis_lob = GaPSE.readxyall(
               "datatest/LDxGNC_SumXiMultipoles/xis_LDxGNC_L$L" * "_lob_noF_specific_ss.txt", comments=true)

          calc_ss_lob, calc_sums_lob, calc_xis_lob = GaPSE.map_sum_ξ_LDxGNC_multipole(COSMO, ss_lob;
               L=L, alg=:lobatto, kwargs...)

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

     end

     @testset "L = 2 no_window" begin
          L = 2

          ##### quad ####
          ss_quad, res_sums_quad, res_xis_quad = GaPSE.readxyall(
               "datatest/LDxGNC_SumXiMultipoles/xis_LDxGNC_L$L" * "_quad_noF_specific_ss.txt", comments=true)

          calc_ss_quad, calc_sums_quad, calc_xis_quad = GaPSE.map_sum_ξ_LDxGNC_multipole(COSMO, ss_quad;
               L=L, alg=:quad, kwargs...)

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


          #### trap ####
          ss_trap, res_sums_trap, res_xis_trap = GaPSE.readxyall(
               "datatest/LDxGNC_SumXiMultipoles/xis_LDxGNC_L$L" * "_trap_noF_specific_ss.txt", comments=true)

          calc_ss_trap, calc_sums_trap, calc_xis_trap = GaPSE.map_sum_ξ_LDxGNC_multipole(COSMO, ss_trap;
               L=L, alg=:trap, kwargs...)

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


          #### lob ####
          ss_lob, res_sums_lob, res_xis_lob = GaPSE.readxyall(
               "datatest/LDxGNC_SumXiMultipoles/xis_LDxGNC_L$L" * "_lob_noF_specific_ss.txt", comments=true)

          calc_ss_lob, calc_sums_lob, calc_xis_lob = GaPSE.map_sum_ξ_LDxGNC_multipole(COSMO, ss_lob;
               L=L, alg=:lobatto, kwargs...)

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

     end

     @testset "L = 3 no_window" begin
          L = 3


          #### trap ####
          ss_trap, res_sums_trap, res_xis_trap = GaPSE.readxyall(
               "datatest/LDxGNC_SumXiMultipoles/xis_LDxGNC_L$L" * "_trap_noF_specific_ss.txt", comments=true)

          calc_ss_trap, calc_sums_trap, calc_xis_trap = GaPSE.map_sum_ξ_LDxGNC_multipole(COSMO, ss_trap;
               L=L, alg=:trap, kwargs...)

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


          #### lob ####
          ss_lob, res_sums_lob, res_xis_lob = GaPSE.readxyall(
               "datatest/LDxGNC_SumXiMultipoles/xis_LDxGNC_L$L" * "_lob_noF_specific_ss.txt", comments=true)

          calc_ss_lob, calc_sums_lob, calc_xis_lob = GaPSE.map_sum_ξ_LDxGNC_multipole(COSMO, ss_lob;
               L=L, alg=:lobatto, kwargs...)

          #println("res_xis_lob[2] = $(res_xis_lob[2]) ;")
          #println("calc_xis_lob[2] = $(calc_xis_lob[2]); ")

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

     end

     @testset "L = 4 no_window" begin
          L = 4

          #### trap ####
          ss_trap, res_sums_trap, res_xis_trap = GaPSE.readxyall(
               "datatest/LDxGNC_SumXiMultipoles/xis_LDxGNC_L$L" * "_trap_noF_specific_ss.txt", comments=true)

          calc_ss_trap, calc_sums_trap, calc_xis_trap = GaPSE.map_sum_ξ_LDxGNC_multipole(COSMO, ss_trap;
               L=L, alg=:trap, kwargs...)

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


          #### lob ####
          ss_lob, res_sums_lob, res_xis_lob = GaPSE.readxyall(
               "datatest/LDxGNC_SumXiMultipoles/xis_LDxGNC_L$L" * "_lob_noF_specific_ss.txt", comments=true)

          calc_ss_lob, calc_sums_lob, calc_xis_lob = GaPSE.map_sum_ξ_LDxGNC_multipole(COSMO, ss_lob;
               L=L, alg=:lobatto, kwargs...)

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

     end
end

println("\nDon't worry, I am fully operative...")

@testset "test map_sum_ξ_LDxGNC_multipole with_window" begin
     RTOL = 1.2e-2

     kwargs = Dict(
          :use_windows => true,
          :enhancer => 1e8,
          :N_trap => 30, :N_lob => 30,
          :atol_quad => 0.0, :rtol_quad => 1e-2,
          :N_χs => 40, :N_χs_2 => 20,
          :pr => false,
     )

     @testset "L = 0 with_window" begin
          L = 0

          ##### quad ####
          ss_quad, res_sums_quad, res_xis_quad = GaPSE.readxyall(
               "datatest/LDxGNC_SumXiMultipoles/xis_LDxGNC_L$L" * "_quad_withF_specific_ss.txt", comments=true)

          calc_ss_quad, calc_sums_quad, calc_xis_quad = GaPSE.map_sum_ξ_LDxGNC_multipole(COSMO, ss_quad;
               L=L, alg=:quad, kwargs...)

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


          #### trap ####
          ss_trap, res_sums_trap, res_xis_trap = GaPSE.readxyall(
               "datatest/LDxGNC_SumXiMultipoles/xis_LDxGNC_L$L" * "_trap_withF_specific_ss.txt", comments=true)

          calc_ss_trap, calc_sums_trap, calc_xis_trap = GaPSE.map_sum_ξ_LDxGNC_multipole(COSMO, ss_trap;
               L=L, alg=:trap, kwargs...)

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


          #### lob ####
          ss_lob, res_sums_lob, res_xis_lob = GaPSE.readxyall(
               "datatest/LDxGNC_SumXiMultipoles/xis_LDxGNC_L$L" * "_lob_withF_specific_ss.txt", comments=true)

          calc_ss_lob, calc_sums_lob, calc_xis_lob = GaPSE.map_sum_ξ_LDxGNC_multipole(COSMO, ss_lob;
               L=L, alg=:lobatto, kwargs...)

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

     end

     @testset "L = 1 with_window" begin
          L = 1


          #### trap ####
          ss_trap, res_sums_trap, res_xis_trap = GaPSE.readxyall(
               "datatest/LDxGNC_SumXiMultipoles/xis_LDxGNC_L$L" * "_trap_withF_specific_ss.txt", comments=true)

          calc_ss_trap, calc_sums_trap, calc_xis_trap = GaPSE.map_sum_ξ_LDxGNC_multipole(COSMO, ss_trap;
               L=L, alg=:trap, kwargs...)

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


          #### lob ####
          ss_lob, res_sums_lob, res_xis_lob = GaPSE.readxyall(
               "datatest/LDxGNC_SumXiMultipoles/xis_LDxGNC_L$L" * "_lob_withF_specific_ss.txt", comments=true)

          calc_ss_lob, calc_sums_lob, calc_xis_lob = GaPSE.map_sum_ξ_LDxGNC_multipole(COSMO, ss_lob;
               L=L, alg=:lobatto, kwargs...)

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

     end

     @testset "L = 2 with_window" begin
          L = 2

          ##### quad ####
          ss_quad, res_sums_quad, res_xis_quad = GaPSE.readxyall(
               "datatest/LDxGNC_SumXiMultipoles/xis_LDxGNC_L$L" * "_quad_withF_specific_ss.txt", comments=true)

          calc_ss_quad, calc_sums_quad, calc_xis_quad = GaPSE.map_sum_ξ_LDxGNC_multipole(COSMO, ss_quad;
               L=L, alg=:quad, kwargs...)

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


          #### trap ####
          ss_trap, res_sums_trap, res_xis_trap = GaPSE.readxyall(
               "datatest/LDxGNC_SumXiMultipoles/xis_LDxGNC_L$L" * "_trap_withF_specific_ss.txt", comments=true)

          calc_ss_trap, calc_sums_trap, calc_xis_trap = GaPSE.map_sum_ξ_LDxGNC_multipole(COSMO, ss_trap;
               L=L, alg=:trap, kwargs...)

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


          #### lob ####
          ss_lob, res_sums_lob, res_xis_lob = GaPSE.readxyall(
               "datatest/LDxGNC_SumXiMultipoles/xis_LDxGNC_L$L" * "_lob_withF_specific_ss.txt", comments=true)

          calc_ss_lob, calc_sums_lob, calc_xis_lob = GaPSE.map_sum_ξ_LDxGNC_multipole(COSMO, ss_lob;
               L=L, alg=:lobatto, kwargs...)

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

     end

     @testset "L = 3 with_window" begin
          L = 3

          #### trap ####
          ss_trap, res_sums_trap, res_xis_trap = GaPSE.readxyall(
               "datatest/LDxGNC_SumXiMultipoles/xis_LDxGNC_L$L" * "_trap_withF_specific_ss.txt", comments=true)

          calc_ss_trap, calc_sums_trap, calc_xis_trap = GaPSE.map_sum_ξ_LDxGNC_multipole(COSMO, ss_trap;
               L=L, alg=:trap, kwargs...)

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


          #### lob ####
          ss_lob, res_sums_lob, res_xis_lob = GaPSE.readxyall(
               "datatest/LDxGNC_SumXiMultipoles/xis_LDxGNC_L$L" * "_lob_withF_specific_ss.txt", comments=true)

          calc_ss_lob, calc_sums_lob, calc_xis_lob = GaPSE.map_sum_ξ_LDxGNC_multipole(COSMO, ss_lob;
               L=L, alg=:lobatto, kwargs...)

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

     end

     @testset "L = 4 with_window" begin
          L = 4
          
          #### trap ####
          ss_trap, res_sums_trap, res_xis_trap = GaPSE.readxyall(
               "datatest/LDxGNC_SumXiMultipoles/xis_LDxGNC_L$L" * "_trap_withF_specific_ss.txt", comments=true)

          calc_ss_trap, calc_sums_trap, calc_xis_trap = GaPSE.map_sum_ξ_LDxGNC_multipole(COSMO, ss_trap;
               L=L, alg=:trap, kwargs...)

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


          #### lob ####
          ss_lob, res_sums_lob, res_xis_lob = GaPSE.readxyall(
               "datatest/LDxGNC_SumXiMultipoles/xis_LDxGNC_L$L" * "_lob_withF_specific_ss.txt", comments=true)

          calc_ss_lob, calc_sums_lob, calc_xis_lob = GaPSE.map_sum_ξ_LDxGNC_multipole(COSMO, ss_lob;
               L=L, alg=:lobatto, kwargs...)

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

     end
end

println("\nJust finished the tests on map_sum_ξ_LDxGNC_multipole!")






