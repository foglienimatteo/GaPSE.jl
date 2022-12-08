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

println("\n\nNow I start the tests of ξ_GNC_multipole.")
println("It will take a while, but do not worry: I'm working.")

@testset "test ξ_GNC_multipole" begin
     name_effect = "auto_doppler"
     func_effect = GaPSE.ξ_GNC_Doppler
     RTOL = 1.2e-2
     #SS_GNC_DOPPLER_WITHOBS = 10 .^ range(-1, log10(2 * COSMO.s_max), length=300);

     @testset "zeros" begin
          @test_throws AssertionError GaPSE.ξ_GNC_multipole(COSMO.s_eff, 10.0, "strange", COSMO)
          @test_throws AssertionError GaPSE.ξ_GNC_multipole(COSMO.s_eff, 10.0, "auto_doppler", COSMO;
               obs=:trap)
          @test_throws AssertionError GaPSE.ξ_GNC_multipole(COSMO.s_eff, 10.0, "auto_doppler", COSMO;
               alg=:noobsvel)

          @test_throws AssertionError GaPSE.ξ_GNC_multipole(COSMO.s_eff, 10.0, "auto_doppler", COSMO;
               L=-1)
          @test_throws AssertionError GaPSE.ξ_GNC_multipole(COSMO.s_eff, 10.0, "auto_doppler", COSMO;
               atol_quad=-0.1)
          @test_throws AssertionError GaPSE.ξ_GNC_multipole(COSMO.s_eff, 10.0, "auto_doppler", COSMO;
               rtol_quad=-0.1)
          @test_throws AssertionError GaPSE.ξ_GNC_multipole(COSMO.s_eff, 10.0, "auto_doppler", COSMO;
               N_trap=2)
          @test_throws AssertionError GaPSE.ξ_GNC_multipole(COSMO.s_eff, 10.0, "auto_doppler", COSMO;
               N_lob=2)
     end


     @testset "no window, with observer terms" begin
          kwargs = Dict(
               :enhancer => 1e8, :use_windows => false,
               :N_lob => 30, :N_trap => 30, :obs => :yes,
               :atol_quad => 0.0, :rtol_quad => 1e-2,
          )

          @testset "L = 0" begin
               L = 0
               ss_lob, xis_lob = GaPSE.readxy("datatest/GNC_doppler_multipoles/xi_GNC_" *
                                              name_effect * "_lobatto_noF_withobs_L$L" * ".txt"; comments=true)
               ss_trap, xis_trap = GaPSE.readxy("datatest/GNC_doppler_multipoles/xi_GNC_" *
                                                name_effect * "_trap_noF_withobs_L$L" * ".txt"; comments=true)

               calc_xis_lob = [GaPSE.ξ_GNC_multipole(COSMO.s_eff, s, name_effect, COSMO;
                    L=L, alg=:lobatto, kwargs...) for s in ss_lob]
               calc_xis_trap = [GaPSE.ξ_GNC_multipole(COSMO.s_eff, s, name_effect, COSMO;
                    L=L, alg=:trap, kwargs...) for s in ss_trap]

               @test all([isapprox(xi, calc_xi, rtol=RTOL) for (xi, calc_xi) in zip(xis_lob, calc_xis_lob)])
               @test all([isapprox(xi, calc_xi, rtol=RTOL) for (xi, calc_xi) in zip(xis_trap, calc_xis_trap)])

               ss_quad, xis_quad = GaPSE.readxy("datatest/GNC_doppler_multipoles/xi_GNC_" *
                                                name_effect * "_quad_noF_withobs_L$L" * ".txt"; comments=true)
               calc_xis_quad = [GaPSE.ξ_GNC_multipole(COSMO.s_eff, s, name_effect, COSMO;
                    L=L, alg=:quad, kwargs...) for s in ss_quad]

               @test all([isapprox(xi, calc_xi, rtol=RTOL) for (xi, calc_xi) in zip(xis_quad, calc_xis_quad)])
          end

          @testset "L = 1" begin
               L = 1
               ss_lob, xis_lob = GaPSE.readxy("datatest/GNC_doppler_multipoles/xi_GNC_" *
                                              name_effect * "_lobatto_noF_withobs_L$L" * ".txt"; comments=true)
               ss_trap, xis_trap = GaPSE.readxy("datatest/GNC_doppler_multipoles/xi_GNC_" *
                                                name_effect * "_trap_noF_withobs_L$L" * ".txt"; comments=true)

               calc_xis_lob = [GaPSE.ξ_GNC_multipole(COSMO.s_eff, s, name_effect, COSMO;
                    L=L, alg=:lobatto, kwargs...) for s in ss_lob]
               calc_xis_trap = [GaPSE.ξ_GNC_multipole(COSMO.s_eff, s, name_effect, COSMO;
                    L=L, alg=:trap, kwargs...) for s in ss_trap]

               @test all([isapprox(xi, calc_xi, rtol=RTOL) for (xi, calc_xi) in zip(xis_lob, calc_xis_lob)])
               @test all([isapprox(xi, calc_xi, rtol=RTOL) for (xi, calc_xi) in zip(xis_trap, calc_xis_trap)])
          end

          @testset "L = 2" begin
               L = 2
               ss_lob, xis_lob = GaPSE.readxy("datatest/GNC_doppler_multipoles/xi_GNC_" *
                                              name_effect * "_lobatto_noF_withobs_L$L" * ".txt"; comments=true)
               ss_trap, xis_trap = GaPSE.readxy("datatest/GNC_doppler_multipoles/xi_GNC_" *
                                                name_effect * "_trap_noF_withobs_L$L" * ".txt"; comments=true)

               calc_xis_lob = [GaPSE.ξ_GNC_multipole(COSMO.s_eff, s, name_effect, COSMO;
                    L=L, alg=:lobatto, kwargs...) for s in ss_lob]
               calc_xis_trap = [GaPSE.ξ_GNC_multipole(COSMO.s_eff, s, name_effect, COSMO;
                    L=L, alg=:trap, kwargs...) for s in ss_trap]

               @test all([isapprox(xi, calc_xi, rtol=RTOL) for (xi, calc_xi) in zip(xis_lob, calc_xis_lob)])
               @test all([isapprox(xi, calc_xi, rtol=RTOL) for (xi, calc_xi) in zip(xis_trap, calc_xis_trap)])


               ss_quad, xis_quad = GaPSE.readxy("datatest/GNC_doppler_multipoles/xi_GNC_" *
                                                name_effect * "_quad_noF_withobs_L$L" * ".txt"; comments=true)
               calc_xis_quad = [GaPSE.ξ_GNC_multipole(COSMO.s_eff, s, name_effect, COSMO;
                    L=L, alg=:quad, kwargs...) for s in ss_quad]

               @test all([isapprox(xi, calc_xi, rtol=RTOL) for (xi, calc_xi) in zip(xis_quad, calc_xis_quad)])
          end

          @testset "L = 3" begin
               L = 3
               ss_lob, xis_lob = GaPSE.readxy("datatest/GNC_doppler_multipoles/xi_GNC_" *
                                              name_effect * "_lobatto_noF_withobs_L$L" * ".txt"; comments=true)
               ss_trap, xis_trap = GaPSE.readxy("datatest/GNC_doppler_multipoles/xi_GNC_" *
                                                name_effect * "_trap_noF_withobs_L$L" * ".txt"; comments=true)

               calc_xis_lob = [GaPSE.ξ_GNC_multipole(COSMO.s_eff, s, name_effect, COSMO;
                    L=L, alg=:lobatto, kwargs...) for s in ss_lob]
               calc_xis_trap = [GaPSE.ξ_GNC_multipole(COSMO.s_eff, s, name_effect, COSMO;
                    L=L, alg=:trap, kwargs...) for s in ss_trap]

               @test all([isapprox(xi, calc_xi, rtol=RTOL) for (xi, calc_xi) in zip(xis_lob, calc_xis_lob)])
               @test all([isapprox(xi, calc_xi, rtol=RTOL) for (xi, calc_xi) in zip(xis_trap, calc_xis_trap)])

          end

          @testset "L = 4" begin
               L = 4
               ss_lob, xis_lob = GaPSE.readxy("datatest/GNC_doppler_multipoles/xi_GNC_" *
                                              name_effect * "_lobatto_noF_withobs_L$L" * ".txt"; comments=true)
               ss_trap, xis_trap = GaPSE.readxy("datatest/GNC_doppler_multipoles/xi_GNC_" *
                                                name_effect * "_trap_noF_withobs_L$L" * ".txt"; comments=true)

               calc_xis_lob = [GaPSE.ξ_GNC_multipole(COSMO.s_eff, s, name_effect, COSMO;
                    L=L, alg=:lobatto, kwargs...) for s in ss_lob]
               calc_xis_trap = [GaPSE.ξ_GNC_multipole(COSMO.s_eff, s, name_effect, COSMO;
                    L=L, alg=:trap, kwargs...) for s in ss_trap]

               @test all([isapprox(xi, calc_xi, rtol=RTOL) for (xi, calc_xi) in zip(xis_lob, calc_xis_lob)])
               @test all([isapprox(xi, calc_xi, rtol=RTOL) for (xi, calc_xi) in zip(xis_trap, calc_xis_trap)])

               #println("xis_trap = $xis_trap;")
               #println("calc_xis_trap = $calc_xis_trap ;")
          end
     end

     @testset "no window, no observer velocity terms" begin
          kwargs = Dict(
               :enhancer => 1e8, :use_windows => false,
               :N_lob => 30, :N_trap => 30, :obs => :noobsvel,
               :atol_quad => 0.0, :rtol_quad => 1e-2,
          )

          @testset "L = 0" begin
               L = 0
               ss_lob, xis_lob = GaPSE.readxy("datatest/GNC_doppler_multipoles/xi_GNC_" *
                                              name_effect * "_lobatto_noF_noobsvel_L$L" * ".txt"; comments=true)
               ss_trap, xis_trap = GaPSE.readxy("datatest/GNC_doppler_multipoles/xi_GNC_" *
                                                name_effect * "_trap_noF_noobsvel_L$L" * ".txt"; comments=true)

               calc_xis_lob = [GaPSE.ξ_GNC_multipole(COSMO.s_eff, s, name_effect, COSMO;
                    L=L, alg=:lobatto, kwargs...) for s in ss_lob]
               calc_xis_trap = [GaPSE.ξ_GNC_multipole(COSMO.s_eff, s, name_effect, COSMO;
                    L=L, alg=:trap, kwargs...) for s in ss_trap]

               @test all([isapprox(xi, calc_xi, rtol=RTOL) for (xi, calc_xi) in zip(xis_lob, calc_xis_lob)])
               @test all([isapprox(xi, calc_xi, rtol=RTOL) for (xi, calc_xi) in zip(xis_trap, calc_xis_trap)])


               ss_quad, xis_quad = GaPSE.readxy("datatest/GNC_doppler_multipoles/xi_GNC_" *
                                                name_effect * "_quad_noF_noobsvel_L$L" * ".txt"; comments=true)
               calc_xis_quad = [GaPSE.ξ_GNC_multipole(COSMO.s_eff, s, name_effect, COSMO;
                    L=L, alg=:quad, kwargs...) for s in ss_quad]

               @test all([isapprox(xi, calc_xi, rtol=RTOL) for (xi, calc_xi) in zip(xis_quad, calc_xis_quad)])
          end

          @testset "L = 1" begin
               L = 1
               ss_lob, xis_lob = GaPSE.readxy("datatest/GNC_doppler_multipoles/xi_GNC_" *
                                              name_effect * "_lobatto_noF_noobsvel_L$L" * ".txt"; comments=true)
               ss_trap, xis_trap = GaPSE.readxy("datatest/GNC_doppler_multipoles/xi_GNC_" *
                                                name_effect * "_trap_noF_noobsvel_L$L" * ".txt"; comments=true)

               calc_xis_lob = [GaPSE.ξ_GNC_multipole(COSMO.s_eff, s, name_effect, COSMO;
                    L=L, alg=:lobatto, kwargs...) for s in ss_lob]
               calc_xis_trap = [GaPSE.ξ_GNC_multipole(COSMO.s_eff, s, name_effect, COSMO;
                    L=L, alg=:trap, kwargs...) for s in ss_trap]

               @test all([isapprox(xi, calc_xi, rtol=RTOL) for (xi, calc_xi) in zip(xis_lob, calc_xis_lob)])
               @test all([isapprox(xi, calc_xi, rtol=RTOL) for (xi, calc_xi) in zip(xis_trap, calc_xis_trap)])
          end

          @testset "L = 2" begin
               L = 2
               ss_lob, xis_lob = GaPSE.readxy("datatest/GNC_doppler_multipoles/xi_GNC_" *
                                              name_effect * "_lobatto_noF_noobsvel_L$L" * ".txt"; comments=true)
               ss_trap, xis_trap = GaPSE.readxy("datatest/GNC_doppler_multipoles/xi_GNC_" *
                                                name_effect * "_trap_noF_noobsvel_L$L" * ".txt"; comments=true)

               calc_xis_lob = [GaPSE.ξ_GNC_multipole(COSMO.s_eff, s, name_effect, COSMO;
                    L=L, alg=:lobatto, kwargs...) for s in ss_lob]
               calc_xis_trap = [GaPSE.ξ_GNC_multipole(COSMO.s_eff, s, name_effect, COSMO;
                    L=L, alg=:trap, kwargs...) for s in ss_trap]

               @test all([isapprox(xi, calc_xi, rtol=RTOL) for (xi, calc_xi) in zip(xis_lob, calc_xis_lob)])
               @test all([isapprox(xi, calc_xi, rtol=RTOL) for (xi, calc_xi) in zip(xis_trap, calc_xis_trap)])

               ss_quad, xis_quad = GaPSE.readxy("datatest/GNC_doppler_multipoles/xi_GNC_" *
                                                name_effect * "_quad_noF_noobsvel_L$L" * ".txt"; comments=true)
               calc_xis_quad = [GaPSE.ξ_GNC_multipole(COSMO.s_eff, s, name_effect, COSMO;
                    L=L, alg=:quad, kwargs...) for s in ss_quad]

               @test all([isapprox(xi, calc_xi, rtol=RTOL) for (xi, calc_xi) in zip(xis_quad, calc_xis_quad)])
               #println("xis_trap = $xis_trap;")
               #println("calc_xis_trap = $calc_xis_trap ;")
          end

          @testset "L = 3" begin
               L = 3
               ss_lob, xis_lob = GaPSE.readxy("datatest/GNC_doppler_multipoles/xi_GNC_" *
                                              name_effect * "_lobatto_noF_noobsvel_L$L" * ".txt"; comments=true)
               ss_trap, xis_trap = GaPSE.readxy("datatest/GNC_doppler_multipoles/xi_GNC_" *
                                                name_effect * "_trap_noF_noobsvel_L$L" * ".txt"; comments=true)

               calc_xis_lob = [GaPSE.ξ_GNC_multipole(COSMO.s_eff, s, name_effect, COSMO;
                    L=L, alg=:lobatto, kwargs...) for s in ss_lob]
               calc_xis_trap = [GaPSE.ξ_GNC_multipole(COSMO.s_eff, s, name_effect, COSMO;
                    L=L, alg=:trap, kwargs...) for s in ss_trap]

               @test all([isapprox(xi, calc_xi, rtol=RTOL) for (xi, calc_xi) in zip(xis_lob, calc_xis_lob)])
               @test all([isapprox(xi, calc_xi, rtol=RTOL) for (xi, calc_xi) in zip(xis_trap, calc_xis_trap)])
          end

          @testset "L = 4" begin
               L = 4
               ss_lob, xis_lob = GaPSE.readxy("datatest/GNC_doppler_multipoles/xi_GNC_" *
                                              name_effect * "_lobatto_noF_noobsvel_L$L" * ".txt"; comments=true)
               ss_trap, xis_trap = GaPSE.readxy("datatest/GNC_doppler_multipoles/xi_GNC_" *
                                                name_effect * "_trap_noF_noobsvel_L$L" * ".txt"; comments=true)

               calc_xis_lob = [GaPSE.ξ_GNC_multipole(COSMO.s_eff, s, name_effect, COSMO;
                    L=L, alg=:lobatto, kwargs...) for s in ss_lob]
               calc_xis_trap = [GaPSE.ξ_GNC_multipole(COSMO.s_eff, s, name_effect, COSMO;
                    L=L, alg=:trap, kwargs...) for s in ss_trap]

               @test all([isapprox(xi, calc_xi, rtol=RTOL) for (xi, calc_xi) in zip(xis_lob, calc_xis_lob)])
               @test all([isapprox(xi, calc_xi, rtol=RTOL) for (xi, calc_xi) in zip(xis_trap, calc_xis_trap)])
          end
     end

     @testset "no window, no observer terms" begin
          kwargs = Dict(
               :enhancer => 1e8, :use_windows => false,
               :N_lob => 30, :N_trap => 30, :obs => :no,
               :atol_quad => 0.0, :rtol_quad => 1e-2,
          )

          @testset "L = 0" begin
               L = 0
               ss_lob, xis_lob = GaPSE.readxy("datatest/GNC_doppler_multipoles/xi_GNC_" *
                                              name_effect * "_lobatto_noF_noobs_L$L" * ".txt"; comments=true)
               ss_trap, xis_trap = GaPSE.readxy("datatest/GNC_doppler_multipoles/xi_GNC_" *
                                                name_effect * "_trap_noF_noobs_L$L" * ".txt"; comments=true)

               calc_xis_lob = [GaPSE.ξ_GNC_multipole(COSMO.s_eff, s, name_effect, COSMO;
                    L=L, alg=:lobatto, kwargs...) for s in ss_lob]
               calc_xis_trap = [GaPSE.ξ_GNC_multipole(COSMO.s_eff, s, name_effect, COSMO;
                    L=L, alg=:trap, kwargs...) for s in ss_trap]

               @test all([isapprox(xi, calc_xi, rtol=RTOL) for (xi, calc_xi) in zip(xis_lob, calc_xis_lob)])
               @test all([isapprox(xi, calc_xi, rtol=RTOL) for (xi, calc_xi) in zip(xis_trap, calc_xis_trap)])


               ss_quad, xis_quad = GaPSE.readxy("datatest/GNC_doppler_multipoles/xi_GNC_" *
                                                name_effect * "_quad_noF_noobs_L$L" * ".txt"; comments=true)
               calc_xis_quad = [GaPSE.ξ_GNC_multipole(COSMO.s_eff, s, name_effect, COSMO;
                    L=L, alg=:quad, kwargs...) for s in ss_quad]

               @test all([isapprox(xi, calc_xi, rtol=RTOL) for (xi, calc_xi) in zip(xis_quad, calc_xis_quad)])
          end

          @testset "L = 1" begin
               L = 1
               ss_lob, xis_lob = GaPSE.readxy("datatest/GNC_doppler_multipoles/xi_GNC_" *
                                              name_effect * "_lobatto_noF_noobs_L$L" * ".txt"; comments=true)
               ss_trap, xis_trap = GaPSE.readxy("datatest/GNC_doppler_multipoles/xi_GNC_" *
                                                name_effect * "_trap_noF_noobs_L$L" * ".txt"; comments=true)

               calc_xis_lob = [GaPSE.ξ_GNC_multipole(COSMO.s_eff, s, name_effect, COSMO;
                    L=L, alg=:lobatto, kwargs...) for s in ss_lob]
               calc_xis_trap = [GaPSE.ξ_GNC_multipole(COSMO.s_eff, s, name_effect, COSMO;
                    L=L, alg=:trap, kwargs...) for s in ss_trap]

               @test all([isapprox(xi, calc_xi, rtol=RTOL) for (xi, calc_xi) in zip(xis_lob, calc_xis_lob)])
               @test all([isapprox(xi, calc_xi, rtol=RTOL) for (xi, calc_xi) in zip(xis_trap, calc_xis_trap)])
          end

          @testset "L = 2" begin
               L = 2
               ss_lob, xis_lob = GaPSE.readxy("datatest/GNC_doppler_multipoles/xi_GNC_" *
                                              name_effect * "_lobatto_noF_noobs_L$L" * ".txt"; comments=true)
               ss_trap, xis_trap = GaPSE.readxy("datatest/GNC_doppler_multipoles/xi_GNC_" *
                                                name_effect * "_trap_noF_noobs_L$L" * ".txt"; comments=true)

               calc_xis_lob = [GaPSE.ξ_GNC_multipole(COSMO.s_eff, s, name_effect, COSMO;
                    L=L, alg=:lobatto, kwargs...) for s in ss_lob]
               calc_xis_trap = [GaPSE.ξ_GNC_multipole(COSMO.s_eff, s, name_effect, COSMO;
                    L=L, alg=:trap, kwargs...) for s in ss_trap]

               @test all([isapprox(xi, calc_xi, rtol=RTOL) for (xi, calc_xi) in zip(xis_lob, calc_xis_lob)])
               @test all([isapprox(xi, calc_xi, rtol=RTOL) for (xi, calc_xi) in zip(xis_trap, calc_xis_trap)])

               ss_quad, xis_quad = GaPSE.readxy("datatest/GNC_doppler_multipoles/xi_GNC_" *
                                                name_effect * "_quad_noF_noobs_L$L" * ".txt"; comments=true)
               calc_xis_quad = [GaPSE.ξ_GNC_multipole(COSMO.s_eff, s, name_effect, COSMO;
                    L=L, alg=:quad, kwargs...) for s in ss_quad]

               @test all([isapprox(xi, calc_xi, rtol=RTOL) for (xi, calc_xi) in zip(xis_quad, calc_xis_quad)])
               #println("xis_trap = $xis_trap;")
               #println("calc_xis_trap = $calc_xis_trap ;")
          end

          @testset "L = 3" begin
               L = 3
               ss_lob, xis_lob = GaPSE.readxy("datatest/GNC_doppler_multipoles/xi_GNC_" *
                                              name_effect * "_lobatto_noF_noobs_L$L" * ".txt"; comments=true)
               ss_trap, xis_trap = GaPSE.readxy("datatest/GNC_doppler_multipoles/xi_GNC_" *
                                                name_effect * "_trap_noF_noobs_L$L" * ".txt"; comments=true)

               calc_xis_lob = [GaPSE.ξ_GNC_multipole(COSMO.s_eff, s, name_effect, COSMO;
                    L=L, alg=:lobatto, kwargs...) for s in ss_lob]
               calc_xis_trap = [GaPSE.ξ_GNC_multipole(COSMO.s_eff, s, name_effect, COSMO;
                    L=L, alg=:trap, kwargs...) for s in ss_trap]

               @test all([isapprox(xi, calc_xi, rtol=RTOL) for (xi, calc_xi) in zip(xis_lob, calc_xis_lob)])
               @test all([isapprox(xi, calc_xi, rtol=RTOL) for (xi, calc_xi) in zip(xis_trap, calc_xis_trap)])
          end

          @testset "L = 4" begin
               L = 4
               ss_lob, xis_lob = GaPSE.readxy("datatest/GNC_doppler_multipoles/xi_GNC_" *
                                              name_effect * "_lobatto_noF_noobs_L$L" * ".txt"; comments=true)
               ss_trap, xis_trap = GaPSE.readxy("datatest/GNC_doppler_multipoles/xi_GNC_" *
                                                name_effect * "_trap_noF_noobs_L$L" * ".txt"; comments=true)

               calc_xis_lob = [GaPSE.ξ_GNC_multipole(COSMO.s_eff, s, name_effect, COSMO;
                    L=L, alg=:lobatto, kwargs...) for s in ss_lob]
               calc_xis_trap = [GaPSE.ξ_GNC_multipole(COSMO.s_eff, s, name_effect, COSMO;
                    L=L, alg=:trap, kwargs...) for s in ss_trap]

               @test all([isapprox(xi, calc_xi, rtol=RTOL) for (xi, calc_xi) in zip(xis_lob, calc_xis_lob)])
               @test all([isapprox(xi, calc_xi, rtol=RTOL) for (xi, calc_xi) in zip(xis_trap, calc_xis_trap)])
          end
     end

     @testset "with window, no observer terms" begin
          kwargs = Dict(
               :enhancer => 1e8, :use_windows => true,
               :N_lob => 30, :N_trap => 30, :obs => :no,
               :atol_quad => 0.0, :rtol_quad => 1e-2,
          )

          @testset "L = 0" begin
               L = 0
               ss_lob, xis_lob = GaPSE.readxy("datatest/GNC_doppler_multipoles/xi_GNC_" *
                                              name_effect * "_lobatto_withF_noobs_L$L" * ".txt"; comments=true)
               ss_trap, xis_trap = GaPSE.readxy("datatest/GNC_doppler_multipoles/xi_GNC_" *
                                                name_effect * "_trap_withF_noobs_L$L" * ".txt"; comments=true)

               calc_xis_lob = [GaPSE.ξ_GNC_multipole(COSMO.s_eff, s, name_effect, COSMO;
                    L=L, alg=:lobatto, kwargs...) for s in ss_lob]
               calc_xis_trap = [GaPSE.ξ_GNC_multipole(COSMO.s_eff, s, name_effect, COSMO;
                    L=L, alg=:trap, kwargs...) for s in ss_trap]

               @test all([isapprox(xi, calc_xi, rtol=RTOL) for (xi, calc_xi) in zip(xis_lob, calc_xis_lob)])
               @test all([isapprox(xi, calc_xi, rtol=RTOL) for (xi, calc_xi) in zip(xis_trap, calc_xis_trap)])


               ss_quad, xis_quad = GaPSE.readxy("datatest/GNC_doppler_multipoles/xi_GNC_" *
                                                name_effect * "_quad_withF_noobs_L$L" * ".txt"; comments=true)
               calc_xis_quad = [GaPSE.ξ_GNC_multipole(COSMO.s_eff, s, name_effect, COSMO;
                    L=L, alg=:quad, kwargs...) for s in ss_quad]

               @test all([isapprox(xi, calc_xi, rtol=RTOL) for (xi, calc_xi) in zip(xis_quad, calc_xis_quad)])
          end

          @testset "L = 1" begin
               L = 1
               ss_lob, xis_lob = GaPSE.readxy("datatest/GNC_doppler_multipoles/xi_GNC_" *
                                              name_effect * "_lobatto_withF_noobs_L$L" * ".txt"; comments=true)
               ss_trap, xis_trap = GaPSE.readxy("datatest/GNC_doppler_multipoles/xi_GNC_" *
                                                name_effect * "_trap_withF_noobs_L$L" * ".txt"; comments=true)

               calc_xis_lob = [GaPSE.ξ_GNC_multipole(COSMO.s_eff, s, name_effect, COSMO;
                    L=L, alg=:lobatto, kwargs...) for s in ss_lob]
               calc_xis_trap = [GaPSE.ξ_GNC_multipole(COSMO.s_eff, s, name_effect, COSMO;
                    L=L, alg=:trap, kwargs...) for s in ss_trap]

               @test all([isapprox(xi, calc_xi, rtol=RTOL) for (xi, calc_xi) in zip(xis_lob, calc_xis_lob)])
               @test all([isapprox(xi, calc_xi, rtol=RTOL) for (xi, calc_xi) in zip(xis_trap, calc_xis_trap)])
          end

          @testset "L = 2" begin
               L = 2
               ss_lob, xis_lob = GaPSE.readxy("datatest/GNC_doppler_multipoles/xi_GNC_" *
                                              name_effect * "_lobatto_withF_noobs_L$L" * ".txt"; comments=true)
               ss_trap, xis_trap = GaPSE.readxy("datatest/GNC_doppler_multipoles/xi_GNC_" *
                                                name_effect * "_trap_withF_noobs_L$L" * ".txt"; comments=true)

               calc_xis_lob = [GaPSE.ξ_GNC_multipole(COSMO.s_eff, s, name_effect, COSMO;
                    L=L, alg=:lobatto, kwargs...) for s in ss_lob]
               calc_xis_trap = [GaPSE.ξ_GNC_multipole(COSMO.s_eff, s, name_effect, COSMO;
                    L=L, alg=:trap, kwargs...) for s in ss_trap]

               @test all([isapprox(xi, calc_xi, rtol=RTOL) for (xi, calc_xi) in zip(xis_lob, calc_xis_lob)])
               @test all([isapprox(xi, calc_xi, rtol=RTOL) for (xi, calc_xi) in zip(xis_trap, calc_xis_trap)])

               ss_quad, xis_quad = GaPSE.readxy("datatest/GNC_doppler_multipoles/xi_GNC_" *
                                                name_effect * "_quad_withF_noobs_L$L" * ".txt"; comments=true)
               calc_xis_quad = [GaPSE.ξ_GNC_multipole(COSMO.s_eff, s, name_effect, COSMO;
                    L=L, alg=:quad, kwargs...) for s in ss_quad]

               @test all([isapprox(xi, calc_xi, rtol=RTOL) for (xi, calc_xi) in zip(xis_quad, calc_xis_quad)])
               #println("xis_trap = $xis_trap;")
               #println("calc_xis_trap = $calc_xis_trap ;")
          end

          @testset "L = 3" begin
               L = 3
               ss_lob, xis_lob = GaPSE.readxy("datatest/GNC_doppler_multipoles/xi_GNC_" *
                                              name_effect * "_lobatto_withF_noobs_L$L" * ".txt"; comments=true)
               ss_trap, xis_trap = GaPSE.readxy("datatest/GNC_doppler_multipoles/xi_GNC_" *
                                                name_effect * "_trap_withF_noobs_L$L" * ".txt"; comments=true)

               calc_xis_lob = [GaPSE.ξ_GNC_multipole(COSMO.s_eff, s, name_effect, COSMO;
                    L=L, alg=:lobatto, kwargs...) for s in ss_lob]
               calc_xis_trap = [GaPSE.ξ_GNC_multipole(COSMO.s_eff, s, name_effect, COSMO;
                    L=L, alg=:trap, kwargs...) for s in ss_trap]

               @test all([isapprox(xi, calc_xi, rtol=RTOL) for (xi, calc_xi) in zip(xis_lob, calc_xis_lob)])
               @test all([isapprox(xi, calc_xi, rtol=RTOL) for (xi, calc_xi) in zip(xis_trap, calc_xis_trap)])
          end

          @testset "L = 4" begin
               L = 4
               ss_lob, xis_lob = GaPSE.readxy("datatest/GNC_doppler_multipoles/xi_GNC_" *
                                              name_effect * "_lobatto_withF_noobs_L$L" * ".txt"; comments=true)
               ss_trap, xis_trap = GaPSE.readxy("datatest/GNC_doppler_multipoles/xi_GNC_" *
                                                name_effect * "_trap_withF_noobs_L$L" * ".txt"; comments=true)

               calc_xis_lob = [GaPSE.ξ_GNC_multipole(COSMO.s_eff, s, name_effect, COSMO;
                    L=L, alg=:lobatto, kwargs...) for s in ss_lob]
               calc_xis_trap = [GaPSE.ξ_GNC_multipole(COSMO.s_eff, s, name_effect, COSMO;
                    L=L, alg=:trap, kwargs...) for s in ss_trap]

               @test all([isapprox(xi, calc_xi, rtol=RTOL) for (xi, calc_xi) in zip(xis_lob, calc_xis_lob)])
               @test all([isapprox(xi, calc_xi, rtol=RTOL) for (xi, calc_xi) in zip(xis_trap, calc_xis_trap)])
          end
     end
end

println("\nJust finished the tests on ξ_GNC_multipole.")


##########################################################################################92

println("Now I work on map_ξ_GNC_multipole...")

@testset "test print_map_ξ_GNC_multipole" begin
     effect = "auto_doppler"

     kwargs = Dict(
          :use_windows => false,
          :pr => false, :alg => :quad,
          :enhancer => 1e8, :N_trap => 30, :N_lob => 30,
          :atol_quad => 0.0, :rtol_quad => 1e-2,
          :N_log => 1000
     )

     @testset "zero" begin
          @test_throws AssertionError GaPSE.print_map_ξ_GNC_multipole(COSMO, "/Users/matteofoglieni/nonexistingdir/file.txt", effect)
          @test_throws AssertionError GaPSE.print_map_ξ_GNC_multipole(COSMO, "/nonexistingdir/here.txt", effect)
          @test_throws AssertionError GaPSE.print_map_ξ_GNC_multipole(COSMO, "/nonexistingdir/", effect)
          @test_throws AssertionError GaPSE.print_map_ξ_GNC_multipole(COSMO, "file", effect)
          @test_throws AssertionError GaPSE.print_map_ξ_GNC_multipole(COSMO, "file.boh", effect)
     end

     @testset "first" begin
          ss, xis = GaPSE.readxy("datatest/GNC_doppler_multipoles/xi_GNC_auto_doppler_L0_first.txt"; comments=true)

          name = "calc_xi_GNC_auto_doppler_L0_first.txt"
          isfile(name) && rm(name)
          GaPSE.print_map_ξ_GNC_multipole(COSMO, name, effect, nothing;
               s1=nothing, L=0, kwargs...)

          calc_table = readdlm(name; comments=true)
          calc_ss, calc_xis = GaPSE.readxy(name; comments=true)

          @test all([isapprox(s, calc_s, rtol=1e-2) for (s, calc_s) in zip(ss, calc_ss)])
          @test all([isapprox(xi, calc_xi, rtol=1e-2) for (xi, calc_xi) in zip(xis, calc_xis)])

          rm(name)
     end

     @testset "second" begin
          ss, xis = GaPSE.readxy("datatest/GNC_doppler_multipoles/xi_GNC_auto_doppler_L0_second.txt"; comments=true)

          name = "calc_xi_GNC_auto_doppler_L0_second.txt"
          isfile(name) && rm(name)
          GaPSE.print_map_ξ_GNC_multipole(COSMO, name, effect, 10 .^ range(0, 3, length=344);
               s1=nothing, L=0, kwargs...)

          calc_table = readdlm(name; comments=true)
          calc_ss, calc_xis = GaPSE.readxy(name; comments=true)

          @test all([isapprox(s, calc_s, rtol=1e-2) for (s, calc_s) in zip(ss, calc_ss)])
          @test all([isapprox(xi, calc_xi, rtol=1e-2) for (xi, calc_xi) in zip(xis, calc_xis)])

          rm(name)
     end

     @testset "third" begin
          ss, xis = GaPSE.readxy("datatest/GNC_doppler_multipoles/xi_GNC_auto_doppler_L0_third.txt"; comments=true)

          name = "calc_xi_GNC_auto_doppler_L0_third.txt"
          isfile(name) && rm(name)
          GaPSE.print_map_ξ_GNC_multipole(COSMO, name, effect, 10 .^ range(0, 3, length=344);
               s1=COSMO.s_eff - 65.0, L=0, kwargs...)

          calc_table = readdlm(name; comments=true)
          calc_ss, calc_xis = GaPSE.readxy(name; comments=true)

          @test all([isapprox(s, calc_s, rtol=1e-2) for (s, calc_s) in zip(ss, calc_ss)])
          @test all([isapprox(xi, calc_xi, rtol=1e-2) for (xi, calc_xi) in zip(xis, calc_xis)])

          rm(name)
     end
end

println("\ndone!")