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



@testset "test ξ_LD_multipole" begin
     name_effect = "auto_doppler"
     func_effect = GaPSE.ξ_LD_Doppler
     RTOL = 1e-2

     kwargs =  Dict(
          :use_windows => false,
          :enhancer => 1e8, :N_μs => 200,
          :μ_atol => 0.0, :μ_rtol => 1e-2
     )


     @testset "zeros" begin
          @test_throws AssertionError GaPSE.ξ_LD_multipole(COSMO.s_eff, 10.0, "strange", COSMO;
               L = 0, kwargs...)
     end

     @testset "L = 0" begin
          L = 0
          table = readdlm("datatest/LD_doppler_multipoles/xi_" * name_effect * "_L$L" * ".txt"; comments = true)
          ss = convert(Vector{Float64}, table[:, 1])
          xis = convert(Vector{Float64}, table[:, 2])

          calc_xis_1 = [GaPSE.ξ_LD_multipole(COSMO.s_eff, s, name_effect, COSMO;
               L = L, kwargs...) for s in ss]
          calc_xis_2 = [GaPSE.ξ_LD_multipole(COSMO.s_eff, s, func_effect, COSMO;
               L = L, kwargs...) for s in ss]

          @test all([isapprox(xi, calc_xi, rtol = RTOL) for (xi, calc_xi) in zip(xis, calc_xis_1)])
          @test all([isapprox(xi, calc_xi, rtol = RTOL) for (xi, calc_xi) in zip(xis, calc_xis_2)])
     end

     @testset "L = 1" begin
          L = 1
          table = readdlm("datatest/LD_doppler_multipoles/xi_" * name_effect * "_L$L" * ".txt"; comments = true)
          ss = convert(Vector{Float64}, table[:, 1])
          xis = convert(Vector{Float64}, table[:, 2])

          calc_xis_1 = [GaPSE.ξ_LD_multipole(COSMO.s_eff, s, name_effect, COSMO;
               L = L, kwargs...) for s in ss]
          calc_xis_2 = [GaPSE.ξ_LD_multipole(COSMO.s_eff, s, func_effect, COSMO;
               L = L, kwargs...) for s in ss]

          @test all([isapprox(xi, calc_xi, rtol = RTOL) for (xi, calc_xi) in zip(xis, calc_xis_1)])
          @test all([isapprox(xi, calc_xi, rtol = RTOL) for (xi, calc_xi) in zip(xis, calc_xis_2)])
     end

     @testset "L = 2" begin
          L = 2
          table = readdlm("datatest/LD_doppler_multipoles/xi_" * name_effect * "_L$L" * ".txt"; comments = true)
          ss = convert(Vector{Float64}, table[:, 1])
          xis = convert(Vector{Float64}, table[:, 2])

          calc_xis_1 = [GaPSE.ξ_LD_multipole(COSMO.s_eff, s, name_effect, COSMO;
               L = L, kwargs...) for s in ss]
          calc_xis_2 = [GaPSE.ξ_LD_multipole(COSMO.s_eff, s, func_effect, COSMO;
               L = L, kwargs...) for s in ss]

          @test all([isapprox(xi, calc_xi, rtol = RTOL) for (xi, calc_xi) in zip(xis, calc_xis_1)])
          @test all([isapprox(xi, calc_xi, rtol = RTOL) for (xi, calc_xi) in zip(xis, calc_xis_2)])
     end

     @testset "L = 3" begin
          L = 3
          table = readdlm("datatest/LD_doppler_multipoles/xi_" * name_effect * "_L$L" * ".txt"; comments = true)
          ss = convert(Vector{Float64}, table[:, 1])
          xis = convert(Vector{Float64}, table[:, 2])

          calc_xis_1 = [GaPSE.ξ_LD_multipole(COSMO.s_eff, s, name_effect, COSMO;
               L = L, kwargs...) for s in ss]
          calc_xis_2 = [GaPSE.ξ_LD_multipole(COSMO.s_eff, s, func_effect, COSMO;
               L = L, kwargs...) for s in ss]

          @test all([isapprox(xi, calc_xi, rtol = RTOL) for (xi, calc_xi) in zip(xis, calc_xis_1)])
          @test all([isapprox(xi, calc_xi, rtol = RTOL) for (xi, calc_xi) in zip(xis, calc_xis_2)])
     end

     @testset "L = 4" begin
          L = 4
          table = readdlm("datatest/LD_doppler_multipoles/xi_" * name_effect * "_L$L" * ".txt"; comments = true)
          ss = convert(Vector{Float64}, table[:, 1])
          xis = convert(Vector{Float64}, table[:, 2])

          calc_xis_1 = [GaPSE.ξ_LD_multipole(COSMO.s_eff, s, name_effect, COSMO;
               L = L, kwargs...) for s in ss]
          calc_xis_2 = [GaPSE.ξ_LD_multipole(COSMO.s_eff, s, func_effect, COSMO;
               L = L, kwargs...) for s in ss]

          @test all([isapprox(xi, calc_xi, rtol = RTOL) for (xi, calc_xi) in zip(xis, calc_xis_1)])
          @test all([isapprox(xi, calc_xi, rtol = RTOL) for (xi, calc_xi) in zip(xis, calc_xis_2)])
     end
end



##########################################################################################92



@testset "test print_map_ξ_LD_multipole" begin
     effect = "auto_doppler"

     kwargs = Dict(
          :use_windows => false,
          :pr => false,
          :enhancer => 1e8, :N_μs => 30,
          :μ_atol => 0.0, :μ_rtol => 1e-2,
          :N_log => 1000
     )

     @testset "zero" begin
          @test_throws AssertionError GaPSE.print_map_ξ_LD_multipole(COSMO, "/Users/matteofoglieni/nonexistingdir/file.txt", effect)
          @test_throws AssertionError GaPSE.print_map_ξ_LD_multipole(COSMO, "/nonexistingdir/here.txt", effect)
          @test_throws AssertionError GaPSE.print_map_ξ_LD_multipole(COSMO, "/nonexistingdir/", effect)
          @test_throws AssertionError GaPSE.print_map_ξ_LD_multipole(COSMO, "file", effect)
          @test_throws AssertionError GaPSE.print_map_ξ_LD_multipole(COSMO, "file.boh", effect)
     end

     @testset "first" begin
          table = readdlm("datatest/LD_doppler_multipoles/xi_auto_doppler_L0_first.txt"; comments = true)
          ss = convert(Vector{Float64}, table[:, 1])
          xis = convert(Vector{Float64}, table[:, 2])

          name = "calc_xi_auto_doppler_L0_first.txt"
          isfile(name) && rm(name)
          GaPSE.print_map_ξ_LD_multipole(COSMO, name, effect, nothing;
               s1 = nothing, L = 0, kwargs...)

          calc_table = readdlm(name; comments = true)
          calc_ss = convert(Vector{Float64}, calc_table[:, 1])
          calc_xis = convert(Vector{Float64}, calc_table[:, 2])

          @test all([isapprox(s, calc_s, rtol = 1e-2) for (s, calc_s) in zip(ss, calc_ss)])
          @test all([isapprox(xi, calc_xi, rtol = 1e-2) for (xi, calc_xi) in zip(xis, calc_xis)])

          rm(name)
     end

     @testset "second" begin
          table = readdlm("datatest/LD_doppler_multipoles/xi_auto_doppler_L0_second.txt"; comments = true)
          ss = convert(Vector{Float64}, table[:, 1])
          xis = convert(Vector{Float64}, table[:, 2])

          name = "calc_xi_auto_doppler_L0_second.txt"
          isfile(name) && rm(name)
          GaPSE.print_map_ξ_LD_multipole(COSMO, name, effect, 10 .^ range(0, 3, length = 344);
               s1 = nothing, L = 0, kwargs...)

          calc_table = readdlm(name; comments = true)
          calc_ss = convert(Vector{Float64}, calc_table[:, 1])
          calc_xis = convert(Vector{Float64}, calc_table[:, 2])

          @test all([isapprox(s, calc_s, rtol = 1e-2) for (s, calc_s) in zip(ss, calc_ss)])
          @test all([isapprox(xi, calc_xi, rtol = 1e-2) for (xi, calc_xi) in zip(xis, calc_xis)])

          rm(name)
     end

     @testset "third" begin
          table = readdlm("datatest/LD_doppler_multipoles/xi_auto_doppler_L0_third.txt"; comments = true)
          ss = convert(Vector{Float64}, table[:, 1])
          xis = convert(Vector{Float64}, table[:, 2])

          name = "calc_xi_auto_doppler_L0_third.txt"
          isfile(name) && rm(name)
          GaPSE.print_map_ξ_LD_multipole(COSMO, name, effect, 10 .^ range(0, 3, length = 344);
               s1 = COSMO.s_eff - 65.0, L = 0, kwargs...)

          calc_table = readdlm(name; comments = true)
          calc_ss = convert(Vector{Float64}, calc_table[:, 1])
          calc_xis = convert(Vector{Float64}, calc_table[:, 2])

          @test all([isapprox(s, calc_s, rtol = 1e-2) for (s, calc_s) in zip(ss, calc_ss)])
          @test all([isapprox(xi, calc_xi, rtol = 1e-2) for (xi, calc_xi) in zip(xis, calc_xis)])

          rm(name)
     end
end
