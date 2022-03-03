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


@testset "test ξ_multipole first method" begin
     effect = "auto_doppler"
     RTOL = 1e-2

     kwargs = Dict(
          :use_windows => false,
          :enhancer => 1e8, :N_μs => 30,
          :μ_atol => 0.0, :μ_rtol => 1e-2,
     )

     #=
     @testset "zeros" begin
          @test_throws AssertionError GaPSE.ξ_multipole(COSMO.s_eff, 10.0, "strange", COSMO;
               L = L, kwargs...)
     end
     =#

     @testset "L = 0" begin
          L = 0
          table = readdlm("datatest/doppler_multipoles/xi_" * effect * "_L$L" * ".txt"; comments = true)
          ss = convert(Vector{Float64}, table[:, 1])
          xis = convert(Vector{Float64}, table[:, 2])

          calc_xis = [GaPSE.ξ_multipole(COSMO.s_eff, s, effect, COSMO;
               L = L, kwargs...) for s in ss]

          @test all([isapprox(xi, calc_xi, rtol = RTOL) for (xi, calc_xi) in zip(xis, calc_xis)])
     end

     @testset "L = 1" begin
          L = 1
          table = readdlm("datatest/doppler_multipoles/xi_" * effect * "_L$L" * ".txt"; comments = true)
          ss = convert(Vector{Float64}, table[:, 1])
          xis = convert(Vector{Float64}, table[:, 2])

          calc_xis = [GaPSE.ξ_multipole(COSMO.s_eff, s, effect, COSMO;
               L = L, kwargs...) for s in ss]

          @test all([isapprox(xi, calc_xi, rtol = RTOL) for (xi, calc_xi) in zip(xis, calc_xis)])
     end

     @testset "L = 2" begin
          L = 2
          table = readdlm("datatest/doppler_multipoles/xi_" * effect * "_L$L" * ".txt"; comments = true)
          ss = convert(Vector{Float64}, table[:, 1])
          xis = convert(Vector{Float64}, table[:, 2])

          calc_xis = [GaPSE.ξ_multipole(COSMO.s_eff, s, effect, COSMO;
               L = L, kwargs...) for s in ss]

          @test all([isapprox(xi, calc_xi, rtol = RTOL) for (xi, calc_xi) in zip(xis, calc_xis)])
     end

     @testset "L = 3" begin
          L = 3
          table = readdlm("datatest/doppler_multipoles/xi_" * effect * "_L$L" * ".txt"; comments = true)
          ss = convert(Vector{Float64}, table[:, 1])
          xis = convert(Vector{Float64}, table[:, 2])

          calc_xis = [GaPSE.ξ_multipole(COSMO.s_eff, s, effect, COSMO;
               L = L, kwargs...) for s in ss]

          @test all([isapprox(xi, calc_xi, rtol = RTOL) for (xi, calc_xi) in zip(xis, calc_xis)])
     end

     @testset "L = 4" begin
          L = 4
          table = readdlm("datatest/doppler_multipoles/xi_" * effect * "_L$L" * ".txt"; comments = true)
          ss = convert(Vector{Float64}, table[:, 1])
          xis = convert(Vector{Float64}, table[:, 2])

          calc_xis = [GaPSE.ξ_multipole(COSMO.s_eff, s, effect, COSMO;
               L = L, kwargs...) for s in ss]

          @test all([isapprox(xi, calc_xi, rtol = RTOL) for (xi, calc_xi) in zip(xis, calc_xis)])
     end
end


@testset "test ξ_multipole second method" begin
     effect = GaPSE.ξ_Doppler
     name_effect = "auto_doppler"
     RTOL = 1e-2

     kwargs = Dict(
          :use_windows => false,
          :enhancer => 1e8, :N_μs => 30,
          :μ_atol => 0.0, :μ_rtol => 1e-2,
     )

     #=
     @testset "zeros" begin
          function a_func(x)
               return x
          end

          @test_throws AssertionError GaPSE.ξ_multipole(COSMO.s_eff, 10.0, a_func, COSMO;
               L = L, kwargs...)
     end
     =#

     @testset "L = 0" begin
          L = 0
          table = readdlm("datatest/doppler_multipoles/xi_" * name_effect * "_L$L" * ".txt"; comments = true)
          ss = convert(Vector{Float64}, table[:, 1])
          xis = convert(Vector{Float64}, table[:, 2])

          calc_xis = [GaPSE.ξ_multipole(COSMO.s_eff, s, effect, COSMO;
               L = L, kwargs...) for s in ss]

          @test all([isapprox(xi, calc_xi, rtol = RTOL) for (xi, calc_xi) in zip(xis, calc_xis)])
     end

     @testset "L = 1" begin
          L = 1
          table = readdlm("datatest/doppler_multipoles/xi_" * name_effect * "_L$L" * ".txt"; comments = true)
          ss = convert(Vector{Float64}, table[:, 1])
          xis = convert(Vector{Float64}, table[:, 2])

          calc_xis = [GaPSE.ξ_multipole(COSMO.s_eff, s, effect, COSMO;
               L = L, kwargs...) for s in ss]

          @test all([isapprox(xi, calc_xi, rtol = RTOL) for (xi, calc_xi) in zip(xis, calc_xis)])
     end

     @testset "L = 2" begin
          L = 2
          table = readdlm("datatest/doppler_multipoles/xi_" * name_effect * "_L$L" * ".txt"; comments = true)
          ss = convert(Vector{Float64}, table[:, 1])
          xis = convert(Vector{Float64}, table[:, 2])

          calc_xis = [GaPSE.ξ_multipole(COSMO.s_eff, s, effect, COSMO;
               L = L, kwargs...) for s in ss]

          @test all([isapprox(xi, calc_xi, rtol = RTOL) for (xi, calc_xi) in zip(xis, calc_xis)])
     end

     @testset "L = 3" begin
          L = 3
          table = readdlm("datatest/doppler_multipoles/xi_" * name_effect * "_L$L" * ".txt"; comments = true)
          ss = convert(Vector{Float64}, table[:, 1])
          xis = convert(Vector{Float64}, table[:, 2])

          calc_xis = [GaPSE.ξ_multipole(COSMO.s_eff, s, effect, COSMO;
               L = L, kwargs...) for s in ss]

          @test all([isapprox(xi, calc_xi, rtol = RTOL) for (xi, calc_xi) in zip(xis, calc_xis)])
     end

     @testset "L = 4" begin
          L = 4
          table = readdlm("datatest/doppler_multipoles/xi_" * name_effect * "_L$L" * ".txt"; comments = true)
          ss = convert(Vector{Float64}, table[:, 1])
          xis = convert(Vector{Float64}, table[:, 2])

          calc_xis = [GaPSE.ξ_multipole(COSMO.s_eff, s, effect, COSMO;
               L = L, kwargs...) for s in ss]

          @test all([isapprox(xi, calc_xi, rtol = RTOL) for (xi, calc_xi) in zip(xis, calc_xis)])
     end
end
