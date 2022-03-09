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



@testset "test constant values" begin
     @test isapprox(GaPSE.ℋ0, HUBBLE_0; rtol = 1e-6)
     z_of_s = Spline1D(COM_DIST, ZS)

     @test isapprox(GaPSE.func_z_eff(S_MAX, S_MIN, z_of_s), Z_EFF; rtol = 1e-6)
end

@testset "test BackgroundData" begin
     BD = GaPSE.BackgroundData(FILE_BACKGROUND, 0.2)
     L = length(CONF_TIME)

     @test all(isapprox(BD.z[begin:L], ZS, rtol = 1e-8))
     @test all(isapprox.(BD.conftime[begin:L], CONF_TIME, rtol = 1e-8))
     @test all(isapprox.(BD.comdist[begin:L], COM_DIST, rtol = 1e-8))
     @test all(isapprox.(BD.angdist[begin:L], ANG_DIST, rtol = 1e-8))
     @test all(isapprox.(BD.lumdist[begin:L], LUM_DIST, rtol = 1e-8))
     @test all(isapprox.(BD.D[begin:L], GROWTH_FACTOR_D, rtol = 1e-8))
     @test all(isapprox.(BD.f[begin:L], GROWTH_FACTOR_F, rtol = 1e-8))
     @test all(isapprox.(BD.ℋ[begin:L], COM_H, rtol = 1e-8))
     @test all(isapprox.(BD.ℋ_p[3:L], COM_H_P[3:end], rtol = 1e-4))
end


@testset "test CosmoParams" begin

     @testset "zeros" begin
          @test_throws AssertionError GaPSE.CosmoParams(0.0, 1.0, π/2.0) 
          @test_throws AssertionError GaPSE.CosmoParams(-1.0, 1.0, π/2.0) 
          @test_throws AssertionError GaPSE.CosmoParams(2.0, 1.0, π/2.0) 
          @test_throws AssertionError GaPSE.CosmoParams(0.5, 1.0, 1.5*π) 

          @test_throws AssertionError GaPSE.CosmoParams(0.5, 1.0, π/2.0; Ω_b = -0.2)
          @test_throws AssertionError GaPSE.CosmoParams(0.5, 1.0, π/2.0; Ω_b = 13) 
          @test_throws AssertionError GaPSE.CosmoParams(0.5, 1.0, π/2.0; Ω_cdm = -0.2)
          @test_throws AssertionError GaPSE.CosmoParams(0.5, 1.0, π/2.0; Ω_cdm = 243)
          @test_throws AssertionError GaPSE.CosmoParams(0.5, 1.0, π/2.0; h_0 = 0.0)
          @test_throws AssertionError GaPSE.CosmoParams(0.5, 1.0, π/2.0; h_0 = 1.5)

          @test_throws AssertionError GaPSE.CosmoParams(0.5, 1.0, π/2.0; 
               IPS_opts = Dict())
          @test_throws AssertionError GaPSE.CosmoParams(0.5, 1.0, π/2.0; 
               IPSTools_opts = Dict())

          @test_throws AssertionError GaPSE.CosmoParams(0.5, 1.0, π/2.0; 
               IPS_opts = Dict(:k_min => true))
          @test_throws AssertionError GaPSE.CosmoParams(0.5, 1.0, π/2.0; 
               IPS_opts = Dict(:N => 12))
          @test_throws AssertionError GaPSE.CosmoParams(0.5, 1.0, π/2.0; 
               IPS_opts = Dict("k_min" => 12))

          @test_throws AssertionError GaPSE.CosmoParams(0.5, 1.0, π/2.0; 
               IPS_opts = Dict(:N => 12.3))
          @test_throws AssertionError GaPSE.CosmoParams(0.5, 1.0, π/2.0; 
               IPS_opts = Dict(:M => 12.3))
          @test_throws AssertionError GaPSE.CosmoParams(0.5, 1.0, π/2.0; 
               IPS_opts = Dict("N" => 12))
     end

     @testset "first" begin
          z_min, z_max, θ_max = 0.05, 0.20, π/2.0
          Ω_b, Ω_cdm, h_0 = 0.023, 0.34, 0.99
          s_lim = 1e-3

          params = GaPSE.CosmoParams(z_min, z_max, θ_max;
               Ω_b = Ω_b, Ω_cdm = Ω_cdm, h_0 = h_0, s_lim = s_lim,
               IPS_opts = Dict{Symbol, Any}(),
               IPSTools_opts = Dict{Symbol, Any}(),
               ) 

          @test params.h_0 ≈ h_0
          @test params.Ω_b ≈ Ω_b
          @test params.Ω_cdm ≈ Ω_cdm
          @test params.Ω_M0 ≈ Ω_b + Ω_cdm
          @test params.s_lim ≈ s_lim

          for k in keys(GaPSE.DEFAULT_IPS_OPTS)
               @test params.IPS[k] ≈ GaPSE.DEFAULT_IPS_OPTS[k]
          end
          for k in keys(GaPSE.DEFAULT_IPSTOOLS_OPTS)
               @test params.IPSTools[k] ≈ GaPSE.DEFAULT_IPSTOOLS_OPTS[k]
          end
     end

     @testset "second" begin
          z_min, z_max, θ_max = 0.05, 0.20, π/2.0
          Ω_b, Ω_cdm, h_0 = 0.023, 0.34, 0.99
          s_lim = 1e-3

          A = Dict(:fit_left_min => 1e-20, :fit_right_min => 0.7)
          B = Dict(:N => 12, :con => false)

          params = GaPSE.CosmoParams(z_min, z_max, θ_max;
               Ω_b = Ω_b, Ω_cdm = Ω_cdm, h_0 = h_0, s_lim = s_lim,
               IPS_opts = A,
               IPSTools_opts = B,
               ) 

          @test params.h_0 ≈ h_0
          @test params.Ω_b ≈ Ω_b
          @test params.Ω_cdm ≈ Ω_cdm
          @test params.Ω_M0 ≈ Ω_b + Ω_cdm
          @test params.s_lim ≈ s_lim

          for k in keys(A)
               @test params.IPS[k] ≈ A[k]
          end
          for k in keys(B)
               @test params.IPSTools[k] ≈ B[k]
          end
          for k in filter(x -> x ∉ keys(A), keys(GaPSE.DEFAULT_IPS_OPTS))
               @test params.IPS[k] ≈ GaPSE.DEFAULT_IPS_OPTS[k]
          end
          for k in filter(x -> x ∉ keys(B), keys(GaPSE.DEFAULT_IPSTOOLS_OPTS))
               @test params.IPSTools[k] ≈ GaPSE.DEFAULT_IPSTOOLS_OPTS[k]
          end

     end

end
