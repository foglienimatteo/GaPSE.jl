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


@testset "test_Cosmology" begin
     z_min, z_max, θ_max = 0.05, 0.20, π / 2.0
     params = GaPSE.CosmoParams(z_min, z_max, θ_max)
     cosmo = GaPSE.Cosmology(params, FILE_BACKGROUND,
          FILE_PS, FILE_F_MAP)

     RTOL = 1e-6

     @test isapprox(cosmo.params.Ω_b , 0.0489, rtol = RTOL)
     @test isapprox(cosmo.params.Ω_cdm ,  0.251020, rtol = RTOL)
     @test isapprox(cosmo.params.Ω_M0 , 0.29992, rtol = RTOL)
     @test isapprox(cosmo.params.h_0, 0.7, rtol = RTOL)

     for k in keys(GaPSE.DEFAULT_IPS_OPTS)
          @test cosmo.params.IPS[k] ≈ GaPSE.DEFAULT_IPS_OPTS[k]
     end
     for k in keys(GaPSE.DEFAULT_IPSTOOLS_OPTS)
          @test cosmo.params.IPSTools[k] ≈ GaPSE.DEFAULT_IPSTOOLS_OPTS[k]
     end


     @test isapprox(cosmo.z_eff, Z_EFF, rtol = RTOL)
     @test isapprox(cosmo.s_min, S_MIN, rtol = RTOL)
     @test isapprox(cosmo.s_max, S_MAX, rtol = RTOL)
     @test isapprox(cosmo.s_eff, S_EFF, rtol = RTOL)

     @test isapprox(cosmo.volume, VOLUME, rtol = RTOL)

     @test cosmo.file_data == FILE_BACKGROUND
     @test cosmo.file_ips == FILE_PS
     @test cosmo.file_windowF == FILE_F_MAP

     @test all(isapprox(cosmo.z_of_s.(COM_DIST), ZS, rtol = RTOL))
     @test all(isapprox(cosmo.s_of_z.(ZS), COM_DIST, rtol = RTOL))
     @test all(isapprox(cosmo.D_of_s.(COM_DIST), GROWTH_FACTOR_D, rtol = RTOL))
     @test all(isapprox(cosmo.f_of_s.(COM_DIST), GROWTH_FACTOR_F, rtol = RTOL))
     @test all(isapprox(cosmo.ℋ_of_s.(COM_DIST), COM_H, rtol = RTOL))
end


