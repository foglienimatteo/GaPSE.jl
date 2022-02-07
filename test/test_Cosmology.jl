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
     z_min = 0.05
     z_max = 0.20
     θ_max = π / 2.0
     params = GaPSE.CosmoParams(z_min, z_max, θ_max;
          k_min = 1e-8, k_max = 10.0,
          Ω_b = 0.0489, Ω_cdm = 0.251020, h_0 = 0.70)
     cosmo = GaPSE.Cosmology(params,
          FILE_BACKGROUND,
          FILE_PS,
          FILE_F_MAP)

     @test all(isapprox(cosmo.z_of_s.(COM_DIST), ZS, rtol = 1e-8))
     @test all(isapprox(cosmo.s_of_z.(ZS), COM_DIST, rtol = 1e-8))
     @test all(isapprox(cosmo.D_of_s.(COM_DIST), GROWTH_FACTOR_D, rtol = 1e-8))
     @test all(isapprox(cosmo.f_of_s.(COM_DIST), GROWTH_FACTOR_F, rtol = 1e-8))
     @test all(isapprox(cosmo.ℋ_of_s.(COM_DIST), COM_H, rtol = 1e-8))
end
