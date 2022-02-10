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

@testset "test struct BackgroundData" begin
     BD = GaPSE.BackgroundData(FILE_BACKGROUND, 0.2)

     @test all(isapprox(BD.z, ZS, rtol = 1e-8))
     @test all(isapprox.(BD.conftime, CONF_TIME, rtol = 1e-8))
     @test all(isapprox.(BD.comdist, COM_DIST, rtol = 1e-8))
     @test all(isapprox.(BD.angdist, ANG_DIST, rtol = 1e-8))
     @test all(isapprox.(BD.lumdist, LUM_DIST, rtol = 1e-8))
     @test all(isapprox.(BD.D, GROWTH_FACTOR_D, rtol = 1e-8))
     @test all(isapprox.(BD.f, GROWTH_FACTOR_F, rtol = 1e-8))
     @test all(isapprox.(BD.ℋ, COM_H, rtol = 1e-8))
     @test all(isapprox.(BD.ℋ_p[3:end], COM_H_P[3:end], rtol = 1e-4))
end



