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
     @test isapprox(GaPSE.ℋ0, 1e5 / 299792458.0; rtol = 1e-6)

     s_min = 148.1920001343431
     s_max = 571.7022420911966
     s_eff = 435.37470960794167
     z_eff = 0.15045636097417317

     z_of_s = Spline1D(com_dist, zs)

     @test isapprox(GaPSE.func_z_eff(s_max, s_min, z_of_s), z_eff; rtol = 1e-6)
end

@testset "test struct Background_functions" begin
     BF = GaPSE.BackgroundData(TEST_FILE, 0.05, 0.2)

     #println([isapprox(BF_z, z, rtol = 1e-8) for (BF_z, z) in zip(BF.zs, zs)])

     @test all(isapprox(BF.z, zs, rtol = 1e-8))
     @test all(isapprox.(BF.conftime, conf_time, rtol = 1e-8))
     @test all(isapprox.(BF.comdist, com_dist, rtol = 1e-8))
     @test all(isapprox.(BF.angdist, ang_dist, rtol = 1e-8))
     @test all(isapprox.(BF.lumdist, lum_dist, rtol = 1e-8))
     @test all(isapprox.(BF.D, growth_factor_D, rtol = 1e-8))
     @test all(isapprox.(BF.f, growth_factor_f, rtol = 1e-8))
     @test all(isapprox.(BF.ℋ, com_H, rtol = 1e-8))
     @test all(isapprox.(BF.ℋ_p, com_H_p, rtol = 1e-8))
end
