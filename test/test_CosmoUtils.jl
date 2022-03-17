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

@testset "test func_z_eff" begin
     z_of_s = Spline1D(COM_DIST, ZS)
     @test isapprox(GaPSE.func_z_eff(S_MAX, S_MIN, z_of_s), Z_EFF; rtol = 1e-6)
end


@testset "V_survey" begin
     @test isapprox(GaPSE.V_survey(1, 2, π / 2), 14.660765716752364, rtol = 1e-4)
     @test isapprox(GaPSE.V_survey(10, 20, π / 2), 14660.765716752363, rtol = 1e-4)

     @test isapprox(GaPSE.V_survey(1, 2, 0.0), 0.0, rtol = 1e-4)
     @test isapprox(GaPSE.V_survey(1, 2, π / 6), 1.9641701671128426, rtol = 1e-4)

     @test isapprox(GaPSE.V_survey(S_MIN, S_MAX, π / 2), 3.845366169354268e8, rtol = 1e-4)
end

@testset "test ϕ" begin
     @test GaPSE.ϕ(0, 1.0, 2.0) ≈ 0.0
     @test GaPSE.ϕ(0.5, 1.0, 2.0) ≈ 0.0
     @test GaPSE.ϕ(1.0, 1.0, 2.0) ≈ 0.0
     @test GaPSE.ϕ(1.0 + 1e-5, 1.0, 2.0) ≈ 1.0
     @test GaPSE.ϕ(1.5, 1.0, 2.0) ≈ 1.0
     @test GaPSE.ϕ(2.0 - 1e-5, 1.0, 2.0) ≈ 1.0
     @test GaPSE.ϕ(2.0, 1.0, 2.0) ≈ 0.0
     @test GaPSE.ϕ(2.5, 1.0, 2.0) ≈ 0.0
end

@testset "test W" begin
     @test GaPSE.W(0, π / 3) ≈ 1.0
     @test GaPSE.W(1.0, π / 3) ≈ 1.0
     @test GaPSE.W(π / 3, π / 3) ≈ 0.0
     @test GaPSE.W(π, π / 3) ≈ 0.0
end
