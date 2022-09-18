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


@testset "TF" begin
     RTOL = 1e-2
     ks, pks = GaPSE.readxy("datatest/Tk.dat")
     my_TF = GaPSE.TF(ks, pks)

     @test my_TF(0.6e-1) ≈ 0.22902775397420208

     ks_true, pks_true = GaPSE.readxy("datatest/PNG/Transfer_Function.dat")
     pks_calc = [my_TF(k) for k in ks_true]

     @test all([isapprox(a, b; rtol=RTOL) for (a, b) in zip(pks_true, pks_calc)])
end


@testset "alpha_bias" begin
     RTOL = 1e-3
     ks, pks = GaPSE.readxy("datatest/Tk.dat")
     my_TF = GaPSE.TF(ks, pks)

     @test isapprox(GaPSE.α_bias(1, my_TF;
               bf=1.0, D=1.0, Ω_M0=0.29992), 1.3942430956475444e-5; rtol=RTOL)
     @test isapprox(GaPSE.α_bias(1, my_TF;
               bf=5.0, D=1.0, Ω_M0=0.29992), 6.971215478237722e-5; rtol=RTOL)
     @test isapprox(GaPSE.α_bias(1, my_TF;
               bf=5.0, D=0.13, Ω_M0=0.29992), 0.0005362473444798248; rtol=RTOL)
     @test isapprox(GaPSE.α_bias(1, my_TF;
               bf=5.0, D=0.13, Ω_M0=0.9), 0.001609171145745006; rtol=RTOL)
     @test isapprox(GaPSE.α_bias(7.934, my_TF;
               bf=5.0, D=0.13, Ω_M0=0.9), 0.0008966355521073972; rtol=RTOL)

     xs_true, ys_true = GaPSE.readxy("datatest/PNG/alpha_bias.dat")
     ys_calc = [GaPSE.α_bias(x, my_TF;
          bf=1.0, D=1.0, Ω_M0=0.29992) for x in xs_true]

     @test all([isapprox(a, b; rtol=RTOL) for (a, b) in zip(ys_true, ys_calc)])
end

@testset "IntegralIPSalpha" begin
     RTOL = 1e-3
     ks, pks = GaPSE.readxy("datatest/Tk.dat")
     my_TF = GaPSE.TF(ks, pks)
     my_IntIPSalpha = GaPSE.IntegralIPSalpha(my_TF, COSMO, 0, 0;
          D=nothing, bf=1.0, N=1024, kmin=1e-6, kmax=1e4, s0=1e-4,
          fit_left_min=nothing, fit_left_max=nothing, p0_left=nothing,
          fit_right_min=nothing, fit_right_max=nothing, p0_right=nothing)

     @test isapprox(my_IntIPSalpha(1.34), 9.484203885952665e-5; rtol=RTOL)

     xs_true, ys_true = GaPSE.readxy("datatest/PNG/Integral_IPS_alpha.dat")
     ys_calc = [my_IntIPSalpha(x) for x in xs_true]

     @test all([isapprox(a, b; rtol=RTOL) for (a, b) in zip(ys_true, ys_calc)])
end

