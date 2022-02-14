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

@testset "derivate_point" begin
     @testset "first" begin
          f(x) = 1.35 + 2.54 * x
          xp, x1, x2 = 2.0, 1.0, 3.0
          yp, y1, y2 = f(xp), f(x1), f(x2)

          @test GaPSE.derivate_point(xp, yp, x1, y1, x2, y2) ≈ 2.54
     end

     @testset "second" begin
          f(x) = 2.634 + 4.65 * x^3.0
          xp, x1, x2 = 3.0, 3.0 - 1e-8, 3.0 + 1e-8
          yp, y1, y2 = f(xp), f(x1), f(x2)
          @test isapprox(GaPSE.derivate_point(xp, yp, x1, y1, x2, y2), 125.55, atol = 1e-4)
     end
end



@testset "mean_spectral_index" begin
     @testset "first" begin
          xs = range(1.0, 10.0, 100)
          ys = [3.54 * x^1.856 for x in xs]

          @test isapprox(GaPSE.mean_spectral_index(xs, ys; con = true), 1.856; atol = 1e-4)
          @test isapprox(GaPSE.mean_spectral_index(xs, ys; con = false), 1.856; atol = 1e-4)
     end

     @testset "second" begin
          xs = range(1.0, 10.0, 100)
          ys = [256.75 - 3.54 * x^1.856 for x in xs]

          @test isapprox(GaPSE.mean_spectral_index(xs, ys; con = true), 1.856; atol = 1e-4)
     end

end

@testset "test func_I04_tilde" begin
     table_ips = readdlm(FILE_PS)
     ks = convert(Vector{Float64}, table_ips[:, 1])
     pks = convert(Vector{Float64}, table_ips[:, 2])
     PK = Spline1D(ks, pks)
     kmin, kmax = 1e-5, 1e3

     ss = [3.0387609674820664]
     res = [-1.7693316431457908]
     my_res = [GaPSE.func_I04_tilde(PK, s, kmin, kmax) for s in ss]

     @test all([isapprox(my, tr, rtol = 1e-3) for (my, tr) in zip(my_res, res)])
end


#=
@testset "test expanded_I04_tilde" begin
     tab_I04_tildes = readdlm("datatest/I40_tilde.txt")
     ss = convert(Vector{Float64}, tab_I04_tildes[2:end, 1])
     I04_tildes = convert(Vector{Float64}, tab_I04_tildes[2:end, 2])

     table_ips = readdlm(FILE_PS)
     ks = convert(Vector{Float64}, table_ips[:, 1])
     pks = convert(Vector{Float64}, table_ips[:, 2])
     PK = Spline1D(ks, pks)
     kmin, kmax = 1e-5, 1e3
     calc_I04_tildes = GaPSE.expanded_I04_tilde(PK, ss; kmin = kmin, kmax = kmax)

     @test all([isapprox(my, tr, rtol = 1e-3) for (my, tr) in zip(calc_I04_tildes, I04_tildes)])
end
=#
