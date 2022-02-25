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



@testset "test derivate_vector" begin
     @testset "zeros" begin
          xs = 1:0.1:2
          ys = [2.5 + 3.65*x^(-3.5) for x in xs]
          @test_throws AssertionError GaPSE.derivate_vector(xs, ys[begin+1:end])
          @test_throws AssertionError GaPSE.derivate_vector(xs[begin+1:end], ys)
          @test_throws AssertionError GaPSE.derivate_vector(xs, ys; N=6)
     end

     @testset "first" begin
          xs = 1:0.1:2
          ys = [2.5 + 3.65*x for x in xs]
          calc = GaPSE.derivate_vector(xs, ys; N=1)
          @test all([isapprox(3.65, c; rtol=1e-4) for c in calc])
     end

     @testset "second" begin
          xs = 1:0.01:2
          ys = [2.5 + 3.65*x^2.0 for x in xs]
          calc = GaPSE.derivate_vector(xs, ys)
          @test all([isapprox(1.0, 3.65*2.0*x/c; rtol=1e-4) 
               for (x,c) in zip(xs[begin+1:end-1], calc[begin+1:end-1])])
     end

     @testset "third" begin
          xs = 1:0.01:2
          ys = [2.5 + 3.65*x^3.5 for x in xs]
          calc = GaPSE.derivate_vector(xs, ys)
          @test all([isapprox(1.0, 3.65*3.5*x^2.5/c; rtol=1e-4) 
               for (x,c) in zip(xs[begin+1:end-1], calc[begin+1:end-1])])
     end

     @testset "fourth" begin
          xs = 1:0.01:2
          ys = [2.5 + 3.65*x^(-3.5) for x in xs]
          calc = GaPSE.derivate_vector(xs, ys)
          @test all([isapprox(1.0, -3.65*3.5*x^(-4.5)/c; rtol=1e-2) 
               for (x,c) in zip(xs[begin+1:end-1], calc[begin+1:end-1])])
     end

     @testset "fifth" begin
          xs = 10 .^ range(4,6, length=100)
          ys = [2.5e99 + 3.65*x^(3.5) for x in xs]
          @test_throws AssertionError GaPSE.derivate_vector(xs, ys; N=1)
          @test_throws AssertionError GaPSE.derivate_vector(xs, ys; N=6)
     end

     @testset "sixth" begin
          xs = 10 .^ range(4,6, length=100)
          ys = [2.5e23 + 3.65*x^(3.5) for x in xs]
          calc = GaPSE.derivate_vector(xs, ys)
          @test all([isapprox(1.0, 3.65*3.5*x^(2.5)/c; rtol=1e-2) 
               for (x,c) in zip(xs[begin+1:end-1], calc[begin+1:end-1])])
     end

     @testset "seventh" begin
          xs = 10 .^ range(-6,-4, length=100)
          ys = [3.6e-22 + 3.65*x^(3.5) for x in xs]
          calc = GaPSE.derivate_vector(xs, ys)
          @test all([isapprox(1.0, 3.65*3.5*x^(2.5)/c; rtol=1e-2) 
               for (x,c) in zip(xs[begin+1:end-1], calc[begin+1:end-1])])
     end

     @testset "eigth" begin
          xs = 10 .^ range(-6,-4, length=100)
          ys = [3.6e-22 + 3.65*x^(-3.5) for x in xs]
          calc = GaPSE.derivate_vector(xs, ys)
          @test all([isapprox(1.0, -3.65*3.5*x^(-4.5)/c; rtol=1e-2) 
               for (x,c) in zip(xs[begin+1:end-1], calc[begin+1:end-1])])
     end
end



@testset "test spectral_index" begin
     @testset "zeros" begin
          xs = 1:0.1:2
          ys = [2.5 + 3.65*x^(-3.5) for x in xs]
          @test_throws AssertionError GaPSE.spectral_index(xs, ys[begin+1:end])
          @test_throws AssertionError GaPSE.spectral_index(xs[begin+1:end], ys)
          @test_throws AssertionError GaPSE.spectral_index(xs, ys; N=6)
     end


     @testset "first" begin
          xs = 10:0.1:20
          ys = [3.65*x^(-1.5) for x in xs]
          calc = GaPSE.spectral_index(xs, ys; N=1, con=false)
          @test all([isapprox(-1.5, c; rtol=1e-3) for c in calc])
     end

     @testset "second" begin
          xs = 10:0.1:20
          ys = [1434.5 + 3.65*x^(-1.5) for x in xs]
          calc = GaPSE.spectral_index(xs, ys; N=1, con=true)
          @test all([isapprox(-1.5, c; rtol=1e-3) for c in calc])
     end

     @testset "third" begin
          xs = 10 .^ range(4,6, length=100)
          ys = [2.5e23 + 3.65*x^(3.5) for x in xs]
          calc = GaPSE.spectral_index(xs, ys; N=1, con=true)
          @test all([isapprox(3.5, c; rtol=1e-3) for c in calc])
     end

     @testset "fourth" begin
          xs = 10 .^ range(4,6, length=100)
          ys = [2.5e99 + 3.65*x^(3.5) for x in xs]
          @test_throws AssertionError GaPSE.spectral_index(xs, ys; N=1, con=true)
          @test_throws AssertionError GaPSE.spectral_index(xs, ys; N=1, con=false)
     end

     @testset "fifth" begin
          xs = 10 .^ range(-6,-4, length=100)
          ys = [3.65*x^(3.5) for x in xs]
          calc = GaPSE.spectral_index(xs, ys; N=1, con=false)
          @test all([isapprox(1.0, 3.5/c; rtol=1e-2) 
               for (x,c) in zip(xs[begin+1:end-1], calc[begin+1:end-1])])
     end

     @testset "sixth" begin
          xs = 10 .^ range(-6,-4, length=100)
          ys = [3.6e-22 + 3.65*x^(3.5) for x in xs]
          calc = GaPSE.spectral_index(xs, ys; N=1, con=true)
          @test all([isapprox(1.0, 3.5/c; rtol=1e-2) 
               for (x,c) in zip(xs[begin+1:end-1], calc[begin+1:end-1])])
     end

     @testset "seventh" begin
          xs = 10 .^ range(-6,-4, length=100)
          ys = [3.65*x^(-3.5) for x in xs]
          calc = GaPSE.spectral_index(xs, ys; N=1, con=false)
          @test all([isapprox(1.0, -3.5/c; rtol=1e-2) 
               for (x,c) in zip(xs[begin+1:end-1], calc[begin+1:end-1])])
     end

     @testset "eigth" begin
          xs = 10 .^ range(-6,-4, length=100)
          ys = [3.6e-22 + 3.65*x^(-3.5) for x in xs]
          calc = GaPSE.spectral_index(xs, ys; N=1, con=true)
          @test all([isapprox(1.0, -3.5/c; rtol=3e-2) 
               for (x,c) in zip(xs[begin+1:end-1], calc[begin+1:end-1])])
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

     @testset "third" begin
          xs = 10 .^ range(4, 6, length=100)
          ys = [3.54e6 * x^4.856 for x in xs]

          @test isapprox(GaPSE.mean_spectral_index(xs, ys; con = false), 4.856; atol = 3e-2)
     end

     @testset "fourth" begin
          xs = 10 .^ range(4, 6, length=100)
          ys = [8.3e30 - 3.54e6 * x^4.856 for x in xs]

          @test isapprox(GaPSE.mean_spectral_index(xs, ys; con = true)/4.856, 1.0; atol = 3e-2)
     end

     @testset "fifth" begin
          xs = 10 .^ range(-6, -4, length=100)
          ys = [3.54e6 * x^(-4.856) for x in xs]

          @test isapprox(GaPSE.mean_spectral_index(xs, ys; con = false)/(-4.856), 1.0; atol = 3e-2)
     end

     @testset "sixth" begin
          xs = 10 .^ range(-6, -4, length=100)
          ys = [8.3e-27 - 3.54e-6 * x^(-4.856) for x in xs]

          @test isapprox(GaPSE.mean_spectral_index(xs, ys; con = true)/(-4.856), 1.0; atol = 3e-2)
     end

end



##########################################################################################92



@testset "test power_law" begin
     x, y = 2.0, 31.0
     si, b, a = 3.0, 2.0, 15.0
     @test isapprox(y, GaPSE.power_law(x, si, b, a), rtol=1e-6)
end


@testset "test power_law_from_data" begin
     @testset "zeros" begin
          si, b, a = 2.69, 3.45, 0.0
          xs = 1:0.1:10
          ys = [a + b * x ^ (si)  for x in xs]
          p0=[1.0, 1.0]
          @test_throws AssertionError GaPSE.power_law_from_data(xs[begin:end-1], ys, p0; con = false)
          @test_throws AssertionError GaPSE.power_law_from_data(xs, ys[begin:end-1], p0; con = false)
          @test_throws AssertionError GaPSE.power_law_from_data(xs, ys, [1.0, 1.0, 1.0]; con = false)
          @test_throws AssertionError GaPSE.power_law_from_data(xs, ys, [1.0, 1.0]; con = true)
     end

     @testset "first" begin
          si, b, a = 2.69, 3.45, 0.0
          xs = 1:0.1:10
          ys = [a + b * x ^ (si)  for x in xs]
          p0 = [1.0, 1.0]
          c_si, c_b, c_a = GaPSE.power_law_from_data(xs, ys, p0; con = false)
          @test isapprox(si, c_si, rtol=1e-2)
          @test isapprox(b, c_b, rtol=1e-2)
          @test isapprox(a, c_a, rtol=1e-2)
     end

     @testset "second" begin
          si, b, a = -2.69, 3.45, 0.0
          xs = 10 .^ range(4, 6, length=100)
          ys = [a + b * x ^ (si)  for x in xs]
          p0 = [-1.0, 1.0]
          c_si, c_b, c_a = GaPSE.power_law_from_data(xs, ys, p0; con = false)
          @test isapprox(si, c_si, rtol=1e-2)
          @test isapprox(b, c_b, rtol=1e-2)
          @test isapprox(a, c_a, rtol=1e-2)
     end

     @testset "third" begin
          si, b, a = -2.69, 3.45, 13245.23
          xs = 10 .^ range(4, 6, length=100)
          ys = [a + b * x ^ (si)  for x in xs]
          p0 =  [-1.0, 1.0, 1.0]
          @test_throws AssertionError GaPSE.power_law_from_data(xs, ys, p0; con = true)
     end

     @testset "fourth" begin
          si, b, a = -2.65, 1.54e3, 2.34563e-15
          xs = 10 .^ range(4, 6, length=100)
          ys = [a + b * x ^ (si)  for x in xs]
          p0 = [-1.0, 1.0, 1.0]
          c_si, c_b, c_a = GaPSE.power_law_from_data(xs, ys, p0; con = true)
          @test isapprox(si, c_si, rtol=1e-2)
          @test isapprox(b, c_b, rtol=1e-2)
          @test isapprox(a, c_a, rtol=1e-2)
     end
end



##########################################################################################92



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


@testset "test expanded_I04_tilde" begin
     tab_I04_tildes = readdlm("datatest/I04_tilde_extended_NO_CONST.txt", comments=true)
     ss = convert(Vector{Float64}, tab_I04_tildes[:, 1])
     I04_tildes = convert(Vector{Float64}, tab_I04_tildes[:, 2])

     table_ips = readdlm(FILE_PS)
     ks = convert(Vector{Float64}, table_ips[:, 1])
     pks = convert(Vector{Float64}, table_ips[:, 2])
     PK = Spline1D(ks, pks)
     kmin, kmax = 1e-5, 1e3
     calc_I04_tildes = GaPSE.expanded_I04_tilde(PK, ss; kmin = kmin, kmax = kmax)

     @test all([isapprox(my, tr, rtol = 1e-3) for (my, tr) in zip(calc_I04_tildes, I04_tildes)])
end

