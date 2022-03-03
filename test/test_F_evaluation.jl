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

@testset "test F" begin
     @test isapprox(GaPSE.F(0, 0)[1], 39.0406; rtol = 1e-2)
     @test isapprox(GaPSE.F(1, 0)[1], 29.25801; rtol = 1e-2)
     @test isapprox(GaPSE.F(2, 0)[1], 24.97067; rtol = 1e-2)
     @test isapprox(GaPSE.F(3, 0)[1], 23.51367; rtol = 1e-2)

     @test isapprox(GaPSE.F(0, -0.8)[1], 38.89266; rtol = 1e-2)
     @test isapprox(GaPSE.F(1, -0.8)[1], 23.35162; rtol = 1e-2)
     @test isapprox(GaPSE.F(2, -0.8)[1], 11.83636; rtol = 1e-2)
     @test isapprox(GaPSE.F(3, -0.8)[1], 10.90119; rtol = 1e-2)

     @test isapprox(GaPSE.F(0, 0.8)[1], 38.89261; rtol = 1e-2)
     @test isapprox(GaPSE.F(1, 0.8)[1], 34.85789; rtol = 1e-2)
     @test isapprox(GaPSE.F(2, 0.8)[1], 33.54063; rtol = 1e-2)
     @test isapprox(GaPSE.F(3, 0.8)[1], 32.91128; rtol = 1e-2)
end


@testset "test F_map first method" begin
     name = "datatest/F_first_method.txt"
     output = "F_first_output.txt"

     @testset "zeros" begin
          @test_throws AssertionError GaPSE.F_map(0.25, 0.25; x1 = -0.5)
          @test_throws AssertionError GaPSE.F_map(0.25, 0.25; x1 = 1.0, x2 = 0.5)
          @test_throws AssertionError GaPSE.F_map(0.25, 0.25; x2 = 11.0)
          @test_throws AssertionError GaPSE.F_map(0.25, 0.25; μ1 = -1.5)
          @test_throws AssertionError GaPSE.F_map(0.25, 0.25; μ1 = -0.9, μ2 = -0.95)
          @test_throws AssertionError GaPSE.F_map(0.25, 0.25; μ2 = 1.5)

          @test_throws AssertionError GaPSE.F_map(0.0, 0.25)
          @test_throws AssertionError GaPSE.F_map(-1.0, 0.25)
          @test_throws AssertionError GaPSE.F_map(2.0, 0.25)
          @test_throws AssertionError GaPSE.F_map(0.25, 0.0)
          @test_throws AssertionError GaPSE.F_map(0.25, -1.0)
          @test_throws AssertionError GaPSE.F_map(0.25, 2.0)
     end

     GaPSE.F_map(0.25, 0.25;
          out = output, x1 = 0, x2 = 3, μ1 = -1, μ2 = 1, 
          rtol = 5e-3, atol = 1e-2)

     @testset "first" begin
          table_output_F = readdlm(output, comments = true)
          output_xs = convert(Vector{Float64}, table_output_F[:, 1])
          output_μs = convert(Vector{Float64}, table_output_F[:, 2])
          output_Fs = convert(Vector{Float64}, table_output_F[:, 3])

          table_F = readdlm(name, comments = true)
          xs = convert(Vector{Float64}, table_F[:, 1])
          μs = convert(Vector{Float64}, table_F[:, 2])
          Fs = convert(Vector{Float64}, table_F[:, 3])

          @test all([x1 ≈ x2 for (x1, x2) in zip(xs, output_xs)])
          @test all([μ1 ≈ μ2 for (μ1, μ2) in zip(μs, output_μs)])
          @test all([F1 ≈ F2 for (F1, F2) in zip(Fs, output_Fs)])
     end

     @testset "second" begin
          table_output_F = GaPSE.WindowF(output)
          output_xs = table_output_F.xs
          output_μs = table_output_F.μs
          output_Fs = table_output_F.Fs

          table_F = GaPSE.WindowF(name)
          xs = table_F.xs
          μs = table_F.μs
          Fs = table_F.Fs

          @test all([x1 ≈ x2 for (x1, x2) in zip(xs, output_xs)])
          @test all([μ1 ≈ μ2 for (μ1, μ2) in zip(μs, output_μs)])
          @test all([F1 ≈ F2 for (F1, F2) in zip(Fs, output_Fs)])
     end

     rm(output)
end


@testset "test F_map second method" begin
     name = "datatest/F_second_method.txt"
     output = "F_second_output.txt"

     calc_xs = [x for x in 0:0.25:3]
     calc_μs = vcat([-1.0, -0.98, -0.95], [μ for μ in -0.9:0.1:0.9], [0.95, 0.98, 1.0])

     @testset "zeros" begin
          @test_throws AssertionError GaPSE.F_map([1.0 for i in 1:10], calc_μs)
          @test_throws AssertionError GaPSE.F_map(calc_xs, [0.5 for i in 1:10])
          @test_throws AssertionError GaPSE.F_map([1.0, 2.0, 100.0], calc_μs)
          @test_throws AssertionError GaPSE.F_map(calc_xs, [-1.5, -0.99, 0.0, 0.99, 1.5])
          @test_throws AssertionError GaPSE.F_map(reverse(calc_xs), calc_μs)
          @test_throws AssertionError GaPSE.F_map(calc_xs, reverse(calc_μs))
     end

     GaPSE.F_map(calc_xs, calc_μs; out = output, rtol = 5e-3, atol = 1e-2)

     @testset "first" begin
          table_output_F = readdlm(output, comments = true)
          output_xs = convert(Vector{Float64}, table_output_F[:, 1])
          output_μs = convert(Vector{Float64}, table_output_F[:, 2])
          output_Fs = convert(Vector{Float64}, table_output_F[:, 3])

          table_F = readdlm(name, comments = true)
          xs = convert(Vector{Float64}, table_F[:, 1])
          μs = convert(Vector{Float64}, table_F[:, 2])
          Fs = convert(Vector{Float64}, table_F[:, 3])

          @test all([x1 ≈ x2 for (x1, x2) in zip(xs, output_xs)])
          @test all([μ1 ≈ μ2 for (μ1, μ2) in zip(μs, output_μs)])
          @test all([F1 ≈ F2 for (F1, F2) in zip(Fs, output_Fs)])
     end

     @testset "second" begin
          table_output_F = GaPSE.WindowF(output)
          output_xs = table_output_F.xs
          output_μs = table_output_F.μs
          output_Fs = table_output_F.Fs

          table_F = GaPSE.WindowF(name)
          xs = table_F.xs
          μs = table_F.μs
          Fs = table_F.Fs

          @test all([x1 ≈ x2 for (x1, x2) in zip(xs, output_xs)])
          @test all([μ1 ≈ μ2 for (μ1, μ2) in zip(μs, output_μs)])
          @test all([F1 ≈ F2 for (F1, F2) in zip(Fs, output_Fs)])
     end

     rm(output)
end