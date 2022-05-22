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

kwargs_map_F_int = Dict(
    :ss_start => 100, :ss_stop => 500, 
    :ss_step => 50, :llim => 0.0, :rlim => Inf, 
    :rtol => 5e-2, :atol => 0.0, :N => 300, :pr => true,
);

calc_μs = vcat([-1.0, -0.98, -0.95], 
    [μ for μ in -0.9:0.3:0.9], 
    [0.95, 0.98, 1.0]);
#=
@testset "test integrated_F_quadgk" begin
     RTOL = 1e-2
     s_min, s_max = 148.1920001465757, 571.7022420258767

     @test isapprox(GaPSE.integrated_F_quadgk(0, 0; kwargs_map_F_int...)[1], 39.0406; rtol = RTOL)
     @test isapprox(GaPSE.integrated_F_quadgk(1, 0; kwargs_map_F_int...)[1], 29.25801; rtol = RTOL)
     @test isapprox(GaPSE.integrated_F_quadgk(2, 0; kwargs_map_F_int...)[1], 25.28027; rtol = RTOL)
     @test isapprox(GaPSE.integrated_F_quadgk(3, 0; kwargs_map_F_int...)[1], 23.51367; rtol = RTOL)

     @test isapprox(GaPSE.integrated_F_quadgk(0, -0.8; kwargs_map_F_int...)[1], 38.89266; rtol = RTOL)
     @test isapprox(GaPSE.integrated_F_quadgk(1, -0.8; kwargs_map_F_int...)[1], 23.35162; rtol = RTOL)
     @test isapprox(GaPSE.integrated_F_quadgk(2, -0.8; kwargs_map_F_int...)[1], 11.83636; rtol = RTOL)
     @test isapprox(GaPSE.integrated_F_quadgk(3, -0.8; kwargs_map_F_int...)[1], 10.90119; rtol = RTOL)

     @test isapprox(GaPSE.integrated_F_quadgk(0, 0.8; kwargs_map_F_int...)[1], 38.89261; rtol = RTOL)
     @test isapprox(GaPSE.integrated_F_quadgk(1, 0.8; kwargs_map_F_int...)[1], 34.85789; rtol = RTOL)
     @test isapprox(GaPSE.integrated_F_quadgk(2, 0.8; kwargs_map_F_int...)[1], 33.54063; rtol = RTOL)
     @test isapprox(GaPSE.integrated_F_quadgk(3, 0.8; kwargs_map_F_int...)[1], 32.91128; rtol = RTOL)
end

@testset "test integrated_F_trap" begin
     RTOL = 1e-4

     @test isapprox(GaPSE.integrated_F_trap(0, 0; kwargs_map_F_int...), 39.40821; rtol = RTOL)
     @test isapprox(GaPSE.integrated_F_trap(1, 0; kwargs_map_F_int...), 29.59887; rtol = RTOL)
     @test isapprox(GaPSE.integrated_F_trap(2, 0; kwargs_map_F_int...), 25.55135; rtol = RTOL)
     @test isapprox(GaPSE.integrated_F_trap(3, 0; kwargs_map_F_int...), 23.77376; rtol = RTOL)

     @test isapprox(GaPSE.integrated_F_trap(0, -0.8; kwargs_map_F_int...), 39.41779; rtol = RTOL)
     @test isapprox(GaPSE.integrated_F_trap(1, -0.8; kwargs_map_F_int...), 23.77100; rtol = RTOL)
     @test isapprox(GaPSE.integrated_F_trap(2, -0.8; kwargs_map_F_int...), 13.87924; rtol = RTOL)
     @test isapprox(GaPSE.integrated_F_trap(3, -0.8; kwargs_map_F_int...), 11.40667; rtol = RTOL)

     @test isapprox(GaPSE.integrated_F_trap(0, 0.8; kwargs_map_F_int...), 39.41779; rtol = RTOL)
     @test isapprox(GaPSE.integrated_F_trap(1, 0.8; kwargs_map_F_int...), 35.42117; rtol = RTOL)
     @test isapprox(GaPSE.integrated_F_trap(2, 0.8; kwargs_map_F_int...), 34.04887; rtol = RTOL)
     @test isapprox(GaPSE.integrated_F_trap(3, 0.8; kwargs_map_F_int...), 33.32257; rtol = RTOL)
end
=#


@testset "test print_map_IntegratedF trap" begin
     name_1 = "datatest/WindowFIntegrated/IntF_trap_1.txt";
     name_2 = "datatest/WindowFIntegrated/IntF_trap_2.txt";
     in = FILE_F_MAP;
     output_1 = "calc_IntF_trap_1.txt";
     output_2 = "calc_IntF_trap_2.txt";
     z_min, z_max = 0.05, 0.20;
     s_min, s_max = 148.1920001465757, 571.7022420258767;
     μs = vcat([-1.0, -0.98, -0.95], [μ for μ in -0.9:0.3:0.9], [0.95, 0.98, 1.0]);

     @testset "zeros" begin
          @test_throws AssertionError GaPSE.print_map_IntegratedF(in, output_1, s_min, s_max, μs; trap = true, ss_start = -1)
          @test_throws AssertionError GaPSE.print_map_IntegratedF(in, output_1, s_min, s_max, μs; trap = true, ss_start = 1, ss_stop = 0.5)
          @test_throws AssertionError GaPSE.print_map_IntegratedF(in, output_1, s_min, s_max, μs; trap = true, ss_step = 1e12)

          @test_throws AssertionError GaPSE.print_map_IntegratedF(in, output_1, -1, 0.0, μs; trap = true)
          @test_throws AssertionError GaPSE.print_map_IntegratedF(in, output_1, 1.0, 0.5, μs; trap = true)
          @test_throws AssertionError GaPSE.print_map_IntegratedF(in, output_1, 1.0, 1.0, μs; trap = true)
          @test_throws AssertionError GaPSE.print_map_IntegratedF(in, output_1, 0.25, 1.0, [-1.5, 0.0, 0.8]; trap = true)
     end

     GaPSE.print_map_IntegratedF(in, output_1, s_min, s_max, μs;
          trap = true, kwargs_map_F_int...)
     GaPSE.print_map_IntegratedF(in, output_2, z_min, z_max, μs, FILE_BACKGROUND;
          trap = true, kwargs_map_F_int...)

     @testset "first" begin
          table_output_F = readdlm(output_1, comments = true)
          output_ss = convert(Vector{Float64}, table_output_F[:, 1])
          output_μs = convert(Vector{Float64}, table_output_F[:, 2])
          output_IFs = convert(Vector{Float64}, table_output_F[:, 3])

          table_F = readdlm(name_1, comments = true)
          ss = convert(Vector{Float64}, table_F[:, 1])
          μs = convert(Vector{Float64}, table_F[:, 2])
          IFs = convert(Vector{Float64}, table_F[:, 3])

          @test all([s1 ≈ s2 for (s1, s2) in zip(ss, output_ss)])
          @test all([μ1 ≈ μ2 for (μ1, μ2) in zip(μs, output_μs)])
          @test all([IF1 ≈ IF2 for (IF1, IF2) in zip(IFs, output_IFs)])
     end

     @testset "second" begin
          table_output_F = readdlm(output_2, comments = true)
          output_ss = convert(Vector{Float64}, table_output_F[:, 1])
          output_μs = convert(Vector{Float64}, table_output_F[:, 2])
          output_IFs = convert(Vector{Float64}, table_output_F[:, 3])

          table_F = readdlm(name_2, comments = true)
          ss = convert(Vector{Float64}, table_F[:, 1])
          μs = convert(Vector{Float64}, table_F[:, 2])
          IFs = convert(Vector{Float64}, table_F[:, 3])

          @test all([s1 ≈ s2 for (s1, s2) in zip(ss, output_ss)])
          @test all([μ1 ≈ μ2 for (μ1, μ2) in zip(μs, output_μs)])
          @test all([IF1 ≈ IF2 for (IF1, IF2) in zip(IFs, output_IFs)])
     end

     rm(output_1)
     rm(output_2)
end


@testset "test print_map_IntegratedF quad" begin
     name_1 = "datatest/WindowFIntegrated/IntF_quad_1.txt";
     name_2 = "datatest/WindowFIntegrated/IntF_quad_2.txt";
     in = FILE_F_MAP;
     output_1 = "calc_IntF_quad_1.txt";
     output_2 = "calc_IntF_quad_2.txt";
     z_min, z_max = 0.05, 0.20;
     s_min, s_max = 148.1920001465757, 571.7022420258767;
     μs = vcat([-1.0, -0.98, -0.95], [μ for μ in -0.9:0.3:0.9], [0.95, 0.98, 1.0]);

     @testset "zeros" begin
          @test_throws AssertionError GaPSE.print_map_IntegratedF(in, output_1, s_min, s_max, μs; trap = false, ss_start = -1)
          @test_throws AssertionError GaPSE.print_map_IntegratedF(in, output_1, s_min, s_max, μs; trap = false, ss_start = 1, ss_stop = 0.5)
          @test_throws AssertionError GaPSE.print_map_IntegratedF(in, output_1, s_min, s_max, μs; trap = false, ss_step = 1e12)

          @test_throws AssertionError GaPSE.print_map_IntegratedF(in, output_1, -1, 0.0, μs; trap = false)
          @test_throws AssertionError GaPSE.print_map_IntegratedF(in, output_1, 1.0, 0.5, μs; trap = false)
          @test_throws AssertionError GaPSE.print_map_IntegratedF(in, output_1, 1.0, 1.0, μs; trap = false)
          @test_throws AssertionError GaPSE.print_map_IntegratedF(in, output_1, 0.25, 1.0, [-1.5, 0.0, 0.8]; trap = false)
     end

     GaPSE.print_map_IntegratedF(in, output_1, s_min, s_max, μs;
          trap = false, kwargs_map_F_int...)
     GaPSE.print_map_IntegratedF(in, output_2, z_min, z_max, μs, FILE_BACKGROUND;
          trap = false, kwargs_map_F_int...)

     @testset "first" begin
          table_output_F = readdlm(output_1, comments = true)
          output_ss = convert(Vector{Float64}, table_output_F[:, 1])
          output_μs = convert(Vector{Float64}, table_output_F[:, 2])
          output_IFs = convert(Vector{Float64}, table_output_F[:, 3])

          table_F = readdlm(name_1, comments = true)
          ss = convert(Vector{Float64}, table_F[:, 1])
          μs = convert(Vector{Float64}, table_F[:, 2])
          IFs = convert(Vector{Float64}, table_F[:, 3])

          @test all([s1 ≈ s2 for (s1, s2) in zip(ss, output_ss)])
          @test all([μ1 ≈ μ2 for (μ1, μ2) in zip(μs, output_μs)])
          @test all([IF1 ≈ IF2 for (IF1, IF2) in zip(IFs, output_IFs)])
     end

     @testset "second" begin
          table_output_F = readdlm(output_2, comments = true)
          output_ss = convert(Vector{Float64}, table_output_F[:, 1])
          output_μs = convert(Vector{Float64}, table_output_F[:, 2])
          output_IFs = convert(Vector{Float64}, table_output_F[:, 3])

          table_F = readdlm(name_2, comments = true)
          ss = convert(Vector{Float64}, table_F[:, 1])
          μs = convert(Vector{Float64}, table_F[:, 2])
          IFs = convert(Vector{Float64}, table_F[:, 3])

          @test all([s1 ≈ s2 for (s1, s2) in zip(ss, output_ss)])
          @test all([μ1 ≈ μ2 for (μ1, μ2) in zip(μs, output_μs)])
          @test all([IF1 ≈ IF2 for (IF1, IF2) in zip(IFs, output_IFs)])
     end

     rm(output_1)
     rm(output_2)
end


#@test 1==2

@testset "test WindowF: first convection" begin
     xs = [0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3]
     μs = [-1, 0, 1, -1, 0, 1, -1, 0, 1, -1, 0, 1]
     Fs = [0, 1, 2, 0, 2, 4, 0, 4, 8, 0, 8, 16]

     unique_xs = [0, 1, 2, 3]
     unique_μs = [-1, 0, 1]
     table_Fs = [0 1 2; 0 2 4; 0 4 8; 0 8 16]

     name = "test_WindowF_fc.txt"
     isfile(name) && rm(name)
     open(name, "w") do io
          println(io, "# line of comment")
          println(io, "# another one")
          for (x, μ, F) in zip(xs, μs, Fs)
               println(io, "$x \t $μ \t $F")
          end
     end

     F_fc = GaPSE.WindowF(name)

     @test size(F_fc.xs) == size(unique_xs)
     @test size(F_fc.μs) == size(unique_μs)
     @test size(F_fc.Fs) == size(table_Fs)
     @test all(F_fc.xs .== unique_xs)
     @test all(F_fc.μs .== unique_μs)
     @test all(F_fc.Fs .== table_Fs)

     rm(name)
end

@testset "test WindowF: second convection" begin
     xs = [0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3]
     μs = [-1, -1, -1, -1, 0, 0, 0, 0, 1, 1, 1, 1]
     Fs = [0, 0, 0, 0, 1, 2, 4, 8, 2, 4, 8, 16]

     unique_xs = [0, 1, 2, 3]
     unique_μs = [-1, 0, 1]
     table_Fs = [0 1 2; 0 2 4; 0 4 8; 0 8 16]

     name = "test_WindowF_sc.txt"
     isfile(name) && rm(name)
     open(name, "w") do io
          println(io, "# line of comment")
          println(io, "# another one")
          for (x, μ, F) in zip(xs, μs, Fs)
               println(io, "$x \t $μ \t $F")
          end
     end

     F_sc = GaPSE.WindowF(name)

     @test size(F_sc.xs) == size(unique_xs)
     @test size(F_sc.μs) == size(unique_μs)
     @test size(F_sc.Fs) == size(table_Fs)
     @test all(F_sc.xs .== unique_xs)
     @test all(F_sc.μs .== unique_μs)
     @test all(F_sc.Fs .== table_Fs)

     rm(name)
end
