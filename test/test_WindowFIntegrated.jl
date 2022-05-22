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

kwargs_F_int_quad = Dict(
    :llim => kwargs_map_F_int[:llim], 
    :rlim => kwargs_map_F_int[:rlim], 
    :rtol => kwargs_map_F_int[:rtol], 
    :atol => kwargs_map_F_int[:atol],
)

kwargs_F_int_trap = Dict(
    :llim => kwargs_map_F_int[:llim], 
    :rlim => kwargs_map_F_int[:rlim], 
    :N => kwargs_map_F_int[:N], 
)

calc_μs = vcat([-1.0, -0.98, -0.95], 
    [μ for μ in -0.9:0.3:0.9], 
    [0.95, 0.98, 1.0]);

windF = GaPSE.WindowF(FILE_F_MAP); 

@testset "test integrated_F_quadgk" begin
     RTOL = 1e-2
     s_min, s_max = 148.1920001465757, 571.7022420258767
     z_min, z_max = 0.05, 0.20

     @test isapprox(GaPSE.integrated_F_quadgk(100, 0, s_min, s_max, windF; 
          kwargs_F_int_quad...), 2.12335e+09; rtol = RTOL)
     @test isapprox(GaPSE.integrated_F_quadgk(150, 0, s_min, s_max, windF; 
          kwargs_F_int_quad...), 1.93936e+09; rtol = RTOL)
     @test isapprox(GaPSE.integrated_F_quadgk(200, 0, s_min, s_max, windF; 
          kwargs_F_int_quad...), 1.68321e+09; rtol = RTOL)
     @test isapprox(GaPSE.integrated_F_quadgk(250, 0, s_min, s_max, windF; 
          kwargs_F_int_quad...), 1.45684e+09; rtol = RTOL)

     @test isapprox(GaPSE.integrated_F_quadgk(100, 0.8, s_min, s_max, windF; 
          kwargs_F_int_quad...), 1.42493e+09; rtol = RTOL)
     @test isapprox(GaPSE.integrated_F_quadgk(150, 0.8, s_min, s_max, windF; 
          kwargs_F_int_quad...), 1.06046e+09; rtol = RTOL)
     @test isapprox(GaPSE.integrated_F_quadgk(200, 0.8, s_min, s_max, windF; 
          kwargs_F_int_quad...), 7.53131e+08; rtol = RTOL)
     @test isapprox(GaPSE.integrated_F_quadgk(250, 0.8, s_min, s_max, windF; 
          kwargs_F_int_quad...), 4.82146e+08; rtol = RTOL)

     @test isapprox(GaPSE.integrated_F_quadgk(100, -0.8, s_min, s_max, windF; 
          kwargs_F_int_quad...), 2.14873e+09; rtol = RTOL)
     @test isapprox(GaPSE.integrated_F_quadgk(150, -0.8, s_min, s_max, windF;  
          kwargs_F_int_quad...), 2.08741e+09; rtol = RTOL)
     @test isapprox(GaPSE.integrated_F_quadgk(200, -0.8, s_min, s_max, windF;  
          kwargs_F_int_quad...), 1.98355e+09; rtol = RTOL)
     @test isapprox(GaPSE.integrated_F_quadgk(250, -0.8, s_min, s_max, windF; 
          kwargs_F_int_quad...), 1.93041e+09; rtol = RTOL)
end

@testset "test integrated_F_trapz" begin
     RTOL = 1e-2
     s_min, s_max = 148.1920001465757, 571.7022420258767
     z_min, z_max = 0.05, 0.20

     @test isapprox(GaPSE.integrated_F_trapz(100, 0, s_min, s_max, windF; 
          kwargs_F_int_trap...), 2.15022e+09; rtol = RTOL)
     @test isapprox(GaPSE.integrated_F_trapz(150, 0, s_min, s_max, windF; 
          kwargs_F_int_trap...), 1.93769e+09; rtol = RTOL)
     @test isapprox(GaPSE.integrated_F_trapz(200, 0, s_min, s_max, windF; 
          kwargs_F_int_trap...), 1.68718e+09; rtol = RTOL)
     @test isapprox(GaPSE.integrated_F_trapz(250, 0, s_min, s_max, windF; 
          kwargs_F_int_trap...), 1.41394e+09; rtol = RTOL)

     @test isapprox(GaPSE.integrated_F_trapz(100, 0.8, s_min, s_max, windF; 
          kwargs_F_int_trap...), 1.44814e+09; rtol = RTOL)
     @test isapprox(GaPSE.integrated_F_trapz(150, 0.8, s_min, s_max, windF; 
          kwargs_F_int_trap...), 1.04627e+09; rtol = RTOL)
     @test isapprox(GaPSE.integrated_F_trapz(200, 0.8, s_min, s_max, windF; 
          kwargs_F_int_trap...), 7.29342e+08; rtol = RTOL)
     @test isapprox(GaPSE.integrated_F_trapz(250, 0.8, s_min, s_max, windF; 
          kwargs_F_int_trap...), 4.84891e+08; rtol = RTOL)

     @test isapprox(GaPSE.integrated_F_trapz(100, -0.8, s_min, s_max, windF; 
          kwargs_F_int_trap...), 2.18180e+09; rtol = RTOL)
     @test isapprox(GaPSE.integrated_F_trapz(150, -0.8, s_min, s_max, windF;  
          kwargs_F_int_trap...), 2.05673e+09; rtol = RTOL)
     @test isapprox(GaPSE.integrated_F_trapz(200, -0.8, s_min, s_max, windF;  
          kwargs_F_int_trap...), 1.93539e+09; rtol = RTOL)
     @test isapprox(GaPSE.integrated_F_trapz(250, -0.8, s_min, s_max, windF; 
          kwargs_F_int_trap...), 1.90523e+09; rtol = RTOL)
end


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


@testset "test WindowFIntegrated: first convection" begin
     ss = 100 .* [0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3] .+ 30
     μs = [-1, 0, 1, -1, 0, 1, -1, 0, 1, -1, 0, 1]
     IFs = [0, 1, 2, 0, 2, 4, 0, 4, 8, 0, 8, 16]

     unique_ss = 100 .* [0, 1, 2, 3] .+ 30
     unique_μs = [-1, 0, 1]
     table_IFs = [0 1 2; 0 2 4; 0 4 8; 0 8 16]

     name = "test_WindowFIntegrated_fc.txt"
     isfile(name) && rm(name)
     open(name, "w") do io
          println(io, "# line of comment")
          println(io, "# another one")
          for (s, μ, IF) in zip(ss, μs, IFs)
               println(io, "$s \t $μ \t $IF")
          end
     end

     F_fc = GaPSE.WindowFIntegrated(name)

     @test size(F_fc.ss) == size(unique_ss)
     @test size(F_fc.μs) == size(unique_μs)
     @test size(F_fc.IFs) == size(table_IFs)
     @test all(F_fc.ss .== unique_ss)
     @test all(F_fc.μs .== unique_μs)
     @test all(F_fc.IFs .== table_IFs)

     rm(name)
end

@testset "test WindowFIntegrated: second convection" begin
     ss = 100 .* [0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3] .+ 30 
     μs = [-1, -1, -1, -1, 0, 0, 0, 0, 1, 1, 1, 1]
     IFs = [0, 0, 0, 0, 1, 2, 4, 8, 2, 4, 8, 16]

     unique_ss = 100 .* [0, 1, 2, 3] .+ 30
     unique_μs = [-1, 0, 1]
     table_IFs = [0 1 2; 0 2 4; 0 4 8; 0 8 16]

     name = "test_WindowFIntegrated_sc.txt"
     isfile(name) && rm(name)
     open(name, "w") do io
          println(io, "# line of comment")
          println(io, "# another one")
          for (s, μ, IF) in zip(ss, μs, IFs)
               println(io, "$s \t $μ \t $IF")
          end
     end

     F_sc = GaPSE.WindowFIntegrated(name)

     @test size(F_sc.ss) == size(unique_ss)
     @test size(F_sc.μs) == size(unique_μs)
     @test size(F_sc.IFs) == size(table_IFs)
     @test all(F_sc.ss .== unique_ss)
     @test all(F_sc.μs .== unique_μs)
     @test all(F_sc.IFs .== table_IFs)

     rm(name)
end


