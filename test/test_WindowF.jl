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

kwargs_F_hcub = Dict(
    :θ_max => π / 2.0,
    :tolerance => 1e-10,
    :rtol => 1e-2,
    :atol => 1e-3,
)

kwargs_F_trap = Dict(
    :θ_max => π / 2.0::Float64,
    :tolerance => 1e-8::Float64,
    :N => 1000::Int64,
    :en => 1.0::Float64,
)

kwargs_print_map_F = Dict(
    :θ_max => π / 2.0,
    :tolerance => 1e-10,
    :rtol => 1e-2,
    :atol => 1e-3,
    :pr => true,
)


##########################################################################################92


@testset "test F_hcub" begin
    RTOL = 1e-2

    @test isapprox(GaPSE.F_hcub(0, 0; kwargs_F_hcub...)[1], 39.0406; rtol=RTOL)
    @test isapprox(GaPSE.F_hcub(1, 0; kwargs_F_hcub...)[1], 29.25801; rtol=RTOL)
    @test isapprox(GaPSE.F_hcub(2, 0; kwargs_F_hcub...)[1], 25.28027; rtol=RTOL)
    @test isapprox(GaPSE.F_hcub(3, 0; kwargs_F_hcub...)[1], 23.51367; rtol=RTOL)

    @test isapprox(GaPSE.F_hcub(0, -0.8; kwargs_F_hcub...)[1], 38.89266; rtol=RTOL)
    @test isapprox(GaPSE.F_hcub(1, -0.8; kwargs_F_hcub...)[1], 23.35162; rtol=RTOL)
    @test isapprox(GaPSE.F_hcub(2, -0.8; kwargs_F_hcub...)[1], 11.83636; rtol=RTOL)
    @test isapprox(GaPSE.F_hcub(3, -0.8; kwargs_F_hcub...)[1], 10.90119; rtol=RTOL)

    @test isapprox(GaPSE.F_hcub(0, 0.8; kwargs_F_hcub...)[1], 38.89261; rtol=RTOL)
    @test isapprox(GaPSE.F_hcub(1, 0.8; kwargs_F_hcub...)[1], 34.85789; rtol=RTOL)
    @test isapprox(GaPSE.F_hcub(2, 0.8; kwargs_F_hcub...)[1], 33.54063; rtol=RTOL)
    @test isapprox(GaPSE.F_hcub(3, 0.8; kwargs_F_hcub...)[1], 32.91128; rtol=RTOL)
end

@testset "test F_trap" begin
    RTOL = 1e-4

    @test isapprox(GaPSE.F_trap(0, 0; kwargs_F_trap...), 39.40821; rtol=RTOL)
    @test isapprox(GaPSE.F_trap(1, 0; kwargs_F_trap...), 29.59887; rtol=RTOL)
    @test isapprox(GaPSE.F_trap(2, 0; kwargs_F_trap...), 25.55135; rtol=RTOL)
    @test isapprox(GaPSE.F_trap(3, 0; kwargs_F_trap...), 23.77376; rtol=RTOL)

    @test isapprox(GaPSE.F_trap(0, -0.8; kwargs_F_trap...), 39.41779; rtol=RTOL)
    @test isapprox(GaPSE.F_trap(1, -0.8; kwargs_F_trap...), 23.77100; rtol=RTOL)
    @test isapprox(GaPSE.F_trap(2, -0.8; kwargs_F_trap...), 13.87924; rtol=RTOL)
    @test isapprox(GaPSE.F_trap(3, -0.8; kwargs_F_trap...), 11.40667; rtol=RTOL)

    @test isapprox(GaPSE.F_trap(0, 0.8; kwargs_F_trap...), 39.41779; rtol=RTOL)
    @test isapprox(GaPSE.F_trap(1, 0.8; kwargs_F_trap...), 35.42117; rtol=RTOL)
    @test isapprox(GaPSE.F_trap(2, 0.8; kwargs_F_trap...), 34.04887; rtol=RTOL)
    @test isapprox(GaPSE.F_trap(3, 0.8; kwargs_F_trap...), 33.32257; rtol=RTOL)
end


##########################################################################################92


@testset "test print_map_F with trap first method" begin
    name = "datatest/WindowF/F_trap_first_method.txt"
    output = "F_trap_first_output.txt"

    @testset "zeros" begin
        @test_throws AssertionError GaPSE.print_map_F(output, 0.25, 0.25; x1=-0.5)
        @test_throws AssertionError GaPSE.print_map_F(output, 0.25, 0.25; x1=1.0, x2=0.5)
        @test_throws AssertionError GaPSE.print_map_F(output, 0.25, 0.25; x2=11.0)
        @test_throws AssertionError GaPSE.print_map_F(output, 0.25, 0.25; μ1=-1.5)
        @test_throws AssertionError GaPSE.print_map_F(output, 0.25, 0.25; μ1=-0.9, μ2=-0.95)
        @test_throws AssertionError GaPSE.print_map_F(output, 0.25, 0.25; μ2=1.5)
        @test_throws AssertionError GaPSE.print_map_F(output, 0.25, 0.25; alg=:try)

        @test_throws AssertionError GaPSE.print_map_F(output, 0.0, 0.25)
        @test_throws AssertionError GaPSE.print_map_F(output, -1.0, 0.25)
        @test_throws AssertionError GaPSE.print_map_F(output, 2.0, 0.25)
        @test_throws AssertionError GaPSE.print_map_F(output, 0.25, 0.0)
        @test_throws AssertionError GaPSE.print_map_F(output, 0.25, -1.0)
        @test_throws AssertionError GaPSE.print_map_F(output, 0.25, 2.0)
    end

    GaPSE.print_map_F(output, 0.25, 0.25;
        alg=:trap, x1=0, x2=3, μ1=-1, μ2=1,
        Fmap_opts=kwargs_F_trap)

    @testset "first" begin
        table_output_F = readdlm(output, comments=true)
        output_xs = convert(Vector{Float64}, table_output_F[:, 1])
        output_μs = convert(Vector{Float64}, table_output_F[:, 2])
        output_Fs = convert(Vector{Float64}, table_output_F[:, 3])

        table_F = readdlm(name, comments=true)
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


@testset "test print_map_F with trap second method" begin
    name = "datatest/WindowF/F_trap_second_method.txt"
    output = "F_trap_second_output.txt"

    calc_xs = [x for x in 0:0.25:3]
    calc_μs = vcat([-1.0, -0.98, -0.95], [μ for μ in -0.9:0.1:0.9], [0.95, 0.98, 1.0])

    @testset "zeros" begin
        @test_throws AssertionError GaPSE.print_map_F(output, [1.0 for i in 1:10], calc_μs)
        @test_throws AssertionError GaPSE.print_map_F(output, calc_xs, [0.5 for i in 1:10])
        @test_throws AssertionError GaPSE.print_map_F(output, [1.0, 2.0, 100.0], calc_μs)
        @test_throws AssertionError GaPSE.print_map_F(output, calc_xs, [-1.5, -0.99, 0.0, 0.99, 1.5])
        @test_throws AssertionError GaPSE.print_map_F(output, reverse(calc_xs), calc_μs)
        @test_throws AssertionError GaPSE.print_map_F(output, calc_xs, reverse(calc_μs))
        @test_throws AssertionError GaPSE.print_map_F(output, calc_xs, calc_μs; alg=:try)
    end

    GaPSE.print_map_F(output, calc_xs, calc_μs;
        alg=:trap, Fmap_opts=kwargs_F_trap)

    @testset "first" begin
        table_output_F = readdlm(output, comments=true)
        output_xs = convert(Vector{Float64}, table_output_F[:, 1])
        output_μs = convert(Vector{Float64}, table_output_F[:, 2])
        output_Fs = convert(Vector{Float64}, table_output_F[:, 3])

        table_F = readdlm(name, comments=true)
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


@testset "test print_map_F with hcub first method" begin
    name = "datatest/WindowF/F_hcub_first_method.txt"
    output = "F_hcub_first_output.txt"

    @testset "zeros" begin
        @test_throws AssertionError GaPSE.print_map_F(output, 0.25, 0.25; x1=-0.5)
        @test_throws AssertionError GaPSE.print_map_F(output, 0.25, 0.25; x1=1.0, x2=0.5)
        @test_throws AssertionError GaPSE.print_map_F(output, 0.25, 0.25; x2=11.0)
        @test_throws AssertionError GaPSE.print_map_F(output, 0.25, 0.25; μ1=-1.5)
        @test_throws AssertionError GaPSE.print_map_F(output, 0.25, 0.25; μ1=-0.9, μ2=-0.95)
        @test_throws AssertionError GaPSE.print_map_F(output, 0.25, 0.25; μ2=1.5)
        @test_throws AssertionError GaPSE.print_map_F(output, 0.25, 0.25; alg=:try)

        @test_throws AssertionError GaPSE.print_map_F(output, 0.0, 0.25)
        @test_throws AssertionError GaPSE.print_map_F(output, -1.0, 0.25)
        @test_throws AssertionError GaPSE.print_map_F(output, 2.0, 0.25)
        @test_throws AssertionError GaPSE.print_map_F(output, 0.25, 0.0)
        @test_throws AssertionError GaPSE.print_map_F(output, 0.25, -1.0)
        @test_throws AssertionError GaPSE.print_map_F(output, 0.25, 2.0)
    end

    GaPSE.print_map_F(output, 0.25, 0.25;
        alg=:hcub, x1=0, x2=3, μ1=-1, μ2=1,
        Fmap_opts=kwargs_F_hcub)

    @testset "first" begin
        table_output_F = readdlm(output, comments=true)
        output_xs = convert(Vector{Float64}, table_output_F[:, 1])
        output_μs = convert(Vector{Float64}, table_output_F[:, 2])
        output_Fs = convert(Vector{Float64}, table_output_F[:, 3])

        table_F = readdlm(name, comments=true)
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


@testset "test print_map_F with trap second method" begin
    name = "datatest/WindowF/F_hcub_second_method.txt"
    output = "F_hcub_second_output.txt"

    calc_xs = [x for x in 0:0.25:3]
    calc_μs = vcat([-1.0, -0.98, -0.95], [μ for μ in -0.9:0.1:0.9], [0.95, 0.98, 1.0])

    @testset "zeros" begin
        @test_throws AssertionError GaPSE.print_map_F(output, [1.0 for i in 1:10], calc_μs)
        @test_throws AssertionError GaPSE.print_map_F(output, calc_xs, [0.5 for i in 1:10])
        @test_throws AssertionError GaPSE.print_map_F(output, [1.0, 2.0, 100.0], calc_μs)
        @test_throws AssertionError GaPSE.print_map_F(output, calc_xs, [-1.5, -0.99, 0.0, 0.99, 1.5])
        @test_throws AssertionError GaPSE.print_map_F(output, reverse(calc_xs), calc_μs)
        @test_throws AssertionError GaPSE.print_map_F(output, calc_xs, reverse(calc_μs))
        @test_throws AssertionError GaPSE.print_map_F(output, calc_xs, calc_μs; alg=:try)
    end

    GaPSE.print_map_F(output, calc_xs, calc_μs;
        alg=:hcub, Fmap_opts=kwargs_F_hcub)

    @testset "first" begin
        table_output_F = readdlm(output, comments=true)
        output_xs = convert(Vector{Float64}, table_output_F[:, 1])
        output_μs = convert(Vector{Float64}, table_output_F[:, 2])
        output_Fs = convert(Vector{Float64}, table_output_F[:, 3])

        table_F = readdlm(name, comments=true)
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


##########################################################################################92


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


##########################################################################################92


@testset "test spline windowF" begin
    wf_trap = GaPSE.WindowF("datatest/WindowF/F_trap_first_method.txt")
    spline_trap(s, μ) = GaPSE.spline_F(s, μ, wf_trap)

    @test isapprox(spline_trap(0.5, -0.8), 2.66590e+01; rtol=1e-5)
    @test isapprox(spline_trap(1.9, -0.8), 1.25007e+01; rtol=1e-5)
    @test isapprox(spline_trap(0.5, 0.1), 3.39167e+01; rtol=1e-5)
    @test isapprox(spline_trap(1.9, 0.1), 2.67966e+01; rtol=1e-5)
end

@testset "test second method print_map_F" begin
    name_trap = "datatest/WindowF/F_trap_first_method.txt"
    orig_wf = GaPSE.WindowF(name_trap)

    print_trap = "test_print_wf.txt"
    isfile(print_trap) && rm(print_trap)
    GaPSE.print_map_F(print_trap, orig_wf)
    other_wf = GaPSE.WindowF(print_trap)

    @test all([x1 ≈ x2 for (x1, x2) in zip(orig_wf.xs, other_wf.xs)])
    @test all([μ1 ≈ μ2 for (μ1, μ2) in zip(orig_wf.μs, other_wf.μs)])
    @test all([IF1 ≈ IF2 for (IF1, IF2) in zip(orig_wf.Fs, other_wf.Fs)])
    rm(print_trap)
end