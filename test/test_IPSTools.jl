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


@testset "test InputPS" begin
     tab_pk = readdlm(FILE_PS, comments = true)
     ks = convert(Vector{Float64}, tab_pk[:, 1])
     k_min, k_max = ks[begin], ks[end]
     pks = convert(Vector{Float64}, tab_pk[:, 2])
     PK = Spline1D(ks, pks)

     @testset "first" begin
          ips = GaPSE.InputPS(FILE_PS; expand = false)
          @test all([k1 ≈ k2 for (k1, k2) in zip(ips.ks, ks)])
          @test all([pk1 ≈ pk2 for (pk1, pk2) in zip(ips.pks, pks)])
     end

     @testset "second" begin
          ips = GaPSE.InputPS(ks, pks; expand = false)
          @test all([k1 ≈ k2 for (k1, k2) in zip(ips.ks, ks)])
          @test all([pk1 ≈ pk2 for (pk1, pk2) in zip(ips.pks, pks)])
     end
end


@testset "test IPSTools" begin
     ips = GaPSE.InputPS(FILE_PS; expand = false)
     k_min, k_max = 1e-6, 10.0

     tools = GaPSE.IPSTools(ips; k_min = k_min, k_max = k_max)
     PK = Spline1D(ips.ks, ips.pks)

     tab_Is = readdlm(FILE_ILN, comments = true)
     true_ss = convert(Vector{Float64}, tab_Is[2:end, 1])
     ss = true_ss[1e-3 .< true_ss .< 1e3]
     I00s = convert(Vector{Float64}, tab_Is[2:end, 2])[1e-3 .< true_ss .< 1e3]
     I20s = convert(Vector{Float64}, tab_Is[2:end, 3])[1e-3 .< true_ss .< 1e3]
     I40s = convert(Vector{Float64}, tab_Is[2:end, 4])[1e-3 .< true_ss .< 1e3]
     I02s = convert(Vector{Float64}, tab_Is[2:end, 5])[1e-3 .< true_ss .< 1e3] ./ ss .^ 2
     I22s = convert(Vector{Float64}, tab_Is[2:end, 6])[1e-3 .< true_ss .< 1e3] ./ ss .^ 2
     I31s = convert(Vector{Float64}, tab_Is[2:end, 7])[1e-3 .< true_ss .< 1e3] ./ ss
     I11s = convert(Vector{Float64}, tab_Is[2:end, 8])[1e-3 .< true_ss .< 1e3] ./ ss

     @testset "test sigmas" begin
          σ_0, _ = quadgk(q -> PK(q) * q^2 / (2 * π^2), k_min, k_max)
          σ_1, _ = quadgk(q -> PK(q) * q / (2 * π^2), k_min, k_max)
          σ_2, _ = quadgk(q -> PK(q) / (2 * π^2), k_min, k_max)
          σ_3, _ = quadgk(q -> PK(q) / (2 * π^2 * q), k_min, k_max)

          @test isapprox(tools.σ_0, σ_0) && isapprox(σ_0, 15.593462966741178)
          @test isapprox(tools.σ_1, σ_1) && isapprox(σ_1, 15.074895881392285)
          @test isapprox(tools.σ_2, σ_2) && isapprox(σ_2, 100.85852368830221)
          @test isapprox(tools.σ_3, σ_3) && isapprox(σ_3, 3734.849690975012)
     end

     @testset "test Iln" begin
          @test all([isapprox(tools.I00(s), i, rtol = 5e-2) for (s, i) in zip(ss, I00s)])
          @test all([isapprox(tools.I20(s), i, rtol = 1e-1) for (s, i) in zip(ss, I20s)])
          @test all([isapprox(tools.I40(s), i, rtol = 1e-2) for (s, i) in zip(ss, I40s)])
          @test all([isapprox(tools.I02(s), i, rtol = 1e-1) for (s, i) in zip(ss, I02s)])
          @test all([isapprox(tools.I22(s), i, rtol = 5e-2) for (s, i) in zip(ss, I22s)])
          @test all([isapprox(tools.I31(s), i, rtol = 1e-2) for (s, i) in zip(ss, I31s)])
          @test all([isapprox(tools.I11(s), i, rtol = 5e-2) for (s, i) in zip(ss, I11s)])
     end

end