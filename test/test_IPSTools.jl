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
     tab_pk = readdlm("datatest/file_pk.txt", comments = true)
     ks = convert(Vector{Float64}, tab_pk[:, 1])
     k_min, k_max = ks[begin], ks[end]
     pks = convert(Vector{Float64}, tab_pk[:, 2])
     PK = Spline1D(ks, pks)

     @testset "first" begin
          ips = GaPSE.InputPS("datatest/file_pk.txt"; expand = false)
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
     ips = GaPSE.InputPS("datatest/file_pk.txt"; expand = false)
     tools = GaPSE.IPSTools(ips)
     kmin, kmax = ips.ks[begin], ips.ks[end]
     PK = Spline1D(ips.ks, ips.pks)

     @testset "test sigmas" begin
          σ_0, _ = quadgk(q -> PK(q) * q^2 / (2 * π^2), kmin, kmax)
          σ_1, _ = quadgk(q -> PK(q) * q / (2 * π^2), kmin, kmax)
          σ_2, _ = quadgk(q -> PK(q) / (2 * π^2), kmin, kmax)
          σ_3, _ = quadgk(q -> PK(q) / (2 * π^2 * q), kmin, kmax)

          @test isapprox(tools.σ_0, σ_0)
          @test isapprox(tools.σ_1, σ_1)
          @test isapprox(tools.σ_2, σ_2)
          @test isapprox(tools.σ_3, σ_3)
     end

end