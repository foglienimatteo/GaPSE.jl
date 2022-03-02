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
          ips = GaPSE.InputPS(FILE_PS; 
               fit_left_min = 1e-6, fit_left_max = 3e-6,
               fit_right_min = 1e1, fit_right_max = 2e1)
          
          @test ips.l_si ≈ 0.9599998007678217
          @test ips.l_b ≈ 3.0123089778963984e6
          @test ips.l_a ≈ 0.0
          @test ips.left ≈ 1e-6

          @test ips.r_si ≈ -2.6531284143857454
          @test ips.r_b ≈ 68.30720789490483
          @test ips.r_a ≈ 0.0
          @test ips.right ≈ 2e1

          @test all([isapprox(ips(k), pk, rtol=1e-3)  for (k, pk) in 
               zip(ks[ks .> 0.7e-7], pks[ks .> 0.7e-7])])
     end

     @testset "second" begin
          ips = GaPSE.InputPS(ks, pks; 
               fit_left_min = 1e-6, fit_left_max = 3e-6,
               fit_right_min = 1e1, fit_right_max = 2e1)

          @test ips.l_si ≈ 0.9599998007678217
          @test ips.l_b ≈ 3.0123089778963984e6
          @test ips.l_a ≈ 0.0
          @test ips.left ≈ 1e-6

          @test ips.r_si ≈ -2.6531284143857454
          @test ips.r_b ≈ 68.30720789490483
          @test ips.r_a ≈ 0.0
          @test ips.right ≈ 2e1

          @test all([isapprox(ips(k), pk, rtol=1e-3)  for (k, pk) in 
               zip(ks[ks .> 0.7e-7], pks[ks .> 0.7e-7])])
     end
end


@testset "test IPSTools" begin
     ips = GaPSE.InputPS(FILE_PS; )
     k_min, k_max = 1e-6, 10.0

     tools = GaPSE.IPSTools(ips; k_min = k_min, k_max = k_max, N = 1024,
          fit_min = 0.05, fit_max = 0.5, con = true)
     #PK = Spline1D(ips.ks, ips.pks)
     PK = ips

     tab_Is = readdlm(FILE_ILN, comments = true)
     ss = convert(Vector{Float64}, tab_Is[:, 1])
     I00s = convert(Vector{Float64}, tab_Is[:, 2])
     I20s = convert(Vector{Float64}, tab_Is[:, 3])
     I40s = convert(Vector{Float64}, tab_Is[:, 4])
     I02s = convert(Vector{Float64}, tab_Is[:, 5])
     I22s = convert(Vector{Float64}, tab_Is[:, 6])
     I31s = convert(Vector{Float64}, tab_Is[:, 7])
     I11s = convert(Vector{Float64}, tab_Is[:, 8])
     I13s = convert(Vector{Float64}, tab_Is[:, 9])
     I04_tildes = convert(Vector{Float64}, tab_Is[:, 10])
     #=
     true_ss = convert(Vector{Float64}, tab_Is[2:end, 1])
     ss = true_ss[1e-3 .< true_ss .< 1e3]
     I00s = convert(Vector{Float64}, tab_Is[2:end, 2])[1e-3 .< true_ss .< 1e3]
     I20s = convert(Vector{Float64}, tab_Is[2:end, 3])[1e-3 .< true_ss .< 1e3]
     I40s = convert(Vector{Float64}, tab_Is[2:end, 4])[1e-3 .< true_ss .< 1e3]
     I02s = convert(Vector{Float64}, tab_Is[2:end, 5])[1e-3 .< true_ss .< 1e3] ./ ss .^ 2
     I22s = convert(Vector{Float64}, tab_Is[2:end, 6])[1e-3 .< true_ss .< 1e3] ./ ss .^ 2
     I31s = convert(Vector{Float64}, tab_Is[2:end, 7])[1e-3 .< true_ss .< 1e3] ./ ss
     I11s = convert(Vector{Float64}, tab_Is[2:end, 8])[1e-3 .< true_ss .< 1e3] ./ ss
     =#


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

     @testset "test I00" begin
          @test tools.I00.l_si ≈ -0.3480008484088528
          @test tools.I00.l_b ≈ 11.575473646894222
          @test tools.I00.l_a ≈  -6.492685627424222
          @test tools.I00.left ≈ 0.05

          @test tools.I00.r_si ≈ 2.8758516589192706
          @test tools.I00.r_b ≈ 8.345990355535678e-29
          @test tools.I00.r_a ≈ 0.0
          @test tools.I00.right ≈ 96466.16199111959

          @test all([isapprox(tools.I00(s), i, rtol = 1e-3) for (s, i) in zip(ss, I00s)])
     end

     @testset "test I20" begin
          @test tools.I20.l_si ≈ -0.3434449376790657
          @test tools.I20.l_b ≈ 1.5505826777939757
          @test tools.I20.l_a ≈  -0.035265713349630266
          @test tools.I20.left ≈ 0.05

          @test tools.I20.r_si ≈ -3.9578746738388397
          @test tools.I20.r_b ≈ 1.1832931457198192e6
          @test tools.I20.r_a ≈ 0.0
          @test tools.I20.right ≈ 96466.16199111959

          @test all([isapprox(tools.I20(s), i, rtol = 1e-3) for (s, i) in zip(ss, I20s)])
     end

     @testset "test I40" begin
          @test tools.I40.l_si ≈ -0.3468022218825587
          @test tools.I40.l_b ≈ 0.7681399976144987
          @test tools.I40.l_a ≈  -0.00028247463623223634
          @test tools.I40.left ≈ 0.05

          @test tools.I40.r_si ≈ -3.9597415622888192
          @test tools.I40.r_b ≈ 6.92151268083521e6
          @test tools.I40.r_a ≈ 0.0
          @test tools.I40.right ≈ 96466.16199111959

          @test all([isapprox(tools.I40(s), i, rtol = 1e-3) for (s, i) in zip(ss, I40s)])
     end

     @testset "test I02" begin
          @test tools.I02.l_si ≈ -2.00164755613052
          @test tools.I02.l_b ≈ 100.44877199997933
          @test tools.I02.l_a ≈  0.0
          @test tools.I02.left ≈ 0.05

          @test tools.I02.r_si ≈ -5.866385247535877
          @test tools.I02.r_b ≈ 1.7556551321550488e14
          @test tools.I02.r_a ≈ 0.0
          @test tools.I02.right ≈ 96466.16199111959

          @test all([isapprox(tools.I02(s), i, rtol = 1e-3) for (s, i) in zip(ss, I02s)])
     end

     @testset "test I22" begin
          @test tools.I22.l_si ≈ -0.34727171081772734
          @test tools.I22.l_b ≈ 0.9412777254767503
          @test tools.I22.l_a ≈ -0.43616927284862367
          @test tools.I22.left ≈ 0.05

          @test tools.I22.r_si ≈ -3.431171897505612
          @test tools.I22.r_b ≈ 851.4607831995294
          @test tools.I22.r_a ≈ 0.0
          @test tools.I22.right ≈ 96466.16199111959

          @test all([isapprox(tools.I22(s), i, rtol = 1e-3) for (s, i) in zip(ss, I22s)])
     end

     @testset "test I31" begin
          @test tools.I31.l_si ≈  -0.3445696461864395
          @test tools.I31.l_b ≈ 0.33123234468023927
          @test tools.I31.l_a ≈ -0.005064156263986402
          @test tools.I31.left ≈ 0.05

          @test tools.I31.r_si ≈ -3.9594640605240348
          @test tools.I31.r_b ≈ 1.157797082289167e6
          @test tools.I31.r_a ≈ 0.0
          @test tools.I31.right ≈ 96466.16199111959

          @test all([isapprox(tools.I31(s), i, rtol = 1e-3) for (s, i) in zip(ss, I31s)])
     end

     @testset "test I11" begin
          @test tools.I11.l_si ≈ -0.3474733407904323
          @test tools.I11.l_b ≈ 4.375198496682168
          @test tools.I11.l_a ≈ -2.17582565657694
          @test tools.I11.left ≈ 0.05

          @test tools.I11.r_si ≈ -2.270802848837033
          @test tools.I11.r_b ≈ 0.002613631263376589
          @test tools.I11.r_a ≈ 0.0
          @test tools.I11.right ≈ 96466.16199111959

          @test all([isapprox(tools.I11(s), i, rtol = 1e-3) for (s, i) in zip(ss, I11s)])
     end

     @testset "test I13" begin
          @test tools.I13.l_si ≈  -2.0001163229203924
          @test tools.I13.l_b ≈ 33.62113907833843
          @test tools.I13.l_a ≈ -0.5498036524608388
          @test tools.I13.left ≈ 0.05

          @test tools.I13.r_si ≈ -4.387393755963048
          @test tools.I13.r_b ≈ 1.5388770882263558e7
          @test tools.I13.r_a ≈ 0.0
          @test tools.I13.right ≈ 96466.16199111959

          @test all([isapprox(tools.I13(s), i, rtol = 1e-3) for (s, i) in zip(ss, I13s)])
     end

     @testset "test I04_tilde" begin
          @test tools.I04_tilde.l_si ≈ -2.0002445476746975
          @test tools.I04_tilde.l_b ≈ -16.80636730509123
          @test tools.I04_tilde.l_a ≈ 0.10615797439188279
          @test tools.I04_tilde.left ≈ 0.1

          @test tools.I04_tilde.r_si ≈ -3.746166614932489
          @test tools.I04_tilde.r_b ≈ -83772.05493095529
          @test tools.I04_tilde.r_a ≈ 0.0
          @test tools.I04_tilde.right ≈ 9888.080416307215

          @test all([isapprox(tools.I04_tilde(s), i, rtol = 1e-3) for (s, i) in zip(ss, I04_tildes)])
     end
end



@testset "V_survey" begin
     @test isapprox(GaPSE.V_survey(1, 2, π/2), 14.660765716752364, rtol=1e-4)
     @test isapprox(GaPSE.V_survey(10, 20, π/2), 14660.765716752363, rtol=1e-4)
     
     @test isapprox(GaPSE.V_survey(1, 2, 0.0), 0.0, rtol=1e-4)
     @test isapprox(GaPSE.V_survey(1, 2, π/6), 1.9641701671128426, rtol=1e-4)

     @test isapprox(GaPSE.V_survey(S_MIN, S_MAX, π/2), 3.845366169354268e8, rtol=1e-4)
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
     @test GaPSE.W(0, π/3) ≈ 1.0
     @test GaPSE.W(1.0, π/3) ≈ 1.0
     @test GaPSE.W(π/3, π/3) ≈ 0.0
     @test GaPSE.W(π, π/3) ≈ 0.0
end

