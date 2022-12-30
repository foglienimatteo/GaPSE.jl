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

#=


@testset "test ξ_from_PS" begin
     RTOL = 1e-3
     kwargs_xis = Dict(:epl => true, :pr => false,
          :N_left => 12, :N_right => 12,
          :p0_left => [-2.0, 1.0], :p0_right => [-2.0, 1.0],
          :N => 300, :int_k_min => 1e-4, :int_k_max => 1e4)

     @testset "monopole" begin
          L = 0
          input = FILE_PS
          true_xi = "datatest/XiMatter/xi_matter_L$L" * ".txt"

          table = readdlm(true_xi; comments=true)
          ss = convert(Vector{Float64}, table[:, 1])
          xis = convert(Vector{Float64}, table[:, 2])

          calc_ss, calc_xis = GaPSE.ξ_from_PS(input; L=L, kwargs_xis...)

          @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(ss, calc_ss)])
          @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(xis, calc_xis)])
     end

     @testset "quadrupole" begin
          L = 2
          input = FILE_PS
          true_xi = "datatest/XiMatter/xi_matter_L$L" * ".txt"

          table = readdlm(true_xi; comments=true)
          ss = convert(Vector{Float64}, table[:, 1])
          xis = convert(Vector{Float64}, table[:, 2])

          calc_ss, calc_xis = GaPSE.ξ_from_PS(input; L=L, kwargs_xis...)

          @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(ss, calc_ss)])
          @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(xis, calc_xis)])
     end

end


@testset "test print_ξ_from_PS first method" begin
     RTOL = 1e-3
     kwargs_xis = Dict(:epl => true, :pr => false,
          :N_left => 12, :N_right => 12,
          :p0_left => [-2.0, 1.0], :p0_right => [-2.0, 1.0],
          :N => 300, :int_k_min => 1e-4, :int_k_max => 1e4)

     @testset "monopole" begin
          L = 0
          input = FILE_PS
          true_xi = "datatest/XiMatter/xi_matter_L$L" * ".txt"
          out_file = "cacl_xi_matter_L$L" * ".txt"

          isfile(out_file) && rm(out_file)

          table = readdlm(true_xi; comments=true)
          ss = convert(Vector{Float64}, table[:, 1])
          xis = convert(Vector{Float64}, table[:, 2])

          GaPSE.print_ξ_from_PS(input, out_file; L=L, kwargs_xis...)
          calc_table = readdlm(out_file; comments=true)
          calc_ss = convert(Vector{Float64}, calc_table[:, 1])
          calc_xis = convert(Vector{Float64}, calc_table[:, 2])

          @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(ss, calc_ss)])
          @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(xis, calc_xis)])

          rm(out_file)
     end

     @testset "quadrupole" begin
          L = 2
          input = FILE_PS
          true_xi = "datatest/XiMatter/xi_matter_L$L" * ".txt"
          out_file = "cacl_xi_matter_L$L" * ".txt"

          isfile(out_file) && rm(out_file)

          table = readdlm(true_xi; comments=true)
          ss = convert(Vector{Float64}, table[:, 1])
          xis = convert(Vector{Float64}, table[:, 2])

          GaPSE.print_ξ_from_PS(input, out_file; L=L, kwargs_xis...)
          calc_table = readdlm(out_file; comments=true)
          calc_ss = convert(Vector{Float64}, calc_table[:, 1])
          calc_xis = convert(Vector{Float64}, calc_table[:, 2])

          @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(ss, calc_ss)])
          @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(xis, calc_xis)])

          rm(out_file)
     end

end

@testset "test print_ξ_from_PS second method" begin
     RTOL = 1e-3
     kwargs_xis = Dict(:epl => true, :pr => false,
          :N_left => 12, :N_right => 12,
          :p0_left => [-2.0, 1.0], :p0_right => [-2.0, 1.0],
          :N => 300, :int_k_min => 1e-4, :int_k_max => 1e4)

     @testset "monopole" begin
          L = 0
          input = FILE_PS
          true_xi = "datatest/XiMatter/xi_matter_L$L" * ".txt"
          out_file = "cacl_xi_matter_L$L" * ".txt"

          isfile(out_file) && rm(out_file)

          table = readdlm(true_xi; comments=true)
          ss = convert(Vector{Float64}, table[:, 1])
          xis = convert(Vector{Float64}, table[:, 2])


          in_table = readdlm(input; comments=true)
          in_ss = convert(Vector{Float64}, in_table[:, 1])
          in_xis = convert(Vector{Float64}, in_table[:, 2])
          GaPSE.print_ξ_from_PS(in_ss, in_xis, out_file; L=L, kwargs_xis...)
          calc_table = readdlm(out_file; comments=true)
          calc_ss = convert(Vector{Float64}, calc_table[:, 1])
          calc_xis = convert(Vector{Float64}, calc_table[:, 2])

          @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(ss, calc_ss)])
          @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(xis, calc_xis)])

          rm(out_file)
     end

     @testset "quadrupole" begin
          L = 2
          input = FILE_PS
          true_xi = "datatest/XiMatter/xi_matter_L$L" * ".txt"
          out_file = "cacl_xi_matter_L$L" * ".txt"

          isfile(out_file) && rm(out_file)

          table = readdlm(true_xi; comments=true)
          ss = convert(Vector{Float64}, table[:, 1])
          xis = convert(Vector{Float64}, table[:, 2])

          in_table = readdlm(input; comments=true)
          in_ss = convert(Vector{Float64}, in_table[:, 1])
          in_xis = convert(Vector{Float64}, in_table[:, 2])
          GaPSE.print_ξ_from_PS(in_ss, in_xis, out_file; L=L, kwargs_xis...)
          calc_table = readdlm(out_file; comments=true)
          calc_ss = convert(Vector{Float64}, calc_table[:, 1])
          calc_xis = convert(Vector{Float64}, calc_table[:, 2])

          @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(ss, calc_ss)])
          @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(xis, calc_xis)])

          rm(out_file)
     end

end

=#