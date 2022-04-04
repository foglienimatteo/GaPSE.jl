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

@testset "test PS_multipole, first method" begin
     RTOL = 1e-3

     @testset "with F" begin
          #=
          # This xi was obtained with the following command:
          GaPSE.print_map_ξ_multipole(cosmo, "xi_auto_doppler_withF_L0.txt", 
               "auto_doppler"; use_windows = true, N_log = 1000, N_μs = 50, pr = false)
          =#
          input = "datatest/power_spectrum/xi_auto_doppler_withF_L0.txt"
          true_pk = "datatest/power_spectrum/ps_auto_doppler_withF_L0.txt"

          table = readdlm(true_pk; comments=true)
          ks = convert(Vector{Float64}, table[:, 1])
          pks = convert(Vector{Float64}, table[:, 2])

          calc_ks, calc_pks = GaPSE.PS_multipole(input; N_left=12, N_right=12,
               L=0, N=1000, pr=false, int_s_min=1e-1, int_s_max=1e3)

          @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(ks, calc_ks)])
          @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(pks, calc_pks)])
     end

     @testset "without F" begin
          #=
          # This xi was obtained with the following command:
          GaPSE.print_map_ξ_multipole(cosmo, "xi_auto_doppler_noF_L0.txt", 
               "auto_doppler"; use_windows = false, N_log = 1000, N_μs = 50, pr = false)
          =#

          input = "datatest/power_spectrum/xi_auto_doppler_noF_L0.txt"
          true_pk = "datatest/power_spectrum/ps_auto_doppler_noF_L0.txt"

          table = readdlm(true_pk; comments=true)
          ks = convert(Vector{Float64}, table[:, 1])
          pks = convert(Vector{Float64}, table[:, 2])

          calc_ks, calc_pks = GaPSE.PS_multipole(input; N_left=12, N_right=12,
               L=0, N=1000, pr=false, int_s_min=1e-1, int_s_max=1e3)

          @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(ks, calc_ks)])
          @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(pks, calc_pks)])
     end
end


@testset "test PS_multipole, second method" begin
     RTOL = 1e-3

     @testset "with F" begin
          #=
          # This xi was obtained with the following command:
          GaPSE.print_map_ξ_multipole(cosmo, "xi_auto_doppler_withF_L0.txt", 
               "auto_doppler"; use_windows = true, N_log = 1000, N_μs = 50, pr = false)
          =#
          input = "datatest/power_spectrum/xi_auto_doppler_withF_L0.txt"
          true_pk = "datatest/power_spectrum/ps_auto_doppler_withF_L0.txt"
          out_file = "calc_pk_auto_doppler_withF_L0.txt"

          isfile(out_file) && rm(out_file)

          table = readdlm(true_pk; comments=true)
          ks = convert(Vector{Float64}, table[:, 1])
          pks = convert(Vector{Float64}, table[:, 2])

          GaPSE.print_PS_multipole(input, out_file; N_left=12, N_right=12,
               L=0, N=1000, pr=false, int_s_min=1e-1, int_s_max=1e3)
          calc_table = readdlm(out_file; comments=true)
          calc_ks = convert(Vector{Float64}, calc_table[:, 1])
          calc_pks = convert(Vector{Float64}, calc_table[:, 2])

          @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(ks, calc_ks)])
          @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(pks, calc_pks)])

          rm(out_file)
     end

     @testset "without F" begin
          #=
          # This xi was obtained with the following command:
          GaPSE.print_map_ξ_multipole(cosmo, "xi_auto_doppler_noF_L0.txt", 
               "auto_doppler"; use_windows = false, N_log = 1000, N_μs = 50, pr = false)
          =#

          input = "datatest/power_spectrum/xi_auto_doppler_noF_L0.txt"
          true_pk = "datatest/power_spectrum/ps_auto_doppler_noF_L0.txt"
          out_file = "calc_pk_auto_doppler_noF_L0.txt"

          isfile(out_file) && rm(out_file)

          table = readdlm(true_pk; comments=true)
          ks = convert(Vector{Float64}, table[:, 1])
          pks = convert(Vector{Float64}, table[:, 2])

          GaPSE.print_PS_multipole(input, out_file; N_left=12, N_right=12,
               L=0, N=1000, pr=false, int_s_min=1e-1, int_s_max=1e3)
          calc_table = readdlm(out_file; comments=true)
          calc_ks = convert(Vector{Float64}, calc_table[:, 1])
          calc_pks = convert(Vector{Float64}, calc_table[:, 2])

          @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(ks, calc_ks)])
          @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(pks, calc_pks)])

          rm(out_file)
     end
end


@testset "test PS_multipole, third method" begin
     RTOL = 1e-3


     @testset "with F" begin
          #=
          # This xi was obtained with the following command:
          GaPSE.print_map_ξ_multipole(cosmo, "xi_auto_doppler_withF_L0.txt", 
               "auto_doppler"; use_windows = true, N_log = 1000, N_μs = 50, pr = false)
          =#
          input = "datatest/power_spectrum/xi_auto_doppler_withF_L0.txt"
          true_pk = "datatest/power_spectrum/ps_auto_doppler_withF_L0.txt"
          out_file = "calc_pk_auto_doppler_withF_L0.txt"

          isfile(out_file) && rm(out_file)

          table = readdlm(true_pk; comments=true)
          ks = convert(Vector{Float64}, table[:, 1])
          pks = convert(Vector{Float64}, table[:, 2])

          xi_table = readdlm(input; comments=true)
          ss = convert(Vector{Float64}, xi_table[:, 1])
          fs = convert(Vector{Float64}, xi_table[:, 2])

          GaPSE.print_PS_multipole(ss, fs, out_file; N_left=12, N_right=12,
               L=0, N=1000, pr=false, int_s_min=1e-1, int_s_max=1e3)
          calc_table = readdlm(out_file; comments=true)
          calc_ks = convert(Vector{Float64}, calc_table[:, 1])
          calc_pks = convert(Vector{Float64}, calc_table[:, 2])

          @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(ks, calc_ks)])
          @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(pks, calc_pks)])

          rm(out_file)
     end

     @testset "without F" begin
          #=
          # This xi was obtained with the following command:
          GaPSE.print_map_ξ_multipole(cosmo, "xi_auto_doppler_noF_L0.txt", 
               "auto_doppler"; use_windows = false, N_log = 1000, N_μs = 50, pr = false)
          =#

          input = "datatest/power_spectrum/xi_auto_doppler_noF_L0.txt"
          true_pk = "datatest/power_spectrum/ps_auto_doppler_noF_L0.txt"
          out_file = "calc_pk_auto_doppler_noF_L0.txt"

          isfile(out_file) && rm(out_file)

          table = readdlm(true_pk; comments=true)
          ks = convert(Vector{Float64}, table[:, 1])
          pks = convert(Vector{Float64}, table[:, 2])

          GaPSE.print_PS_multipole(input, out_file; N_left=12, N_right=12,
               L=0, N=1000, pr=false, int_s_min=1e-1, int_s_max=1e3)
          calc_table = readdlm(out_file; comments=true)
          calc_ks = convert(Vector{Float64}, calc_table[:, 1])
          calc_pks = convert(Vector{Float64}, calc_table[:, 2])

          @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(ks, calc_ks)])
          @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(pks, calc_pks)])

          rm(out_file)
     end
end


##########################################################################################92

#=
@testset "test print_all_PS_multipole, first method" begin
     RTOL = 1e-3

     #=
     # This xi was obtained with the following command:
     GaPSE.print_map_ξ_multipole(cosmo, "xi_auto_doppler_withF_L0.txt", 
          "auto_doppler"; use_windows = true, N_log = 1000, N_μs = 50, pr = false)
     =#
     input = "datatest/power_spectrum/xi_auto_doppler_withF_L0.txt"
     true_pk = "datatest/power_spectrum/ps_auto_doppler_withF_L0.txt"
     out_file = "calc_pk_auto_doppler_withF_L0.txt"

     isfile(out_file) && rm(out_file)

     table = readdlm(true_pk; comments=true)
     ks = convert(Vector{Float64}, table[:, 1])
     pks = convert(Vector{Float64}, table[:, 2])

     GaPSE.print_PS_multipole(input, out_file; N_left=12, N_right=12,
          L=0, N=1000, pr=false, int_s_min=1e-1, int_s_max=1e3)
     calc_table = readdlm(out_file; comments=true)
     calc_ks = convert(Vector{Float64}, calc_table[:, 1])
     calc_pks = convert(Vector{Float64}, calc_table[:, 2])

     @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(ks, calc_ks)])
     @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(pks, calc_pks)])

     rm(out_file)
    
end
=#
