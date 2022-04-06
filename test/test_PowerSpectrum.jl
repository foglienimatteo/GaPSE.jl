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
# The TPCF files in this tests were obtained with the following code lines: 
kwargs_xis = Dict(
     :pr => false,
     :enhancer => 1e8, :N_μs => 30,
     :μ_atol => 0.0, :μ_rtol => 1e-2,
     :N_log => 100, :s1 => nothing,
);
common = PATH_TO_GAPSE * "test/datatest/power_spectrum/"

GaPSE.print_map_ξ_multipole(cosmo,  common*"xi_auto_doppler_withF_L0.txt", 
    "auto_doppler", 10 .^ range(0,3,length=300); use_windows = false, 
    L = 0, kwargs_xis...);
GaPSE.print_map_ξ_multipole(cosmo,  common*"xi_auto_doppler_withF_L0.txt", 
    "auto_doppler", 10 .^ range(0,3,length=300); use_windows = false, 
    L = 2, kwargs_xis...);

GaPSE.print_map_ξ_multipole(cosmo,  common*"xi_auto_doppler_withF_L0.txt", 
    "auto_doppler", 10 .^ range(0,3,length=300); use_windows = true, 
    L = 0, kwargs_xis...);
GaPSE.print_map_ξ_multipole(cosmo,  common*"xi_auto_doppler_withF_L2.txt", 
    "auto_doppler", 10 .^ range(0,3,length=300); use_windows = true, 
    L = 2, kwargs_xis...);
=#

#=
# The Power Spectrum files in this tests were obtained with the following code lines: 

kwargs = Dict(:epl=>true, :pr=>false,
     :N_left=>12, :N_right=>12,
     :p0_left=>[-2.0, 1.0], :p0_right=>[-2.0, 1.0],
    :N => 300, :int_s_min => 1e-4, :int_s_max => 1e4)
common = PATH_TO_GAPSE * "test/datatest/power_spectrum/"

GaPSE.print_PS_multipole(ss, xis_L0, common*"ps_auto_doppler_noF_L0.txt"; L=0, kwargs...);
GaPSE.print_PS_multipole(ss, xis_L2, common*"ps_auto_doppler_noF_L2.txt"; L=2, kwargs...);

GaPSE.print_PS_multipole(ss, xis_L0_wF, common*"ps_auto_doppler_withF_L0.txt"; L=0, kwargs...);
GaPSE.print_PS_multipole(ss, xis_L2_wF, common*"ps_auto_doppler_withF_L2.txt"; L=2, kwargs...);
=#



@testset "test PS_multipole" begin
     RTOL = 1e-3
     kwargs_ps = Dict(:epl => true, :pr => false,
          :N_left => 12, :N_right => 12,
          :p0_left => [-2.0, 1.0], :p0_right => [-2.0, 1.0],
          :N => 300, :int_s_min => 1e-4, :int_s_max => 1e4)

     @testset "with F" begin
          @testset "monopole" begin
               L = 0
               input = "datatest/power_spectrum/xi_auto_doppler_withF_L$L" * ".txt"
               true_pk = "datatest/power_spectrum/ps_auto_doppler_withF_L$L" * ".txt"

               table = readdlm(true_pk; comments=true)
               ks = convert(Vector{Float64}, table[:, 1])
               pks = convert(Vector{Float64}, table[:, 2])

               calc_ks, calc_pks = GaPSE.PS_multipole(input; L=L, kwargs_ps...)

               println("calc_ks = $calc_ks ;")
               println("calc_pks = $calc_pks ;")
               println("pks = $pks ;")

               @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(ks, calc_ks)])
               @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(pks, calc_pks)])
          end

          @testset "quadrupole" begin
               L = 2
               input = "datatest/power_spectrum/xi_auto_doppler_withF_L$L" * ".txt"
               true_pk = "datatest/power_spectrum/ps_auto_doppler_withF_L$L" * ".txt"

               table = readdlm(true_pk; comments=true)
               ks = convert(Vector{Float64}, table[:, 1])
               pks = convert(Vector{Float64}, table[:, 2])

               calc_ks, calc_pks = GaPSE.PS_multipole(input; L=L, kwargs_ps...)

               @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(ks, calc_ks)])
               @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(pks, calc_pks)])
          end
     end

     @testset "without F" begin
          @testset "monopole" begin
               L = 0
               input = "datatest/power_spectrum/xi_auto_doppler_noF_L$L" * ".txt"
               true_pk = "datatest/power_spectrum/ps_auto_doppler_noF_L$L" * ".txt"

               table = readdlm(true_pk; comments=true)
               ks = convert(Vector{Float64}, table[:, 1])
               pks = convert(Vector{Float64}, table[:, 2])

               calc_ks, calc_pks = GaPSE.PS_multipole(input; L=L, kwargs_ps...)

               @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(ks, calc_ks)])
               @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(pks, calc_pks)])
          end

          @testset "quadrupole" begin
               L = 2
               input = "datatest/power_spectrum/xi_auto_doppler_noF_L$L" * ".txt"
               true_pk = "datatest/power_spectrum/ps_auto_doppler_noF_L$L" * ".txt"

               table = readdlm(true_pk; comments=true)
               ks = convert(Vector{Float64}, table[:, 1])
               pks = convert(Vector{Float64}, table[:, 2])

               calc_ks, calc_pks = GaPSE.PS_multipole(input; L=L, kwargs_ps...)

               @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(ks, calc_ks)])
               @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(pks, calc_pks)])
          end
     end
end


@testset "test print_PS_multipole" begin
     RTOL = 1e-3
     kwargs_ps = Dict(:epl => true, :pr => false,
          :N_left => 12, :N_right => 12,
          :p0_left => [-2.0, 1.0], :p0_right => [-2.0, 1.0],
          :N => 300, :int_s_min => 1e-4, :int_s_max => 1e4)

     @testset "with F" begin
          @testset "monopole" begin
               L = 0
               input = "datatest/power_spectrum/xi_auto_doppler_withF_L$L" * ".txt"
               true_pk = "datatest/power_spectrum/ps_auto_doppler_withF_L$L" * ".txt"
               out_file = "calc_pk_auto_doppler_withF_L$L" * ".txt"

               isfile(out_file) && rm(out_file)

               table = readdlm(true_pk; comments=true)
               ks = convert(Vector{Float64}, table[:, 1])
               pks = convert(Vector{Float64}, table[:, 2])

               GaPSE.print_PS_multipole(input, out_file; L=L, kwargs_ps...)
               calc_table = readdlm(out_file; comments=true)
               calc_ks = convert(Vector{Float64}, calc_table[:, 1])
               calc_pks = convert(Vector{Float64}, calc_table[:, 2])

               @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(ks, calc_ks)])
               @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(pks, calc_pks)])

               rm(out_file)
          end

          @testset "quadrupole" begin
               L = 2
               input = "datatest/power_spectrum/xi_auto_doppler_withF_L$L" * ".txt"
               true_pk = "datatest/power_spectrum/ps_auto_doppler_withF_L$L" * ".txt"
               out_file = "calc_pk_auto_doppler_withF_L$L" * ".txt"

               isfile(out_file) && rm(out_file)

               table = readdlm(true_pk; comments=true)
               ks = convert(Vector{Float64}, table[:, 1])
               pks = convert(Vector{Float64}, table[:, 2])

               GaPSE.print_PS_multipole(input, out_file; L=L, kwargs_ps...)
               calc_table = readdlm(out_file; comments=true)
               calc_ks = convert(Vector{Float64}, calc_table[:, 1])
               calc_pks = convert(Vector{Float64}, calc_table[:, 2])

               @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(ks, calc_ks)])
               @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(pks, calc_pks)])

               rm(out_file)
          end
     end

     @testset "without F" begin
          @testset "monopole" begin
               L = 0
               input = "datatest/power_spectrum/xi_auto_doppler_noF_L$L" * ".txt"
               true_pk = "datatest/power_spectrum/ps_auto_doppler_noF_L$L" * ".txt"
               out_file = "calc_pk_auto_doppler_noF_L$L" * ".txt"

               isfile(out_file) && rm(out_file)

               table = readdlm(true_pk; comments=true)
               ks = convert(Vector{Float64}, table[:, 1])
               pks = convert(Vector{Float64}, table[:, 2])

               GaPSE.print_PS_multipole(input, out_file; L=L, kwargs_ps...)
               calc_table = readdlm(out_file; comments=true)
               calc_ks = convert(Vector{Float64}, calc_table[:, 1])
               calc_pks = convert(Vector{Float64}, calc_table[:, 2])

               @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(ks, calc_ks)])
               @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(pks, calc_pks)])

               rm(out_file)
          end

          @testset "quadrupole" begin
               L = 2
               input = "datatest/power_spectrum/xi_auto_doppler_noF_L$L" * ".txt"
               true_pk = "datatest/power_spectrum/ps_auto_doppler_noF_L$L" * ".txt"
               out_file = "calc_pk_auto_doppler_noF_L$L" * ".txt"

               isfile(out_file) && rm(out_file)
               "datatest/power_spectrum/x"

               table = readdlm(true_pk; comments=true)
               ks = convert(Vector{Float64}, table[:, 1])
               pks = convert(Vector{Float64}, table[:, 2])

               GaPSE.print_PS_multipole(input, out_file; L=L, kwargs_ps...)
               calc_table = readdlm(out_file; comments=true)
               calc_ks = convert(Vector{Float64}, calc_table[:, 1])
               calc_pks = convert(Vector{Float64}, calc_table[:, 2])

               @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(ks, calc_ks)])
               @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(pks, calc_pks)])

               rm(out_file)
          end
     end
end


