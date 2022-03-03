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


@testset "test integral_on_mu" begin
     name_effect = "auto_doppler"
     func_effect = GaPSE.ξ_Doppler
     RTOL = 1e-2

     kwargs = Dict(
          :enhancer => 1e8, :N_μs => 30,
          :μ_atol => 0.0, :μ_rtol => 1e-2,
     )

     ss = [10, 100, 300, 700, 1000]

     @testset "zeros" begin
          a_func(x) = x
          @test_throws AssertionError GaPSE.integral_on_mu(COSMO.s_eff, 10, "strange", COSMO;
               L = 0, use_windows = false, SPLINE = true, kwargs...)
          @test_throws AssertionError GaPSE.integral_on_mu(COSMO.s_eff, 10, a_func, COSMO;
               L = 0, use_windows = false, SPLINE = true, kwargs...)
     end

     @testset "L = 0, no_window, no_spline" begin
          L, use_windows, SPLINE = 0, false, false
          res = [6.216729464626079e-5, 1.204973036205085e-5, 1.7985341392845314e-6,
               2.811035752494547e-7, 7.40858576296966e-8]

          calc_ints_1 = [GaPSE.integral_on_mu(COSMO.s_eff, s, name_effect, COSMO;
               L = L, use_windows = use_windows, SPLINE = SPLINE, kwargs...) for s in ss]
          calc_ints_2 = [GaPSE.integral_on_mu(COSMO.s_eff, s, func_effect, COSMO;
               L = L, use_windows = use_windows, SPLINE = SPLINE, kwargs...) for s in ss]

          @test all([isapprox(a, r, rtol = RTOL) for (a, r) in zip(calc_ints_1, res)])
          @test all([isapprox(a, r, rtol = RTOL) for (a, r) in zip(calc_ints_2, res)])
     end

     @testset "L = 1, no_window, no_spline" begin
          L, use_windows, SPLINE = 1, false, false
          res = [-5.311205477414576e-7, -1.159934989165109e-6, -6.334076305557245e-7,
               -2.0381833943775459e-7, -6.206602978481535e-8]

          calc_ints_1 = [GaPSE.integral_on_mu(COSMO.s_eff, s, name_effect, COSMO;
               L = L, use_windows = use_windows, SPLINE = SPLINE, kwargs...) for s in ss]
          calc_ints_2 = [GaPSE.integral_on_mu(COSMO.s_eff, s, func_effect, COSMO;
               L = L, use_windows = use_windows, SPLINE = SPLINE, kwargs...) for s in ss]

          @test all([isapprox(a, r, rtol = RTOL) for (a, r) in zip(calc_ints_1, res)])
          @test all([isapprox(a, r, rtol = RTOL) for (a, r) in zip(calc_ints_2, res)])
     end

     @testset "L = 0, yes_window, no_spline" begin
          L, use_windows, SPLINE = 0, true, false
          res = [0.002406339467257516, 0.0004378931763537511, 5.214481768705607e-5,
               2.2114227295248647e-6, 1.8237097173450833e-9]

          calc_ints_1 = [GaPSE.integral_on_mu(COSMO.s_eff, s, name_effect, COSMO;
               L = L, use_windows = use_windows, SPLINE = SPLINE, kwargs...) for s in ss]
          calc_ints_2 = [GaPSE.integral_on_mu(COSMO.s_eff, s, func_effect, COSMO;
               L = L, use_windows = use_windows, SPLINE = SPLINE, kwargs...) for s in ss]

          @test all([isapprox(a, r, rtol = RTOL) for (a, r) in zip(calc_ints_1, res)])
          @test all([isapprox(a, r, rtol = RTOL) for (a, r) in zip(calc_ints_2, res)])
     end

     @testset "L = 1, yes_window, no_spline" begin
          L, use_windows, SPLINE = 1, true, false
          res = [-2.02440383383452e-5, -4.1180157984924956e-5, -1.8920053056318384e-5,
               -1.7050686653742866e-6, -1.8124991003555014e-9]

          calc_ints_1 = [GaPSE.integral_on_mu(COSMO.s_eff, s, name_effect, COSMO;
               L = L, use_windows = use_windows, SPLINE = SPLINE, kwargs...) for s in ss]
          calc_ints_2 = [GaPSE.integral_on_mu(COSMO.s_eff, s, func_effect, COSMO;
               L = L, use_windows = use_windows, SPLINE = SPLINE, kwargs...) for s in ss]

          @test all([isapprox(a, r, rtol = RTOL) for (a, r) in zip(calc_ints_1, res)])
          @test all([isapprox(a, r, rtol = RTOL) for (a, r) in zip(calc_ints_2, res)])
     end

     @testset "L = 0, no_window, yes_spline" begin
          L, use_windows, SPLINE = 0, false, true
          res = [6.216729467077863e-5, 1.2049734975168255e-5, 1.7982097197486115e-6,
               2.811064352363306e-7, 7.408582429145476e-8]

          calc_ints_1 = [GaPSE.integral_on_mu(COSMO.s_eff, s, name_effect, COSMO;
               L = L, use_windows = use_windows, SPLINE = SPLINE, kwargs...) for s in ss]
          calc_ints_2 = [GaPSE.integral_on_mu(COSMO.s_eff, s, func_effect, COSMO;
               L = L, use_windows = use_windows, SPLINE = SPLINE, kwargs...) for s in ss]

          @test all([isapprox(a, r, rtol = RTOL) for (a, r) in zip(calc_ints_1, res)])
          @test all([isapprox(a, r, rtol = RTOL) for (a, r) in zip(calc_ints_2, res)])
     end
end




##########################################################################################92



@testset "test ξ_multipole" begin
     name_effect = "auto_doppler"
     func_effect = GaPSE.ξ_Doppler
     RTOL = 1e-2

     kwargs = Dict(
          :use_windows => false,
          :enhancer => 1e8, :N_μs => 30,
          :μ_atol => 0.0, :μ_rtol => 1e-2,
     )


     @testset "zeros" begin
          @test_throws AssertionError GaPSE.ξ_multipole(COSMO.s_eff, 10.0, "strange", COSMO;
               L = 0, kwargs...)
     end

     @testset "L = 0" begin
          L = 0
          table = readdlm("datatest/doppler_multipoles/xi_" * name_effect * "_L$L" * ".txt"; comments = true)
          ss = convert(Vector{Float64}, table[:, 1])
          xis = convert(Vector{Float64}, table[:, 2])

          calc_xis_1 = [GaPSE.ξ_multipole(COSMO.s_eff, s, name_effect, COSMO;
               L = L, kwargs...) for s in ss]
          calc_xis_2 = [GaPSE.ξ_multipole(COSMO.s_eff, s, func_effect, COSMO;
               L = L, kwargs...) for s in ss]

          @test all([isapprox(xi, calc_xi, rtol = RTOL) for (xi, calc_xi) in zip(xis, calc_xis_1)])
          @test all([isapprox(xi, calc_xi, rtol = RTOL) for (xi, calc_xi) in zip(xis, calc_xis_2)])
     end

     @testset "L = 1" begin
          L = 1
          table = readdlm("datatest/doppler_multipoles/xi_" * name_effect * "_L$L" * ".txt"; comments = true)
          ss = convert(Vector{Float64}, table[:, 1])
          xis = convert(Vector{Float64}, table[:, 2])

          calc_xis_1 = [GaPSE.ξ_multipole(COSMO.s_eff, s, name_effect, COSMO;
               L = L, kwargs...) for s in ss]
          calc_xis_2 = [GaPSE.ξ_multipole(COSMO.s_eff, s, func_effect, COSMO;
               L = L, kwargs...) for s in ss]

          @test all([isapprox(xi, calc_xi, rtol = RTOL) for (xi, calc_xi) in zip(xis, calc_xis_1)])
          @test all([isapprox(xi, calc_xi, rtol = RTOL) for (xi, calc_xi) in zip(xis, calc_xis_2)])
     end

     @testset "L = 2" begin
          L = 2
          table = readdlm("datatest/doppler_multipoles/xi_" * name_effect * "_L$L" * ".txt"; comments = true)
          ss = convert(Vector{Float64}, table[:, 1])
          xis = convert(Vector{Float64}, table[:, 2])

          calc_xis_1 = [GaPSE.ξ_multipole(COSMO.s_eff, s, name_effect, COSMO;
               L = L, kwargs...) for s in ss]
          calc_xis_2 = [GaPSE.ξ_multipole(COSMO.s_eff, s, func_effect, COSMO;
               L = L, kwargs...) for s in ss]

          @test all([isapprox(xi, calc_xi, rtol = RTOL) for (xi, calc_xi) in zip(xis, calc_xis_1)])
          @test all([isapprox(xi, calc_xi, rtol = RTOL) for (xi, calc_xi) in zip(xis, calc_xis_2)])
     end

     @testset "L = 3" begin
          L = 3
          table = readdlm("datatest/doppler_multipoles/xi_" * name_effect * "_L$L" * ".txt"; comments = true)
          ss = convert(Vector{Float64}, table[:, 1])
          xis = convert(Vector{Float64}, table[:, 2])

          calc_xis_1 = [GaPSE.ξ_multipole(COSMO.s_eff, s, name_effect, COSMO;
               L = L, kwargs...) for s in ss]
          calc_xis_2 = [GaPSE.ξ_multipole(COSMO.s_eff, s, func_effect, COSMO;
               L = L, kwargs...) for s in ss]

          @test all([isapprox(xi, calc_xi, rtol = RTOL) for (xi, calc_xi) in zip(xis, calc_xis_1)])
          @test all([isapprox(xi, calc_xi, rtol = RTOL) for (xi, calc_xi) in zip(xis, calc_xis_2)])
     end

     @testset "L = 4" begin
          L = 4
          table = readdlm("datatest/doppler_multipoles/xi_" * name_effect * "_L$L" * ".txt"; comments = true)
          ss = convert(Vector{Float64}, table[:, 1])
          xis = convert(Vector{Float64}, table[:, 2])

          calc_xis_1 = [GaPSE.ξ_multipole(COSMO.s_eff, s, name_effect, COSMO;
               L = L, kwargs...) for s in ss]
          calc_xis_2 = [GaPSE.ξ_multipole(COSMO.s_eff, s, func_effect, COSMO;
               L = L, kwargs...) for s in ss]

          @test all([isapprox(xi, calc_xi, rtol = RTOL) for (xi, calc_xi) in zip(xis, calc_xis_1)])
          @test all([isapprox(xi, calc_xi, rtol = RTOL) for (xi, calc_xi) in zip(xis, calc_xis_2)])
     end
end



##########################################################################################92




@testset "test print_map_int_on_mu" begin
     effect = "auto_doppler"

     kwargs = Dict(
          :use_windows => false,
          :pr => false,
          :enhancer => 1e8, :N_μs => 30,
          :μ_atol => 0.0, :μ_rtol => 1e-2,
          :N_log => 1000
     )

     @testset "first" begin
          table = readdlm("datatest/int_on_mu/doppler_map_int_on_mu_first.txt"; comments = true)
          ss = convert(Vector{Float64}, table[:, 1])
          xis = convert(Vector{Float64}, table[:, 2])

          name = "calc_doppler_int_on_mu_first.txt"
          isfile(name) && rm(name)
          GaPSE.print_map_int_on_mu(COSMO, name, effect, nothing;
               s_1 = nothing, L = 0, kwargs...)

          calc_table = readdlm(name; comments = true)
          calc_ss = convert(Vector{Float64}, calc_table[:, 1])
          calc_xis = convert(Vector{Float64}, calc_table[:, 2])

          @test all([isapprox(s, calc_s, rtol = 1e-2) for (s, calc_s) in zip(ss, calc_ss)])
          @test all([isapprox(xi, calc_xi, rtol = 1e-2) for (xi, calc_xi) in zip(xis, calc_xis)])

          rm(name)
     end

     @testset "second" begin
          table = readdlm("datatest/int_on_mu/doppler_map_int_on_mu_second.txt"; comments = true)
          ss = convert(Vector{Float64}, table[:, 1])
          xis = convert(Vector{Float64}, table[:, 2])

          name = "calc_doppler_int_on_mu_second.txt"
          isfile(name) && rm(name)
          GaPSE.print_map_int_on_mu(COSMO, name, effect, 10 .^ range(0, 3, length = 344);
               s_1 = nothing, L = 0, kwargs...)

          calc_table = readdlm(name; comments = true)
          calc_ss = convert(Vector{Float64}, calc_table[:, 1])
          calc_xis = convert(Vector{Float64}, calc_table[:, 2])

          @test all([isapprox(s, calc_s, rtol = 1e-2) for (s, calc_s) in zip(ss, calc_ss)])
          @test all([isapprox(xi, calc_xi, rtol = 1e-2) for (xi, calc_xi) in zip(xis, calc_xis)])

          rm(name)
     end

     @testset "third" begin
          table = readdlm("datatest/int_on_mu/doppler_map_int_on_mu_third.txt"; comments = true)
          ss = convert(Vector{Float64}, table[:, 1])
          xis = convert(Vector{Float64}, table[:, 2])

          name = "calc_doppler_int_on_mu_third.txt"
          isfile(name) && rm(name)
          GaPSE.print_map_int_on_mu(COSMO, name, effect, 10 .^ range(0, 3, length = 344);
               s_1 = COSMO.s_eff - 65.0, L = 0, kwargs...)

          calc_table = readdlm(name; comments = true)
          calc_ss = convert(Vector{Float64}, calc_table[:, 1])
          calc_xis = convert(Vector{Float64}, calc_table[:, 2])

          @test all([isapprox(s, calc_s, rtol = 1e-2) for (s, calc_s) in zip(ss, calc_ss)])
          @test all([isapprox(xi, calc_xi, rtol = 1e-2) for (xi, calc_xi) in zip(xis, calc_xis)])

          rm(name)
     end
end

