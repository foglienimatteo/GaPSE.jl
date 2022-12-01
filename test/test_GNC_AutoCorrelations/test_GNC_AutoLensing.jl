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

@testset "test xi auto_lensing L = 0 no observer terms" begin
     effect = "auto_lensing"
     L = 0

     ss = ss_GNC_L0_noF_noobs
     xis = all_xis_GNC_L0_noF_noobs[GaPSE.INDEX_GR_EFFECT_GNC[effect]]

     name = "calc_xi_" * effect * "_GNC_L$L" * ".txt"
     isfile(name) && rm(name)
     GaPSE.print_map_ξ_GNC_multipole(COSMO, name, effect, SS_GNC;
          L = L, alg = :quad, obs = :no, GaPSE.specif_kwargs_GNC(effect, KWARGS_GNC)...)

     calc_table = readdlm(name; comments = true)
     calc_ss = convert(Vector{Float64}, calc_table[:, 1])
     calc_xis = convert(Vector{Float64}, calc_table[:, 2])

     @test all([isapprox(s, calc_s, rtol = 1e-2) for (s, calc_s) in zip(ss, calc_ss)])
     @test all([isapprox(xi, calc_xi, rtol = 1e-2) for (xi, calc_xi) in zip(xis, calc_xis)])

     rm(name)

     #println("xis = $xis ;")
     #println("calc_xis = $calc_xis ;")
end

@testset "test xi auto_lensing L = 0 no observer velocity terms" begin
     effect = "auto_lensing"
     L = 0

     ss = ss_GNC_L0_noF_noobsvel
     xis = all_xis_GNC_L0_noF_noobsvel[GaPSE.INDEX_GR_EFFECT_GNC[effect]]

     name = "calc_xi_" * effect * "_GNC_L$L" * ".txt"
     isfile(name) && rm(name)
     GaPSE.print_map_ξ_GNC_multipole(COSMO, name, effect, SS_GNC;
          L = L, alg = :quad, obs = :noobsvel, GaPSE.specif_kwargs_GNC(effect, KWARGS_GNC)...)

     calc_table = readdlm(name; comments = true)
     calc_ss = convert(Vector{Float64}, calc_table[:, 1])
     calc_xis = convert(Vector{Float64}, calc_table[:, 2])

     @test all([isapprox(s, calc_s, rtol = 1e-2) for (s, calc_s) in zip(ss, calc_ss)])
     @test all([isapprox(xi, calc_xi, rtol = 1e-2) for (xi, calc_xi) in zip(xis, calc_xis)])

     rm(name)

     #println("xis = $xis ;")
     #println("calc_xis = $calc_xis ;")
end


@testset "test xi auto_lensing L = 0 with observer terms" begin
     effect = "auto_lensing"
     L = 0

     ss = ss_GNC_L0_noF_withobs
     xis = all_xis_GNC_L0_noF_withobs[GaPSE.INDEX_GR_EFFECT_GNC[effect]]

     name = "calc_xi_" * effect * "_GNC_L$L" * ".txt"
     isfile(name) && rm(name)
     GaPSE.print_map_ξ_GNC_multipole(COSMO, name, effect, SS_GNC;
          L = L, alg = :quad, obs = :yes, GaPSE.specif_kwargs_GNC(effect, KWARGS_GNC)...)

     calc_table = readdlm(name; comments = true)
     calc_ss = convert(Vector{Float64}, calc_table[:, 1])
     calc_xis = convert(Vector{Float64}, calc_table[:, 2])

     @test all([isapprox(s, calc_s, rtol = 1e-2) for (s, calc_s) in zip(ss, calc_ss)])
     @test all([isapprox(xi, calc_xi, rtol = 1e-2) for (xi, calc_xi) in zip(xis, calc_xis)])

     rm(name)

     #println("xis = $xis ;")
     #println("calc_xis = $calc_xis ;")
end
