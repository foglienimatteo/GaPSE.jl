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

@testset "test xi GNCxLD doppler_lensing L = 0" begin
     effect = "doppler_lensing"
     L = 0

     ss = ss_GNCxLD_L0_noF
     xis = all_xis_GNCxLD_L0_noF[GaPSE.INDEX_GR_EFFECT_GNCxLD[effect]]

     name = "calc_xi_" * effect * "_GNCxLD_L$L" * ".txt"
     isfile(name) && rm(name)
     GaPSE.print_map_ξ_GNCxLD_multipole(COSMO, name, effect, 10 .^ range(-1, 3, length=KWARGS_GNCxLD[:N_log]);
          L = L, GaPSE.specif_kwargs_GNCxLD(effect, KWARGS_GNCxLD)...)

     calc_table = readdlm(name; comments = true)
     calc_ss = convert(Vector{Float64}, calc_table[:, 1])
     calc_xis = convert(Vector{Float64}, calc_table[:, 2])

     @test all([isapprox(s, calc_s, rtol = 1e-2) for (s, calc_s) in zip(ss, calc_ss)])
     @test all([isapprox(xi, calc_xi, rtol = 1e-2) for (xi, calc_xi) in zip(xis, calc_xis)])

     rm(name)
end

@testset "test xi LDxGNC lensing_doppler L = 0" begin
     effect = "lensing_doppler"
     L = 0

     ss = ss_LDxGNC_L0_noF
     xis = all_xis_LDxGNC_L0_noF[GaPSE.INDEX_GR_EFFECT_LDxGNC[effect]]

     name = "calc_xi_" * effect * "_LDxGNC_L$L" * ".txt"
     isfile(name) && rm(name)
     GaPSE.print_map_ξ_LDxGNC_multipole(COSMO, name, effect, 10 .^ range(-1, 3, length=KWARGS_LDxGNC[:N_log]);
          L = L, GaPSE.specif_kwargs_LDxGNC(effect, KWARGS_LDxGNC)...)

     calc_table = readdlm(name; comments = true)
     calc_ss = convert(Vector{Float64}, calc_table[:, 1])
     calc_xis = convert(Vector{Float64}, calc_table[:, 2])

     @test all([isapprox(s, calc_s, rtol = 1e-2) for (s, calc_s) in zip(ss, calc_ss)])
     @test all([isapprox(xi, calc_xi, rtol = 1e-2) for (xi, calc_xi) in zip(xis, calc_xis)])

     rm(name)
end

