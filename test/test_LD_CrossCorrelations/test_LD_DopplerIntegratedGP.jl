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

@testset "test xi doppler_integratedgp L = 0" begin
    effect = "doppler_integratedgp"
    L = 0

    ss = ss_LD_L0_noF
    xis = all_xis_LD_L0_noF[GaPSE.INDEX_GR_EFFECT_LD[effect]]

    name = "calc_xi_LD_" * effect * "_LD_L$L" * ".txt"
    isfile(name) && rm(name)
    GaPSE.print_map_ξ_LD_multipole(COSMO, name, effect, SS_LD;
        L = L, alg = :quad, GaPSE.specif_kwargs_LD(effect, KWARGS_LD)...)

    calc_table = readdlm(name; comments = true)
    calc_ss = convert(Vector{Float64}, calc_table[:, 1])
    calc_xis = convert(Vector{Float64}, calc_table[:, 2])

    @test all([isapprox(s, calc_s, rtol = 1e-2) for (s, calc_s) in zip(ss, calc_ss)])
    @test all([isapprox(xi, calc_xi, rtol = 1e-2) for (xi, calc_xi) in zip(xis, calc_xis)])

    rm(name)
end

@testset "test xi integratedgp_doppler L = 0" begin
    effect = "integratedgp_doppler"
    L = 0

    ss = ss_LD_L0_noF
    xis = all_xis_LD_L0_noF[GaPSE.INDEX_GR_EFFECT_LD[effect]]

    name = "calc_xi_LD_" * effect * "_LD_L$L" * ".txt"
    isfile(name) && rm(name)
    GaPSE.print_map_ξ_LD_multipole(COSMO, name, effect, SS_LD;
        L = L, alg = :quad, GaPSE.specif_kwargs_LD(effect, KWARGS_LD)...)

    calc_table = readdlm(name; comments = true)
    calc_ss = convert(Vector{Float64}, calc_table[:, 1])
    calc_xis = convert(Vector{Float64}, calc_table[:, 2])

    @test all([isapprox(s, calc_s, rtol = 1e-2) for (s, calc_s) in zip(ss, calc_ss)])
    @test all([isapprox(xi, calc_xi, rtol = 1e-2) for (xi, calc_xi) in zip(xis, calc_xis)])

    rm(name)
end