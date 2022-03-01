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

using Pkg
Pkg.activate(normpath(@__DIR__))
using GaPSE
using DelimitedFiles

FILE_NAME = split(PROGRAM_FILE, "/")[end]

#main(x::Union{String, Float64, Int64}...) = main([string(var) for var in [x...]])
function main()
     FILE_F_MAP = "data/F_REFERENCE.txt"
     #FILE_PS = "data/WideA_ZA_pk.dat"
     FILE_PS = "test/datatest/file_pk.txt"
     FILE_BACKGROUND = "data/WideA_ZA_background.dat"

     #casto_tab = readdlm("data/tab_kappa_terms.dat", comments = true)
     #casto_ss = convert(Vector{Float64}, casto_tab[5:end, 1])


     z_min = 0.05
     z_max = 0.20
     θ_max = π / 2.0
     params = GaPSE.CosmoParams(z_min, z_max, θ_max;
          k_min = 1e-8, k_max = 10.0,
          Ω_b = 0.0489, Ω_cdm = 0.251020, h_0 = 0.70)
     cosmo = GaPSE.Cosmology(params,
          FILE_BACKGROUND,
          FILE_PS,
          FILE_F_MAP)
     GaPSE.parameters_used(stdout, cosmo)

     #=
     GaPSE.print_map_map_ξ_multipole(cosmo, "outputs/xi_doppler_L0.txt", "auto_doppler";
          μ_atol = 1e-4, μ_rtol = 1e-3, use_windows = true)

     GaPSE.print_PS_multipole("outputs/xi_doppler.txt", "outputs/P_doppler.txt", nothing)
     =#

     #GaPSE.print_map_ξ_multipole(cosmo, 
     #     "outputs/xi_lensing_L0.txt", "auto_lensing"; use_windows = false, 
     #     N_χs = 50, Δχ_min = 1e-4, enhancer=1e10, μ_rtol=1e-2)

     xs = [x for x in 0:0.25:3]
     μs = vcat([-1.0, -0.98, -0.95], [μ for μ in -0.9:0.1:0.9], [0.95, 0.98, 1.0])
     #μs = vcat([μ for μ in -1:0.01:-0.91], [μ for μ in -0.9:0.1:0.9], [μ for μ in 0.91:0.01:1.0])
     GaPSE.F_map(xs, μs; out = "outputs/F.txt", rtol=5e-3, atol=1e-2)

end

if (ARGS == String[])
     #println("\nwithout input arguments/commands, show this help message and exit\n")
     #main(["--help"])
     main()
else
     #main(ARGS)
     return 0
end
