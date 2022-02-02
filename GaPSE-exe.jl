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
     FILE_PS = "data/WideA_ZA_pk.dat"
     FILE_BACKGROUND = "data/WideA_ZA_background.dat"

     z_min = 0.05
     z_max = 0.20
     BD = GaPSE.BackgroundData(FILE_BACKGROUND, z_min, z_max; h=0.7)
     BS = GaPSE.BackgroundSplines(BD)

     #xs = [x for x in 0:0.1:3]
     #μs = vcat([μ for μ in -1:0.01:-0.91], [μ for μ in -0.9:0.1:0.9], [μ for μ in 0.91:0.01:1.0])
     #GaPSE.F_map(xs, μs; out = "data/F_REFERENCE.txt", rtol=5e-3, atol=1e-2)

     #GaPSE.print_map_int_on_mu_lensing("outputs/xi_lensing.txt"; 
     #     χ_atol = 5e-3, χ_rtol=1e-3, atol = 1e-3, rtol = 1e-3, tol=1)
     #GaPSE.print_PS_multipole("outputs/P_lensing.txt", "outputs/xi_lensing.txt")

     casto_tab = readdlm("data/tab_kappa_terms.dat", comments = true)
     casto_ss = convert(Vector{Float64}, casto_tab[5:end, 1])
     #GaPSE.print_map_int_on_mu_lensing("outputs/xi_lensing.txt", casto_ss;
     #     χ_atol = 1e-4, χ_rtol = 5e-3, atol = 1e-4, rtol = 1e-3, tol = 1.0, Δχ_min = 1.0,
     #     use_windows = false)
     #
     GaPSE.print_map_int_on_mu("outputs/xi_lensing.txt", "auto_lensing", casto_ss;
          χ_atol = 1e-4, χ_rtol = 1e-2, μ_atol = 1e-4, μ_rtol = 1e-2, Δχ_min = 1.0,
          use_windows = false
     )
     #GaPSE.print_PS_multipole("outputs/P_lensing.txt", "outputs/xi_lensing.txt")
     #GaPSE.print_map_int_on_mu("outputs/xi_doppler.txt", "auto_doppler", casto_ss;
     #     μ_atol = 1e-4, μ_rtol = 1e-3, use_windows = true)
     #GaPSE.print_PS_multipole("outputs/P_doppler.txt", "outputs/xi_doppler.txt")



     #GaPSE.integral_on_mu_lensing(GaPSE.s_eff, 1.0)

     #s1 = 435.3747096079416
     #s2 = 436.3725783407533
     #s2 = 480.0
     #y = 0.9999999887939074
     #y=0.86
     #t1 = time()
     #GaPSE.ξ_lensing(s1, s2, y; enhancer = 1e10, rtol=1e-2, atol=1e-4, Δχ_min=1)
     #t2 = time()
     #println("time = $(t2-t1)")

     #GaPSE.print_map_int_on_mu_doppler("xi_doppler.txt",
     #     atol = 1e-3, rtol = 1e-2, tol = 0.1, enhancer=1)
     #GaPSE.print_PS("P_doppler.txt", "xi_doppler.txt")
     #GaPSE.print_PS("P_doppler.txt", "auto_doppler")

end

if (ARGS == String[])
     #println("\nwithout input arguments/commands, show this help message and exit\n")
     #main(["--help"])
     main()
else
     #main(ARGS)
     return 0
end
