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
     # Have you followed the "ipynb/TUTORIAL.ipynb" file? 
     # All this stuff is there explained!
     # If you already know it, just wrote the code you want execute and
     # then run in the command line:
     #
     #    $ julia GaPSE-exe.jl
     #
     # Have a nice day!

     FILE_F_MAP = "data/F_REFERENCE.txt"
     FILE_PS = "data/WideA_ZA_pk.dat"
     FILE_BACKGROUND = "data/WideA_ZA_background.dat"

     z_min, z_max, θ_max = 0.05, 0.20, π / 2.0

     params = GaPSE.CosmoParams(z_min, z_max, θ_max;
          Ω_b = 0.0489, Ω_cdm = 0.251020, h_0 = 0.70, s_lim = 1e-2,
          IPS_opts = Dict(
               :fit_left_min => 1e-6, :fit_left_max => 3e-6,
               :fit_right_min => 1e1, :fit_right_max => 2e1),
          IPSTools_opts = Dict(
               :N => 1024, :fit_min => 0.05, :fit_max => 0.5,
               :con => true, :k_min => 1e-8, :k_max => 10.0)
     )

     cosmo = GaPSE.Cosmology(params, FILE_BACKGROUND, FILE_PS, FILE_F_MAP)

     GaPSE.parameters_used(stdout, cosmo)

     
     GaPSE.print_map_ξ_multipole(cosmo, "my_first_doppler.txt", 
          "auto_doppler", 10 .^ range(0, 3, length = 1000); 
          L = 0, use_windows = false);

     GaPSE.print_PS_multipole("my_first_doppler.txt",
          "my_first_ps_doppler.txt";
          L = 0, N = 1024, int_s_min = 1e-3, int_s_max = 1e3)

     #=
     GaPSE.print_map_sum_ξ_multipole(cosmo, "my_first_all_xis.txt",
          10 .^ range(0, 3, length = 100), use_windows = false);
     =#

     #=
     xs = [x for x in 0:0.25:3]
     μs = vcat([-1.0, -0.98, -0.95], [μ for μ in -0.9:0.1:0.9], [0.95, 0.98, 1.0])
     #μs = vcat([μ for μ in -1:0.01:-0.91], [μ for μ in -0.9:0.1:0.9], [μ for μ in 0.91:0.01:1.0])
     GaPSE.F_map(xs, μs; out = "outputs/F.txt", rtol = 5e-3, atol = 1e-2)
     =#

end

if (ARGS == String[])
     #println("\nwithout input arguments/commands, show this help message and exit\n")
     #main(["--help"])
     main()
else
     #main(ARGS)
     return 0
end
