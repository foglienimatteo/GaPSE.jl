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


     ########## This is the basic configuration for a Cosmology ###########

     
     # This is the directory name where to put the files computed in the 
     # following steps; makes sure that it's name ends with "/" !
     DIR = "dir_name/";
     isdir(DIR) || mkdir(DIR)

     PATH_TO_GAPSE = "./"
     FILE_F_MAP = PATH_TO_GAPSE * "data/F_REFERENCE_pi2.txt";

     #=
     xs = [x for x in 0:0.02:5]
     μs = union(
          [μ for μ in range(-1.0, -0.95, length = 100)], 
          [μ for μ in range(-0.95, 0.95, length = 100)],
          [μ for μ in range(0.95, 1.0, length = 100)]
     )
     GaPSE.print_map_F(FILE_F_MAP, 
     xs, μs; 
     trap = true,
     Fmap_opts = Dict(
               :θ_max => π / 2.0, :tolerance => 1e-8, 
               :N => 1000, :pr => true
          )
     )
     =#

     FILE_BACKGROUND = PATH_TO_GAPSE * "data/WideA_ZA_background.dat"

     z_min, z_max, θ_max = 1.0, 1.5, π / 2.0
     WFI_opts = Dict(
          :ss_start => 0.0, :ss_stop => 0.0,
          :ss_step => 100, :llim => 0.0, :rlim => Inf,
          :rtol => 5e-2, :atol => 0.0, :N => 1000, #:pr => true,
     )

     params = GaPSE.CosmoParams(z_min, z_max, θ_max;
          Ω_b=0.0489, Ω_cdm=0.251020, h_0=0.70, s_lim=1e-2,
          s_b=0.0, 𝑓_evo=0.0, b=1.5,
          IPS_opts=Dict(
               :fit_left_min => 1e-6, :fit_left_max => 3e-6,
               :fit_right_min => 1e1, :fit_right_max => 2e1),
          IPSTools_opts=Dict(
               :N => 1024, :fit_min => 0.05, :fit_max => 0.5,
               :con => true, :k_min => 1e-8, :k_max => 10.0),
          WFI_opts=WFI_opts
     )

     FILE_IF_MAP = PATH_TO_GAPSE * "data/IntegrF_REFERENCE_pi2_z005020.txt"

     #=
     new_calc_μs = union([μ for μ in -1.0:0.0001:(-0.98)], 
          [μ for μ in -0.98:0.01:0.98], 
          [μ for μ in 0.98:0.0001:1.0]);

     GaPSE.print_map_IntegratedF(
          FILE_F_MAP, FILE_IF_MAP, 
          z_min, z_max, new_calc_μs, FILE_BACKGROUND;
          trap = true, WFI_opts...
          )
     =#

     FILE_PS = PATH_TO_GAPSE * "test/datatest/file_pk.txt"
     cosmo = GaPSE.Cosmology(params, FILE_BACKGROUND, FILE_PS, FILE_F_MAP, FILE_IF_MAP)

     GaPSE.parameters_used(stdout, cosmo)



     ###### Here is the code you want to run ##########

     GaPSE.print_map_sum_ξ_GNC_multipole(
          cosmo, DIR * "GNC_sb0_fevo0_L0_noF_noobsvel_quad.txt",
          10 .^ range(0, log10(2 * cosmo.s_max), length=100);
          L=0, use_windows=false, alg=:quad, obs=:noobsvel,
          N_lob=1000, N_trap=5, atol_quad=0.0, rtol_quad=1e-2,
          N_χs=100, N_χs_2=60
     );

     ps_kwargs(twofast::Bool) = twofast ?
          Dict(
               :twofast => true, :epl => true, :pr => false,
               :N_left => 12, :N_right => 12,
               :p0_left => [-2.0, 1.0], :p0_right => [-2.0, 1.0],
               :int_s_min => 1e0, :int_s_max => 1200.0
          ) : Dict(
               :twofast => false, :pr => true, :ν => 1.5, 
               :n_extrap_low => 0, :n_extrap_high => 0, 
               :n_pad => 500, :cut_first_n => 6
          )

     GaPSE.print_all_PS_multipole(DIR * "GNC_sb0_fevo0_L0_noF_noobsvel_quad.txt",
          DIR * "PS_sb0_fevo0_L0_noF_noobsvel_quad.txt", "GNC"; L=0, ps_kwargs(false)...)


end

if (ARGS == String[])
     #println("\nwithout input arguments/commands, show this help message and exit\n")
     #main(["--help"])
     main()
else
     #main(ARGS)
     return 0
end
