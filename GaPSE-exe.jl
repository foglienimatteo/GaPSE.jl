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
     z_min, z_max, θ_max = 1.0, 1.5, π / 2.0;

     #=
     kwargs_map_F_hcub = Dict(
          :θ_max => θ_max, :tolerance => 1e-10, 
          :rtol => 1e-2, :atol => 1e-3, :pr => true,
     );

     kwargs_map_F_trap = Dict(
          :θ_max => θ_max, :tolerance => 1e-10, 
          :N => 1000, :pr => true,
     );

     xs = [x for x in 0:0.02:5]
     μs = union(
          [μ for μ in range(-1.0, -0.98, length = 50)], 
          [μ for μ in range(-0.98, 0.98, length = 102)],
          [μ for μ in range(0.98, 1.0, length = 50)]);
     GaPSE.print_map_F(FILE_F_MAP, xs, μs; 
          alg = :trap, Fmap_opts = kwargs_map_F_trap # we recommend to use :trap
          #alg = :hcub, Fmap_opts = kwargs_map_F_hcub # but you can use also :hcub if you prefer
     )
     =#

     FILE_BACKGROUND = PATH_TO_GAPSE * "data/WideA_ZA_background.dat"

     # NOTE: the bias-like parameters are now specified per-species: `b1`, `s_b1`
     # and `𝑓_evo1` refer to the first species, while `b2`, `s_b2` and `𝑓_evo2`
     # refer to the second one. If you leave the latter to `nothing` (the
     # default), they are set equal to the former ones, i.e. you are doing an
     # auto-correlation of a single species.
     #WFI_opts = Dict(
     #     :ss_start => 0.0, :ss_stop => 0.0,
     #     :ss_step => 100, :llim => 0.0, :rlim => Inf,
     #     :rtol => 5e-2, :atol => 0.0, :N => 1000, #:pr => true,
     #)
     params = GaPSE.CosmoParams(z_min, z_max, θ_max;
          Ω_b=0.0489, Ω_cdm=0.251020, h_0=0.70, s_lim=1e-2,
          b1=1.5, s_b1=0.0, 𝑓_evo1=0.0,
          b2=nothing, s_b2=nothing, 𝑓_evo2=nothing,
          IPS_opts=Dict(
               :fit_left_min => 1e-6, :fit_left_max => 3e-6,
               :fit_right_min => 1e1, :fit_right_max => 2e1),
          IPSTools_opts=Dict(
               :N => 1024, :fit_min => 0.05, :fit_max => 0.5,
               :con => true, :k_min => 1e-8, :k_max => 10.0),
          #WFI_opts=WFI_opts
     )

     # This integrated window function map must be computed for the same
     # (z_min, z_max, θ_max) of `params`; the one below matches z_min=1.0,
     # z_max=1.5 and θ_max=π/2.
     FILE_IF_MAP = PATH_TO_GAPSE * "data/IntegrF_REFERENCE_pi2_z115.txt"

     #=
     double_z_max = GaPSE.corresponding_redshift(z_max, 3.0, FILE_BACKGROUND)
     calc_zs = [z for z in range(0.0, double_z_max, length=100)]

     calc_μs = union(
          [μ for μ in range(-1.0, -0.98, length = 50)], 
          [μ for μ in range(-0.98, 0.98, length = 102)],
          [μ for μ in range(0.98, 1.0, length = 50)]);

     GaPSE.print_map_IntegratedF(
          z_min, z_max, calc_zs, calc_μs,
          PATH_TO_GAPSE * "data/F_REFERENCE_pi2.txt", 
          PATH_TO_GAPSE * "data/IntegrF_REFERENCE_pi2_z115.txt", 
          FILE_BACKGROUND;
          alg = :trap, Dict(
               :llim => nothing, :rlim => nothing, 
               :rtol => 1e-2, :atol => 0.0, 
               :N => 1000, :pr => true,
          )...
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

     function ps_kwargs(alg::Symbol = :fftlog) 
          if alg == :twofast
               return Dict(
                    :alg => :twofast, :L => 0, :pr => true,
                    :cut_first_n => 6, :cut_last_n => 4,
                    :epl => true, 
                    :N_left => 12, :N_right => 12,
                    :p0_left => [-2.0, 1.0], :p0_right => [-2.0, 1.0],
                    :int_s_min => 1e0, :int_s_max => 1200.0,
                    #cut_first_n => 6, cut_last_n => 3,
               ) 
          elseif alg == :fftlog
          
               return Dict(
                    :alg => :fftlog, :L => 0,
                    :cut_first_n => 6, :cut_last_n => 4,
                    :ν => 1.5, :n_pad => 500,
                    :n_extrap_low => 0, :n_extrap_high => 0,  
               )
          else
               throw(AssertionError("alg = :$alg is not a valid algorithm for PS_multipole"))
          end
     end

     GaPSE.print_all_PS_multipole(DIR * "GNC_sb0_fevo0_L0_noF_noobsvel_quad.txt",
          DIR * "PS_sb0_fevo0_L0_noF_noobsvel_quad.txt", "GNC"; L=0, ps_kwargs(:fftlog)...)


end

if (ARGS == String[])
     #println("\nwithout input arguments/commands, show this help message and exit\n")
     #main(["--help"])
     main()
else
     #main(ARGS)
     return 0
end
