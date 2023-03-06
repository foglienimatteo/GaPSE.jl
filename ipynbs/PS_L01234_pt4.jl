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

PATH_TO_GAPSE = "../"
Pkg.activate(normpath(PATH_TO_GAPSE)) # or @__DIR__
using GaPSE


const global B = 1.5     # Galaxy bias
const global S_B = 0.0   # Magnification bias
const global F_EVO = 0.0 # Evolution bias
const global filenames_appendix = "_b1p5-sb0-fevo0"
COMPUTE_XIS_GNC = true # set this to true if you want to compute the TPCFs of the GNC!
# This is the directory name where to put the files computed in the 
# following steps; makes sure that it's name ends with "/" !
DIR = "PS_L01234/"
@assert isdir(DIR) "ERROR: DIR=$DIR DOESN'T EXIST!!!"

using DelimitedFiles, Dierckx, Printf

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

     ps_kwargs(alg::Symbol=:fftlog) = alg == :twofast ?
                                      Dict(
          :alg => :twofast, :epl => true, :pr => false,
          :N_left => 12, :N_right => 12,
          :p0_left => [-2.0, 1.0], :p0_right => [-2.0, 1.0],
          :int_s_min => 1e0, :int_s_max => 1200.0,
          :cut_first_n => 0, :cut_last_n => 0
     ) : alg == :fftlog ?
                                      Dict(
          :alg => :fftlog, :pr => true, :ν => 1.5,
          :n_extrap_low => 0, :n_extrap_high => 0,
          :n_pad => 500, :cut_first_n => 0, :cut_last_n => 0,
     ) : throw(AssertionError("alg = :fftlog (recommended) or alg = :twofast !"))
     tf = :fftlog

     FILE_F_MAP = PATH_TO_GAPSE * "data/NEW_F_pi2.txt"
     #FILE_F_MAP =  PATH_TO_GAPSE * "data/F_REFERENCE_pi2.txt";
     #=
     kwargs_map_F_hcub = Dict(
          :θ_max => π / 2.0, :tolerance => 1e-10, 
          :rtol => 1e-2, :atol => 1e-3, :pr => true,
     );

     kwargs_map_F_trap = Dict(
          :θ_max => π / 2.0, :tolerance => 1e-10, 
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

     z_min, z_max, θ_max = 1.0, 1.5, π / 2.0

     FILE_BACKGROUND = PATH_TO_GAPSE * "data/WideA_ZA_background.dat"

     #=
     WFI_opts = Dict(
         :llim => nothing, :rlim => nothing, 
         :rtol => 1e-2, :atol => 0.0, 
         :N => 1000, :pr => true,
     )
     =#

     params = GaPSE.CosmoParams(z_min, z_max, θ_max;
          Ω_b=0.0489, Ω_cdm=0.251020, h_0=0.70, s_lim=1e-2,
          s_b=S_B, 𝑓_evo=F_EVO, b=B,
          IPS_opts=Dict(
               :fit_left_min => 1e-6, :fit_left_max => 3e-6,
               :fit_right_min => 1e1, :fit_right_max => 2e1),
          IPSTools_opts=Dict(
               :N => 1024, :fit_min => 0.05, :fit_max => 0.5,
               :con => true, :k_min => 1e-8, :k_max => 10.0)
          #=
          WFI_opts = Dict(
              :llim => nothing, :rlim => nothing, 
              :rtol => 1e-2, :atol => 0.0, 
              :N => 1000, :pr => true,)
          =#
     )

     #FILE_F_MAP =  PATH_TO_GAPSE * "data/NEW_F_pi2.txt";
     FILE_IF_MAP = PATH_TO_GAPSE * "data/NEW_IntegrF_pi2_z115.txt"

     #FILE_F_MAP = PATH_TO_GAPSE*"data/F_REFERENCE_pi2.txt";
     #FILE_IF_MAP = PATH_TO_GAPSE*"data/IntegrF_REFERENCE_pi2_z115.txt";

     #=
     calc_μs = union(
     [μ for μ in range(-1.0, -0.98, length = 50)], 
     [μ for μ in range(-0.98, 0.98, length = 102)],
     [μ for μ in range(0.98, 1.0, length = 50)]);

     GaPSE.print_map_IntegratedF(
     z_min, z_max, calc_μs,
     FILE_F_MAP, FILE_IF_MAP, 
     FILE_BACKGROUND;
     alg = :trap, N_ss = 200, m = 2.1,
     Dict(
          :llim => nothing, :rlim => nothing, 
          :rtol => 1e-2, :atol => 0.0, 
          :N => 1000, :pr => true,
     )...
     )
     =#


     FILE_PS = PATH_TO_GAPSE * "test/datatest/file_pk.txt"
     cosmo = GaPSE.Cosmology(params, FILE_BACKGROUND, FILE_PS, FILE_F_MAP, FILE_IF_MAP)

     GaPSE.parameters_used(stdout, cosmo)



     ######## Computation of the GNC TPCFs and Newtonian TPCF (and their PS) #########

     VEC_L = [4]
     vec_name_xis_GNC_noF_noobs_file = [
          "xis_GNC_L$L" * "_withF_noobs" * filenames_appendix * ".txt" for L in VEC_L]
     vec_name_ps_GNC_noF_noobs_file = [
          "ps_GNC_L$L" * "_withF_noobs" * filenames_appendix * ".txt" for L in VEC_L]


     if COMPUTE_XIS_GNC == true
          for (i, L) in enumerate(VEC_L)
               GaPSE.print_map_sum_ξ_GNC_multipole(
                    cosmo, DIR * vec_name_xis_GNC_noF_noobs_file[i],
                    10 .^ range(0, log10(2 * cosmo.s_max), length=500);
                    use_windows=true, L=L, alg=:quad, obs=:no,
                    single=true, enhancer=1e8,
                    N_trap=200, N_lob=500, atol_quad=0.0, rtol_quad=1e-2,
                    N_χs=100, N_χs_2=60
               )
          end
     end


     for (i, L) in enumerate(VEC_L)
          GaPSE.print_all_PS_multipole(DIR * vec_name_xis_GNC_noF_noobs_file[i],
               DIR * vec_name_ps_GNC_noF_noobs_file[i], "GNC"; L=L, ps_kwargs(tf)...)
     end

end







if (ARGS == String[])
     #println("\nwithout input arguments/commands, show this help message and exit\n")
     #main(["--help"])
     main()
else
     #main(ARGS)
     return 0
end
