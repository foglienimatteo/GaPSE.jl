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

PATH_TO_GAPSE ="../"
#PATH_TO_GAPSE = "/dss/dsshome1/08/di75tom/julia_with_salvo/GaPSE.jl/"
Pkg.activate(normpath(PATH_TO_GAPSE)) # or @__DIR__
using GaPSE


const global N_LOB = 1000
const global N_TRAP = 1000
const global ALG = :lobatto
const global N_CHIS = 100
const global N_CHIS_2 = 50
const global SUIT_SAMPLING = true
const global LENGTH = 500
const global OM = true
COMPUTE_XIS_GNC = true # set this to true if you want to compute the TPCFs of the GNC!
# This is the directory name where to put the files computed in the 
# following steps; makes sure that it's name ends with "/" !
#DIR = "/dss/dsshome1/08/di75tom/julia_with_salvo/"
DIR="./"
@assert isdir(DIR) "ERROR: DIR=$DIR DOESN'T EXIST!!!"

using DelimitedFiles, Dierckx, Printf

FILE_NAME = split(PROGRAM_FILE, "/")[end]

#main(x::Union{String, Float64, Int64}...) = main([string(var) for var in [x...]])
function main()

     # run in the command line:
     #
     #    $ julia gapse-computations.jl


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
          :n_extrap_low => 300, :n_extrap_high => 300,
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

     params = GaPSE.CosmoParams(z_min, z_max, θ_max;
          Ω_b=0.0489, Ω_cdm=0.251020, 
          h_0 = 0.70, s_lim = 1e-2, z_spline_lim = 1000.0,
          b1=1.5, s_b1=0.0, 𝑓_evo1=0.0, 
          b2 = nothing, s_b2 = nothing, 𝑓_evo2 = nothing,
          IPS_opts=Dict(
               :fit_left_min => 1e-6, :fit_left_max => 3e-6,
               :fit_right_min => 1e1, :fit_right_max => 2e1),
          IPSTools_opts=Dict(
               :N => 1024, :fit_min => 0.05, :fit_max => 0.5,
               :con => true, :k_min => 1e-8, :k_max => 10.0)
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

     name_xis_GNC_L0_noF_noobsvel_file = DIR * "xis_GNC_L0_noF_noobsvel_GenWin.txt"
     name_xis_GNC_L1_noF_noobsvel_file = DIR * "xis_GNC_L1_noF_noobsvel_GenWin.txt"
     name_xis_GNC_L2_noF_noobsvel_file = DIR * "xis_GNC_L2_noF_noobsvel_GenWin.txt"
     name_xis_GNC_L3_noF_noobsvel_file = DIR * "xis_GNC_L3_noF_noobsvel_GenWin.txt"
     name_xis_GNC_L4_noF_noobsvel_file = DIR * "xis_GNC_L4_noF_noobsvel_GenWin.txt"


     if (COMPUTE_XIS_GNC == true)
          GaPSE.print_map_sum_ξ_GNC_multipole(
               cosmo, name_xis_GNC_L0_noF_noobsvel_file,
               10 .^ range(0, log10(2 * cosmo.s_max), length=LENGTH);
               use_windows=false, L=0, alg=ALG, obs=:noobsvel,
               single=true, enhancer=1e8,
               N_trap=N_TRAP, N_lob=N_LOB, atol_quad=0.0, rtol_quad=1e-2,
               N_χs=N_CHIS, N_χs_2=N_CHIS_2, suit_sampling=SUIT_SAMPLING
          )
     end

     if (COMPUTE_XIS_GNC == true && OM == false)
          GaPSE.print_map_sum_ξ_GNC_multipole(
               cosmo, name_xis_GNC_L1_noF_noobsvel_file,
               10 .^ range(0, log10(2 * cosmo.s_max), length=LENGTH);
               use_windows=false, L=1, alg=ALG, obs=:noobsvel,
               single=true, enhancer=1e8,
               N_trap=N_TRAP, N_lob=N_LOB, atol_quad=0.0, rtol_quad=1e-2,
               N_χs=N_CHIS, N_χs_2=N_CHIS_2, suit_sampling=SUIT_SAMPLING
          )
     end

     if (COMPUTE_XIS_GNC == true && OM == false)
          GaPSE.print_map_sum_ξ_GNC_multipole(
               cosmo, name_xis_GNC_L2_noF_noobsvel_file,
               10 .^ range(0, log10(2 * cosmo.s_max), length=LENGTH);
               use_windows=false, L=2, alg=ALG, obs=:noobsvel,
               single=true, enhancer=1e8,
               N_trap=N_TRAP, N_lob=N_LOB, atol_quad=0.0, rtol_quad=1e-2,
               N_χs=N_CHIS, N_χs_2=N_CHIS_2, suit_sampling=SUIT_SAMPLING
          )
     end

     if (COMPUTE_XIS_GNC == true && OM == false)
          GaPSE.print_map_sum_ξ_GNC_multipole(
               cosmo, name_xis_GNC_L3_noF_noobsvel_file,
               10 .^ range(0, log10(2 * cosmo.s_max), length=LENGTH);
               use_windows=false, L=3, alg=ALG, obs=:noobsvel,
               single=true, enhancer=1e8,
               N_trap=N_TRAP, N_lob=N_LOB, atol_quad=0.0, rtol_quad=1e-2,
               N_χs=N_CHIS, N_χs_2=N_CHIS_2, suit_sampling=SUIT_SAMPLING
          )
     end

     if (COMPUTE_XIS_GNC == true && OM == false)
          GaPSE.print_map_sum_ξ_GNC_multipole(
               cosmo, name_xis_GNC_L4_noF_noobsvel_file,
               10 .^ range(0, log10(2 * cosmo.s_max), length=LENGTH);
               use_windows=false, L=4, alg=ALG, obs=:noobsvel,
               single=true, enhancer=1e8,
               N_trap=N_TRAP, N_lob=N_LOB, atol_quad=0.0, rtol_quad=1e-2,
               N_χs=N_CHIS, N_χs_2=N_CHIS_2, suit_sampling=SUIT_SAMPLING
          )
     end


     ######## create_file_for_XiMultipoles function #########



     if OM == false
          name_file_ximultipoles = DIR * "xis_GNC_L01234_noF_noobsvel.txt"

          vector_xis_names = [
               name_xis_GNC_L0_noF_noobsvel_file,
               name_xis_GNC_L1_noF_noobsvel_file,
               name_xis_GNC_L2_noF_noobsvel_file,
               name_xis_GNC_L3_noF_noobsvel_file,
               name_xis_GNC_L4_noF_noobsvel_file,
          ]

          GaPSE.create_file_for_XiMultipoles(name_file_ximultipoles, vector_xis_names, 2, "GNC")

          #run(`head -n 18 $(name_file_ximultipoles)`)

          ps_kwargs_2(alg::Symbol=:fftlog) = alg == :twofast ?
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

          ks, pks = GaPSE.PS_multipole_GenWin(name_file_ximultipoles, name_file_Qmultipoles; L=0, ps_kwargs_2(:fftlog)...)

          out = DIR * "results.txt"
          isfile(out) && run(`rm $out`)
          open(out, "w") do io
               println(io, BRAND)

               println(io, "# Power Spectrum GenWin monopole computation.")
               println(io, "#\n# For this PS_multipole computation we set: ")
               println(io, "# \t algorithm chosen for the computation: :$ALG")
               println(io, "# \t N_LOB = $N_LOB \t N_TRAP = $N_TRAP \t LENGTH = $LENGTH")
               println(io, "# \t N_LOB = $N_LOB \t N_TRAP = $N_TRAP \t LENGTH = $LENGTH")
               #println(io, "# computational time needed (in s) : $(@sprintf("%.4f", time_2-time_1))")
               print(io, "# kwards passed to \"print_PS_multipole\": ")

               if isempty(ps_kwargs_2)
                    println(io, "none")
               else
                    print(io, "\n")
                    for key in keys(ps_kwargs_2)
                         println(io, "# \t\t$(key) = $(ps_kwargs_2[key])")
                    end
               end
               println(io, "# ")
               println(io, "# k [h_0/Mpc] \t \t  P [(Mpc/h_0)^3]")
               for (k, pk) in zip(ks, pks)
                    println(io, "$k \t " * GaPSE.number_to_string(pk))
               end
          end
     end

     #=
     GaPSE.print_all_PS_multipole(DIR * name_xis_GNC_L0_noF_noobsvel_file,
          DIR * name_ps_GNC_L0_noF_noobsvel_file, "GNC"; L=0, ps_kwargs(tf)...)
     ks_GNC_L0_noF_noobsvel, pks_sum_GNC_L0_noF_noobsvel, pks_all_GNC_L0_noF_noobsvel =
          GaPSE.readxyall(DIR * name_ps_GNC_L0_noF_noobsvel_file)
     spline_GNCsum_pks_L0_noF_noobsvel = Spline1D(ks_GNC_L0_noF_noobsvel, pks_sum_GNC_L0_noF_noobsvel; bc="error")

     GaPSE.print_all_PS_multipole(DIR * name_xis_GNC_L0_noF_noobs_file,
          DIR * name_ps_GNC_L0_noF_noobs_file, "GNC"; L=0, ps_kwargs(tf)...)
     ks_GNC_L0_noF_noobs, pks_sum_GNC_L0_noF_noobs, pks_all_GNC_L0_noF_noobs =
          GaPSE.readxyall(DIR * name_ps_GNC_L0_noF_noobs_file)
     spline_GNCsum_pks_L0_noF_noobs = Spline1D(ks_GNC_L0_noF_noobs, pks_sum_GNC_L0_noF_noobs; bc="error")

     GaPSE.print_all_PS_multipole(DIR * name_xis_GNC_L0_withF_noobsvel_file,
          DIR * name_ps_GNC_L0_withF_noobsvel_file, "GNC"; L=0, ps_kwargs(tf)...)
     ks_GNC_L0_withF_noobsvel, pks_sum_GNC_L0_withF_noobsvel, pks_all_GNC_L0_withF_noobsvel =
          GaPSE.readxyall(DIR * name_ps_GNC_L0_withF_noobsvel_file)
     spline_GNCsum_pks_L0_withF_noobsvel = Spline1D(ks_GNC_L0_withF_noobsvel, pks_sum_GNC_L0_withF_noobsvel; bc="error")

     GaPSE.print_all_PS_multipole(DIR * name_xis_GNC_L0_withF_noobs_file,
          DIR * name_ps_GNC_L0_withF_noobs_file, "GNC"; L=0, ps_kwargs(tf)...)
     ks_GNC_L0_withF_noobs, pks_sum_GNC_L0_withF_noobs, pks_all_GNC_L0_withF_noobs =
          GaPSE.readxyall(DIR * name_ps_GNC_L0_withF_noobs_file)
     spline_GNCsum_pks_L0_withF_noobs = Spline1D(ks_GNC_L0_withF_noobs, pks_sum_GNC_L0_withF_noobs; bc="error")


     ks_Newtonian_L0_noF_noobsvel, pks_Newtonian_L0_noF_noobsvel =
          ks_GNC_L0_noF_noobsvel, pks_all_GNC_L0_noF_noobsvel[GaPSE.INDEX_GR_EFFECT_GNC["auto_newton"]]
     spline_Newt_pks_L0_noF_noobsvel = Spline1D(ks_Newtonian_L0_noF_noobsvel, pks_Newtonian_L0_noF_noobsvel; bc="error")
     ks_Newtonian_L0_noF_noobs, pks_Newtonian_L0_noF_noobs =
          ks_GNC_L0_noF_noobs, pks_all_GNC_L0_noF_noobs[GaPSE.INDEX_GR_EFFECT_GNC["auto_newton"]]
     spline_Newt_pks_L0_noF_noobs = Spline1D(ks_Newtonian_L0_noF_noobs, pks_Newtonian_L0_noF_noobs; bc="error")

     ks_Newtonian_L0_withF_noobsvel, pks_Newtonian_L0_withF_noobsvel =
          ks_GNC_L0_withF_noobsvel, pks_all_GNC_L0_withF_noobsvel[GaPSE.INDEX_GR_EFFECT_GNC["auto_newton"]]
     spline_Newt_pks_L0_withF_noobsvel = Spline1D(ks_Newtonian_L0_withF_noobsvel, pks_Newtonian_L0_withF_noobsvel; bc="error")
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
