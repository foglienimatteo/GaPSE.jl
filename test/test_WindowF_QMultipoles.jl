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

@testset "test WindowFIntegrated_multipole" begin
     windowfint = GaPSE.WindowFIntegrated("datatest/IntegrF_REFERENCE_pi2_z115.txt")
     results = [4.1857800000750543e11, 3.889433126687778e11, 1.3718401137777765e11, 4.302441469372041e10]

     for (s, res) in zip([1.0, 100.0, 1000.0, 2500.0], results)
          val = GaPSE.WindowFIntegrated_multipole(s, windowfint;
               s_min=2312.5774317555542, s_max=3054.483020067565, 
               L=0, alg=:lobatto, N_lob=100, N_trap=200, atol_quad=0.0, rtol_quad=1e-2, enhancer=1e6)
          @test isapprox(val, res, rtol=1e-4)
     end
end

@testset "test print_map_WindowFIntegrated_multipole" begin
     windowfint = GaPSE.WindowFIntegrated("datatest/IntegrF_REFERENCE_pi2_z115.txt")
     name_file_Qmultipoles = "datatest/PowerSpectrumGenWin/Ql_multipoles_of_F_REFERENCE_z115.txt"
     calc_name_file_Qmultipoles = "calc_Ql_multipoles_of_F_REFERENCE_z115.txt"
     isfile(calc_name_file_Qmultipoles) && rm(calc_name_file_Qmultipoles)

     GaPSE.print_map_WindowFIntegrated_multipole(windowfint, calc_name_file_Qmultipoles,
          "datatest/WideA_ZA_background.dat"; z_min=1.0, z_max=1.5, st=1.0, N=1000, pr=true, L_max=4,
          alg=:lobatto, N_lob=100, N_trap=200, atol_quad=0.0, rtol_quad=1e-2, enhancer=1e6, h_0=0.7)

     calc_ss, calc_all_Q_l1 = GaPSE.readxall(calc_name_file_Qmultipoles)
     ss, all_Q_l1 = GaPSE.readxall(name_file_Qmultipoles)

     @test all([isapprox(t, c, rtol=1e-4) for (t,c) in zip(ss, calc_ss)])
     for (qls,calc_qls) in zip(all_Q_l1, calc_all_Q_l1)
          @test all([isapprox(t, c, rtol=1e-4) for (t,c) in zip(qls,calc_qls) ])
     end

     rm(calc_name_file_Qmultipoles)
end
