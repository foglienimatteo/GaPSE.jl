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


@testset "test create_file_for_XiMultipoles" begin
     RTOL = 1e-5
     
     name = "test/datatest/PowerSpectrumGenWin/xis_LD_L0123_lob_noF_specific_ss.txt"
     calc_name = "spec_ss.txt"
     names = "datatest/PowerSpectrumGenWin/" .* [
          "xis_LD_L0_lob_noF_specific_ss.txt",
          "xis_LD_L1_lob_noF_specific_ss.txt",
          "xis_LD_L2_lob_noF_specific_ss.txt",
          "xis_LD_L3_lob_noF_specific_ss.txt",
          ]

     GaPSE.create_file_for_XiMultipoles(calc_name, names, "auto_doppler", "LD")
     calc_ss, calc_all_xis = GaPSE.readxall(name)

     ss, all_xis = GaPSE.readxall(name)
     @test all(isapprox(t,c; rtol = RTOL) for (t,c) in zip(ss, calc_ss))
     @test all(isapprox(t,c; rtol = RTOL) for (t,c) in zip(all_xis, calc_all_xis))

end
