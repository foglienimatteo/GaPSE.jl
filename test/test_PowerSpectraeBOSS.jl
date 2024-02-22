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


@testset "test print_PS_multipole_GenWin with eBOSS" begin
    RTOL = 1e-2

    name = "datatest/PowerSpectrumeBOSS/PS_eBOSS_until_L8.txt"
    name_file_eBOSS_multipoles = "datatest/PowerSpectrumeBOSS/eBOSS_window_NGC.txt"
    calc_name = "calc_ps_eboss_L8.txt"
    isfile(calc_name) && rm(calc_name)

    ps_kwargs(alg::Symbol = :fftlog) = alg == :twofast ?
    Dict(
        :alg => :twofast, :epl => true, :pr => false, 
        :N_left => 12, :N_right => 12,
        :p0_left => [-2.0, 1.0], :p0_right => [-2.0, 1.0],
        :int_s_min => 1e0, :int_s_max => 1200.0,
        :cut_first_n=>0, :cut_last_n => 0
    ) : alg == :fftlog ?
    Dict(
        :alg => :fftlog, :pr=>true, :ν => 1.5, 
        :n_extrap_low => 0, :n_extrap_high => 0, 
        :n_pad => 500, :cut_first_n=>0, :cut_last_n => 0,
    ) : throw(AssertionError("alg = :fftlog (recommended) or alg = :twofast !"));
    tf = :fftlog;

    xi_filenames = "datatest/PowerSpectrumeBOSS/" .* [
        "all_xis_GNC_L0_noF_noobsvel_eBOSS.txt",
        "all_xis_GNC_L1_noF_noobsvel_eBOSS.txt",
        "all_xis_GNC_L2_noF_noobsvel_eBOSS.txt",
        "all_xis_GNC_L3_noF_noobsvel_eBOSS.txt",
        "all_xis_GNC_L4_noF_noobsvel_eBOSS.txt",
        "all_xis_GNC_L5_noF_noobsvel_eBOSS.txt",
        "all_xis_GNC_L6_noF_noobsvel_eBOSS.txt",
        "all_xis_GNC_L7_noF_noobsvel_eBOSS.txt",
        "all_xis_GNC_L8_noF_noobsvel_eBOSS.txt"
    ]

    calc_name_file_ximultipoles = "calc_xis_eboss_L012345678.txt"
    isfile(calc_name_file_ximultipoles) && rm(calc_name_file_ximultipoles)
    GaPSE.create_file_for_XiMultipoles(calc_name_file_ximultipoles, xi_filenames, 2, "GNC")

    isfile(calc_name) && rm(calc_name)
    GaPSE.print_PS_multipole_GenWin(calc_name_file_ximultipoles, name_file_eBOSS_multipoles, calc_name;
        L=0, ps_kwargs(:fftlog)...)
    calc_ks_fftlog, calc_pks_fftlog = GaPSE.readxy(calc_name)

    ks_fftlog, pks_fftlog = GaPSE.readxy(name)
    @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(calc_ks_fftlog, ks_fftlog)])
    @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(calc_pks_fftlog, pks_fftlog)])

    rm(calc_name)
    rm(calc_name_file_ximultipoles)
end


