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


@testset "TF" begin
RTOL = 1e-2
ks, pks = GaPSE.readxy("datatest/Tk.dat")
my_TF = GaPSE.TF(ks, pks)

@test my_TF(0.6e-1) ≈ 0.22902775397420208

ks_true, pks_true = GaPSE.readxy("datatest/PNG/Transfer_Function.dat")
pks_calc = [my_TF(k) for k in ks_true]

@test all([isapprox(a, b; rtol=RTOL) for (a, b) in zip(pks_true, pks_calc)])
end


@testset "alpha_bias" begin
RTOL = 1e-3
ks, pks = GaPSE.readxy("datatest/Tk.dat")
my_TF = GaPSE.TF(ks, pks)

@test isapprox(GaPSE.α_bias(1, my_TF;
        bf=1.0, D=1.0, Ω_M0=0.29992), 1.3942430956475444e-5; rtol=RTOL)
@test isapprox(GaPSE.α_bias(1, my_TF;
        bf=5.0, D=1.0, Ω_M0=0.29992), 6.971215478237722e-5; rtol=RTOL)
@test isapprox(GaPSE.α_bias(1, my_TF;
        bf=5.0, D=0.13, Ω_M0=0.29992), 0.0005362473444798248; rtol=RTOL)
@test isapprox(GaPSE.α_bias(1, my_TF;
        bf=5.0, D=0.13, Ω_M0=0.9), 0.001609171145745006; rtol=RTOL)
@test isapprox(GaPSE.α_bias(7.934, my_TF;
        bf=5.0, D=0.13, Ω_M0=0.9), 0.0008966355521073972; rtol=RTOL)

xs_true, ys_true = GaPSE.readxy("datatest/PNG/alpha_bias.dat")
ys_calc = [GaPSE.α_bias(x, my_TF;
    bf=1.0, D=1.0, Ω_M0=0.29992) for x in xs_true]

@test all([isapprox(a, b; rtol=RTOL) for (a, b) in zip(ys_true, ys_calc)])
end


@testset "IntegralIPSalpha" begin
RTOL = 1e-3
ks, pks = GaPSE.readxy("datatest/Tk.dat")
my_TF = GaPSE.TF(ks, pks)
my_IntIPSalpha = GaPSE.IntegralIPSalpha(my_TF, COSMO, 0, 0;
    D=nothing, bf=1.0, N=1024, kmin=1e-6, kmax=1e4, s0=1e-4,
    fit_left_min=nothing, fit_left_max=nothing, p0_left=nothing,
    fit_right_min=nothing, fit_right_max=nothing, p0_right=nothing)

@test isapprox(my_IntIPSalpha(1.34), 9.484203885952665e-5; rtol=RTOL)

xs_true, ys_true = GaPSE.readxy("datatest/PNG/Integral_IPS_alpha.dat")
ys_calc = [my_IntIPSalpha(x) for x in xs_true]

@test all([isapprox(a, b; rtol=RTOL) for (a, b) in zip(ys_true, ys_calc)])
end



##########################################################################################92

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


@testset "test CosmoPNGParams" begin
# The following D_eff is the linear growth factor in s_eff
# for 1.0 < z < 1.5, so s_eff = 2725.87720083235
D_eff = 0.5487092330211136

@testset "zeros" begin
    @test_throws AssertionError GaPSE.CosmoPNGParams(-1.0)

    @test_throws AssertionError GaPSE.CosmoPNGParams(D_eff; bf=-1.0)

    @test_throws AssertionError GaPSE.CosmoPNGParams(D_eff; flm_0=-1.0)
    @test_throws AssertionError GaPSE.CosmoPNGParams(D_eff; flM_0=-1.0)
    @test_throws AssertionError GaPSE.CosmoPNGParams(D_eff; flm_0=1.0, flM_0=0.5)
    @test_throws AssertionError GaPSE.CosmoPNGParams(D_eff; kmin_0=-1.0)
    @test_throws AssertionError GaPSE.CosmoPNGParams(D_eff; kmax_0=-1.0)
    @test_throws AssertionError GaPSE.CosmoPNGParams(D_eff; kmin_0=1.0, kmax_0=0.5)
    @test_throws AssertionError GaPSE.CosmoPNGParams(D_eff; s0_0=1.0, kmin_0=0.1, kmax_0=0.5)
    @test_throws AssertionError GaPSE.CosmoPNGParams(D_eff; N_0=3)

    @test_throws AssertionError GaPSE.CosmoPNGParams(D_eff; flm_2=-1.0)
    @test_throws AssertionError GaPSE.CosmoPNGParams(D_eff; flM_2=-1.0)
    @test_throws AssertionError GaPSE.CosmoPNGParams(D_eff; flm_2=1.0, flM_2=0.5)
    @test_throws AssertionError GaPSE.CosmoPNGParams(D_eff; kmin_2=-1.0)
    @test_throws AssertionError GaPSE.CosmoPNGParams(D_eff; kmax_2=-1.0)
    @test_throws AssertionError GaPSE.CosmoPNGParams(D_eff; kmin_2=1.0, kmax_2=0.5)
    @test_throws AssertionError GaPSE.CosmoPNGParams(D_eff; s0_2=1.0, kmin_2=0.1, kmax_2=0.5)
    @test_throws AssertionError GaPSE.CosmoPNGParams(D_eff; N_2=3)
end

@testset "first" begin
    bf = 1.0
    flm_0, flM_0, s0_0 = 5e-2, 1e-1, 1e-4
    kmin_0, kmax_0, N_0 = 1e-6, 1e4, 1024
    flm_2, flM_2, s0_2 = 5e-1, 1e0, 1e-4
    kmin_2, kmax_2, N_2 = 1e-6, 1e4, 1024

    params = GaPSE.CosmoPNGParams(D_eff;
        bf=bf,
        flm_0=flm_0, flM_0=flM_0, s0_0=s0_0,
        kmin_0=kmin_0, kmax_0=kmax_0, N_0=N_0,
        flm_2=flm_2, flM_2=flM_2, s0_2=s0_2,
        kmin_2=kmin_2, kmax_2=kmax_2, N_2=N_2
    )

    @test params.D ≈ D_eff
    @test params.bf ≈ bf

    @test params.flm_0 ≈ flm_0
    @test params.flM_0 ≈ flM_0
    @test params.s0_0 ≈ s0_0
    @test params.kmin_0 ≈ kmin_0
    @test params.kmax_0 ≈ kmax_0
    @test params.N_0 ≈ N_0

    @test params.flm_2 ≈ flm_2
    @test params.flM_2 ≈ flM_2
    @test params.s0_2 ≈ s0_2
    @test params.kmin_2 ≈ kmin_2
    @test params.kmax_2 ≈ kmax_2
    @test params.N_2 ≈ N_2
end

@testset "second" begin
    bf = 2.0
    flm_0, flM_0, s0_0 = 4e-2, 2e-1, 1e-3
    kmin_0, kmax_0, N_0 = 1e-5, 1e3, 24
    flm_2, flM_2, s0_2 = 1e-1, 1e1, 1e-3
    kmin_2, kmax_2, N_2 = 1e-8, 1e8, 2048

    params = GaPSE.CosmoPNGParams(D_eff;
        bf=bf,
        flm_0=flm_0, flM_0=flM_0, s0_0=s0_0,
        kmin_0=kmin_0, kmax_0=kmax_0, N_0=N_0,
        flm_2=flm_2, flM_2=flM_2, s0_2=s0_2,
        kmin_2=kmin_2, kmax_2=kmax_2, N_2=N_2
    )

    @test params.D ≈ D_eff
    @test params.bf ≈ bf

    @test params.flm_0 ≈ flm_0
    @test params.flM_0 ≈ flM_0
    @test params.s0_0 ≈ s0_0
    @test params.kmin_0 ≈ kmin_0
    @test params.kmax_0 ≈ kmax_0
    @test params.N_0 ≈ N_0

    @test params.flm_2 ≈ flm_2
    @test params.flM_2 ≈ flM_2
    @test params.s0_2 ≈ s0_2
    @test params.kmin_2 ≈ kmin_2
    @test params.kmax_2 ≈ kmax_2
    @test params.N_2 ≈ N_2
end
end




@testset "test_CosmoPNG" begin
RTOL = 1e-2

@testset "first" begin
    my_cosmopngparams = GaPSE.CosmoPNGParams(
        COSMO.D_of_s(COSMO.s_eff);
        bf=1.0,
        Dict(
            :flm_0 => 5e-2, :flM_0 => 1e-1, :s0_0 => 1e-4,
            :kmin_0 => 1e-6, :kmax_0 => 1e4, :N_0 => 1024,
            :flm_2 => 5e-1, :flM_2 => 1e0, :s0_2 => 1e-4,
            :kmin_2 => 1e-6, :kmax_2 => 1e4, :N_2 => 1024
        )...
    )

    my_cosmopng = GaPSE.CosmoPNG(my_cosmopngparams, COSMO, "datatest/Tk.dat")

    ks_true, pks_true = GaPSE.readxy("datatest/PNG/Transfer_Function.dat")
    my_TF = GaPSE.TF(ks_true, pks_true)

    IntIPSalpha_J0 = GaPSE.IntegralIPSalpha(my_TF, COSMO, 0, 0;
        D=nothing, bf=1.0, N=1024, kmin=1e-6, kmax=1e4, s0=1e-4,
        fit_left_min=nothing, fit_left_max=nothing, p0_left=nothing,
        fit_right_min=nothing, fit_right_max=nothing, p0_right=nothing)
    IntIPSalpha_J2 = GaPSE.IntegralIPSalpha(my_TF, COSMO, 2, 0;
        D=nothing, bf=1.0, N=1024, kmin=1e-6, kmax=1e4, s0=1e-4,
        fit_left_min=nothing, fit_left_max=nothing, p0_left=nothing,
        fit_right_min=nothing, fit_right_max=nothing, p0_right=nothing)

    xs = 10.0 .^ range(-4, 3, length=300)

    ys_1_true = [IntIPSalpha_J0(k) for k in xs]
    ys_2_true = [IntIPSalpha_J2(k) for k in xs]

    pks_calc = [my_cosmopng.tf(k) for k in ks_true]
    ys_1_calc = [my_cosmopng.J0(k) for k in xs]
    ys_2_calc = [my_cosmopng.J2(k) for k in xs]

    @test isapprox(IntIPSalpha_J0(1.34), 9.484203885952665e-5; rtol=RTOL)
    @test isapprox(my_cosmopng.J0(1.34), 9.484203885952665e-5; rtol=RTOL)
    @test isapprox(IntIPSalpha_J2(1.34), 1.689394925744043e-5; rtol=RTOL)
    @test isapprox(my_cosmopng.J2(1.34), 1.689394925744043e-5; rtol=RTOL)

    @test all([isapprox(a, b; rtol=RTOL) for (a, b) in zip(pks_true, pks_calc)])
    @test my_cosmopng.file_TF == "datatest/Tk.dat"
    @test all([isapprox(a, b; rtol=RTOL) for (a, b) in zip(ys_1_true, ys_1_calc)])
    @test all([isapprox(a, b; rtol=RTOL) for (a, b) in zip(ys_2_true, ys_2_calc)])

end

@testset "second" begin
    flm_0, flM_0 = 5e-3, 1e-2
    flm_2, flM_2 = 4e-3, 1e-1

    my_cosmopngparams = GaPSE.CosmoPNGParams(
        COSMO.D_of_s(COSMO.s_eff);
        bf=1.0,
        Dict(
            :flm_0 => flm_0, :flM_0 => flM_0, :s0_0 => 1e-4,
            :kmin_0 => 1e-6, :kmax_0 => 1e4, :N_0 => 1024,
            :flm_2 => flm_2, :flM_2 => flM_2, :s0_2 => 1e-4,
            :kmin_2 => 1e-6, :kmax_2 => 1e4, :N_2 => 1024
        )...
    )

    my_cosmopng = GaPSE.CosmoPNG(my_cosmopngparams, COSMO, "datatest/Tk.dat";
        comments=true)

    ks_true, pks_true = GaPSE.readxy("datatest/PNG/Transfer_Function.dat")
    my_TF = GaPSE.TF(ks_true, pks_true)

    IntIPSalpha_J0 = GaPSE.IntegralIPSalpha(my_TF, COSMO, 0, 0;
        D=nothing, bf=1.0, N=1024, kmin=1e-6, kmax=1e4, s0=1e-4,
        fit_left_min=flm_0, fit_left_max=flM_0, p0_left=nothing,
        fit_right_min=nothing, fit_right_max=nothing, p0_right=nothing)
    IntIPSalpha_J2 = GaPSE.IntegralIPSalpha(my_TF, COSMO, 2, 0;
        D=nothing, bf=1.0, N=1024, kmin=1e-6, kmax=1e4, s0=1e-4,
        fit_left_min=flm_2, fit_left_max=flM_2, p0_left=nothing,
        fit_right_min=nothing, fit_right_max=nothing, p0_right=nothing)


    xs = 10.0 .^ range(-6, 3, length=300)

    ys_1_true = [IntIPSalpha_J0(k) for k in xs]
    ys_2_true = [IntIPSalpha_J2(k) for k in xs]

    pks_calc = [my_cosmopng.tf(k) for k in ks_true]
    ys_1_calc = [my_cosmopng.J0(k) for k in xs]
    ys_2_calc = [my_cosmopng.J2(k) for k in xs]

    @test isapprox(IntIPSalpha_J0(1.34), 9.484203885952665e-5; rtol=RTOL)
    @test isapprox(my_cosmopng.J0(1.34), 9.484203885952665e-5; rtol=RTOL)
    @test isapprox(IntIPSalpha_J2(1.34), 1.689394925744043e-5; rtol=RTOL)
    @test isapprox(my_cosmopng.J2(1.34), 1.689394925744043e-5; rtol=RTOL)

    @test all([isapprox(a, b; rtol=RTOL) for (a, b) in zip(pks_true, pks_calc)])
    @test my_cosmopng.file_TF == "datatest/Tk.dat"
    @test all([isapprox(a, b; rtol=RTOL) for (a, b) in zip(ys_1_true, ys_1_calc)])
    @test all([isapprox(a, b; rtol=RTOL) for (a, b) in zip(ys_2_true, ys_2_calc)])
end

end

##########################################################################################92


@testset "test map_ξ_S_multipole" begin
RTOL = 1e-3
kwargs_xis_PP = Dict(
    :pr => false, :enhancer => 1e8,
    :atol_quad => 0.0, :rtol_quad => 1e-2,
    :N_log => 100,
)

cosmopngparams = GaPSE.CosmoPNGParams(
    COSMO.D_of_s(COSMO.s_eff);
    bf=1.0,
    Dict(
        :flm_0 => 5e-3, :flM_0 => 1e-2, :s0_0 => 1e-4,
        :kmin_0 => 1e-6, :kmax_0 => 1e4, :N_0 => 1024,
        :flm_2 => 5e-3, :flM_2 => 1e-2, :s0_2 => 1e-4,
        :kmin_2 => 1e-6, :kmax_2 => 1e4, :N_2 => 1024
    )...
)
cosmopng = GaPSE.CosmoPNG(cosmopngparams, COSMO, "datatest/Tk.dat")

@testset "with F" begin
    @testset "monopole" begin
        L = 0
        true_xi = "datatest/PNG/xi_s_withF_L$L" * ".txt"

        table = readdlm(true_xi; comments=true)
        ss = convert(Vector{Float64}, table[:, 1])
        xis = convert(Vector{Float64}, table[:, 2])

        calc_ss, calc_xis = GaPSE.map_ξ_S_multipole(COSMO, cosmopng,
            10 .^ range(0, 3, length=300); use_windows=true,
            L=L, kwargs_xis_PP...)

        @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(ss, calc_ss)])
        @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(xis, calc_xis)])
    end

    @testset "quadrupole" begin
        L = 2
        true_xi = "datatest/PNG/xi_s_withF_L$L" * ".txt"

        table = readdlm(true_xi; comments=true)
        ss = convert(Vector{Float64}, table[:, 1])
        xis = convert(Vector{Float64}, table[:, 2])

        calc_ss, calc_xis = GaPSE.map_ξ_S_multipole(COSMO, cosmopng,
            10 .^ range(0, 3, length=300); use_windows=true,
            L=L, kwargs_xis_PP...)

        @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(ss, calc_ss)])
        @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(xis, calc_xis)])
    end

end

@testset "without F" begin
    @testset "monopole" begin
        L = 0
        true_xi = "datatest/PNG/xi_s_noF_L$L" * ".txt"

        table = readdlm(true_xi; comments=true)
        ss = convert(Vector{Float64}, table[:, 1])
        xis = convert(Vector{Float64}, table[:, 2])

        calc_ss, calc_xis = GaPSE.map_ξ_S_multipole(COSMO, cosmopng,
            10 .^ range(0, 3, length=300); use_windows=false,
            L=L, kwargs_xis_PP...)

        @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(ss, calc_ss)])
        @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(xis, calc_xis)])
    end

    @testset "quadrupole" begin
        L = 2
        true_xi = "datatest/PNG/xi_s_noF_L$L" * ".txt"

        table = readdlm(true_xi; comments=true)
        ss = convert(Vector{Float64}, table[:, 1])
        xis = convert(Vector{Float64}, table[:, 2])

        calc_ss, calc_xis = GaPSE.map_ξ_S_multipole(COSMO, cosmopng,
            10 .^ range(0, 3, length=300); use_windows=false,
            L=L, kwargs_xis_PP...)

        @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(ss, calc_ss)])
        @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(xis, calc_xis)])
    end

end
end



##########################################################################################92



@testset "test print_map_ξ_S_multipole" begin
RTOL = 1e-3
kwargs_xis_PP = Dict(
    :pr => false, :enhancer => 1e8,
    :atol_quad => 0.0, :rtol_quad => 1e-2,
    :N_log => 100,
)
cosmopngparams = GaPSE.CosmoPNGParams(
    COSMO.D_of_s(COSMO.s_eff);
    bf=1.0,
    Dict(
        :flm_0 => 5e-3, :flM_0 => 1e-2, :s0_0 => 1e-4,
        :kmin_0 => 1e-6, :kmax_0 => 1e4, :N_0 => 1024,
        :flm_2 => 5e-3, :flM_2 => 1e-2, :s0_2 => 1e-4,
        :kmin_2 => 1e-6, :kmax_2 => 1e4, :N_2 => 1024
    )...
)
cosmopng = GaPSE.CosmoPNG(cosmopngparams, COSMO, "datatest/Tk.dat")

@testset "with F" begin
    @testset "monopole" begin
        L = 0
        name = "calc_xi_s_withF_L$L" * ".txt"
        true_xi = "datatest/PNG/xi_s_withF_L$L" * ".txt"

        isfile(name) && rm(name)

        table = readdlm(true_xi; comments=true)
        ss = convert(Vector{Float64}, table[:, 1])
        xis = convert(Vector{Float64}, table[:, 2])

        GaPSE.print_map_ξ_S_multipole(COSMO, cosmopng, name,
            10 .^ range(0, 3, length=300); use_windows=true,
            L=L, kwargs_xis_PP...)

        calc_table = readdlm(true_xi; comments=true)
        calc_ss = convert(Vector{Float64}, calc_table[:, 1])
        calc_xis = convert(Vector{Float64}, calc_table[:, 2])

        @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(ss, calc_ss)])
        @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(xis, calc_xis)])

        rm(name)
    end

    @testset "quadrupole" begin
        L = 2
        name = "calc_xi_s_withF_L$L" * ".txt"
        true_xi = "datatest/PNG/xi_s_withF_L$L" * ".txt"

        isfile(name) && rm(name)

        table = readdlm(true_xi; comments=true)
        ss = convert(Vector{Float64}, table[:, 1])
        xis = convert(Vector{Float64}, table[:, 2])

        GaPSE.print_map_ξ_S_multipole(COSMO, cosmopng, name,
            10 .^ range(0, 3, length=300); use_windows=true,
            L=L, kwargs_xis_PP...)

        calc_table = readdlm(true_xi; comments=true)
        calc_ss = convert(Vector{Float64}, calc_table[:, 1])
        calc_xis = convert(Vector{Float64}, calc_table[:, 2])

        @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(ss, calc_ss)])
        @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(xis, calc_xis)])

        rm(name)
    end

end

@testset "without F" begin
    @testset "monopole" begin
        L = 0
        name = "calc_xi_s_noF_L$L" * ".txt"
        true_xi = "datatest/PNG/xi_s_noF_L$L" * ".txt"

        isfile(name) && rm(name)

        table = readdlm(true_xi; comments=true)
        ss = convert(Vector{Float64}, table[:, 1])
        xis = convert(Vector{Float64}, table[:, 2])

        GaPSE.print_map_ξ_S_multipole(COSMO, cosmopng, name,
            10 .^ range(0, 3, length=300); use_windows=false,
            L=L, kwargs_xis_PP...)

        calc_table = readdlm(true_xi; comments=true)
        calc_ss = convert(Vector{Float64}, calc_table[:, 1])
        calc_xis = convert(Vector{Float64}, calc_table[:, 2])

        @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(ss, calc_ss)])
        @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(xis, calc_xis)])

        rm(name)
    end

    @testset "quadrupole" begin
        L = 2
        name = "calc_xi_s_noF_L$L" * ".txt"
        true_xi = "datatest/PNG/xi_s_noF_L$L" * ".txt"

        isfile(name) && rm(name)

        table = readdlm(true_xi; comments=true)
        ss = convert(Vector{Float64}, table[:, 1])
        xis = convert(Vector{Float64}, table[:, 2])

        GaPSE.print_map_ξ_S_multipole(COSMO, cosmopng, name,
            10 .^ range(0, 3, length=300); use_windows=false,
            L=L, kwargs_xis_PP...)

        calc_table = readdlm(true_xi; comments=true)
        calc_ss = convert(Vector{Float64}, calc_table[:, 1])
        calc_xis = convert(Vector{Float64}, calc_table[:, 2])

        @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(ss, calc_ss)])
        @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(xis, calc_xis)])

        rm(name)
    end

end
end
