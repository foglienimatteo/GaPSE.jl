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

#run(`wget https://zenodo.org/record/6021744/files/test_FFTLog.tar.xz\?download=1`)
#run(`mv test_FFTLog.tar.xz\?download=1 test_FFTLog.tar.xz`)
#run(`tar xvf test_FFTLog.tar.xz`)

#input_path = pwd()*"/test_FFTLog/"
input_path = "./datatest/FFTLog/"

read_file = readdlm(input_path * "Pk_test")
k = read_file[:, 1]
Pk = read_file[:, 2]

function f(k, a=1)
    return exp(-k^2.0 * a^2 / 2)
end

function F(r, a=1)
    return exp(-r^2.0 / (2.0 * a^2))
end


function LogSpaced(min::T, max::T, n::I) where {T,I}
    logmin = log10(min)
    logmax = log10(max)
    logarray = Array(LinRange(logmin, logmax, n))
    return exp10.(logarray)
end

#=
@testset "BeyondLimber checks" begin
    FFTTest = GaPSE.FFTLog.SingleBesselPlan(x=k, ν=1.01, n_extrap_low=1500, n_extrap_high=1500,
        n_pad=2000)
    Ell = Array([1.0, 2.0])
    GaPSE.FFTLog.prepare_FFTLog!(FFTTest, Ell)
    Y = GaPSE.FFTLog.get_y(FFTTest)
    FY = GaPSE.FFTLog.evaluate_FFTLog(FFTTest, Pk)
    Fr = npzread(input_path * "check_by_Fr.py.npy")
    r = npzread(input_path * "check_by_r.py.npy")

    @test isapprox(Fr, FY[1, :], rtol=1e-5)
    @test isapprox(r, Y[1, :], rtol=1e-5)

    FFTTest = GaPSE.FFTLog.SingleBesselPlan(x=k, ν=1.01, n_extrap_low=1500, n_extrap_high=1500,
        n_pad=2000, n=1)
    Ell = Array([1.0, 2.0])
    GaPSE.FFTLog.prepare_FFTLog!(FFTTest, Ell)
    Y = GaPSE.FFTLog.get_y(FFTTest)
    FY = GaPSE.FFTLog.evaluate_FFTLog(FFTTest, Pk)
    Fr1 = npzread(input_path * "check_by_Fr1.py.npy")
    r1 = npzread(input_path * "check_by_r1.py.npy")

    @test isapprox(Fr1, FY[1, :], rtol=1e-8)
    @test isapprox(r1, Y[1, :], rtol=1e-8)

    FFTTest = GaPSE.FFTLog.SingleBesselPlan(x=k, ν=1.01, n_extrap_low=1500, n_extrap_high=1500,
        n_pad=2000, n=2)
    Ell = Array([1.0, 2.0])
    GaPSE.FFTLog.prepare_FFTLog!(FFTTest, Ell)
    Y = GaPSE.FFTLog.get_y(FFTTest)
    FY = GaPSE.FFTLog.evaluate_FFTLog(FFTTest, Pk)
    Fr2 = npzread(input_path * "check_by_Fr2.py.npy")
    r2 = npzread(input_path * "check_by_r2.py.npy")

    @test isapprox(Fr2, FY[1, :], rtol=1e-8)
    @test isapprox(r2, Y[1, :], rtol=1e-8)

    FFTTest = GaPSE.FFTLog.HankelPlan(x=k, n_extrap_low=1500, ν=1.01, n_extrap_high=1500,
        n_pad=0)
    Ell = Array([0.0, 2.0])
    GaPSE.FFTLog.prepare_Hankel!(FFTTest, Ell)
    Y = GaPSE.FFTLog.get_y(FFTTest)
    FY = GaPSE.FFTLog.evaluate_Hankel(FFTTest, Pk .* (k .^ (-2)))
    Fr = npzread(input_path * "check_by_hank_Fr.py.npy")
    r = npzread(input_path * "check_by_hank_r.py.npy")
    @test isapprox(Fr, FY[1, :], rtol=1e-8)
    @test isapprox(r, Y[1, :], rtol=1e-8)
end
=#


@testset "Analytical Hankel test" begin
    GaPSE.FFTLog.set_num_threads(Threads.nthreads())
    k = LogSpaced(10^(-5), 10.0, 2^10)
    fk = f.(k)
    FFTTest = GaPSE.FFTLog.HankelPlan(x=k, n_extrap_low=1500, ν=1.01, n_extrap_high=1500,
        n_pad=500)
    Ell = Array([0.0])
    GaPSE.FFTLog.prepare_Hankel!(FFTTest, Ell)
    Y = GaPSE.FFTLog.get_y(FFTTest)
    FY = zeros(size(Y))
    FY = GaPSE.FFTLog.evaluate_Hankel(FFTTest, fk)
    FY_mul = zeros(size(FY))
    Fr_analytical = F.(Y[1, :])
    @test isapprox(Fr_analytical, FY[1, :], rtol=1e-10)
    @test isapprox(FY, GaPSE.FFTLog.mul!(FY_mul, FFTTest, fk), rtol=1e-10)
    @test isapprox(FY, FFTTest * fk, rtol=1e-10)
end

@testset "mul! operator test" begin
    k = LogSpaced(10^(-5), 10.0, 2^10)
    fk = f.(k)
    FFTTest = GaPSE.FFTLog.SingleBesselPlan(x=k, n_extrap_low=1500, ν=1.01, n_extrap_high=1500,
        n_pad=500)
    Ell = Array([0.0])
    GaPSE.FFTLog.prepare_FFTLog!(FFTTest, Ell)
    Y = GaPSE.FFTLog.get_y(FFTTest)
    FY = zeros(size(Y))
    fk_k2 = fk .* (k .^ 2)
    FY = GaPSE.FFTLog.evaluate_FFTLog(FFTTest, fk_k2)
    FY_mul = zeros(size(FY))
    GaPSE.FFTLog.mul!(FY_mul, FFTTest, fk_k2)
    @test isapprox(FY_mul, FY, rtol=1e-10)
    @test isapprox(FY_mul, FFTTest * fk_k2, rtol=1e-10)
end

