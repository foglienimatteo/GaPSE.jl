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


@testset "test  specif_kwargs_LD" begin
    kwargs = Dict(:N_χs_2 => 200 , :en => 1e6, :N_χs=> 150)
    @test_throws AssertionError GaPSE.specif_kwargs_LD("ciao", kwargs)
    @test_throws AssertionError GaPSE.specif_kwargs_LD("auto_newton", kwargs)
    @test_throws AssertionError GaPSE.specif_kwargs_LD("newton_localgp", kwargs)

    @test GaPSE.specif_kwargs_LD("auto_doppler", kwargs) == Dict()
    @test GaPSE.specif_kwargs_LD("auto_lensing", kwargs) == Dict(:N_χs_2 => 200 , :en => 1e6)
    @test GaPSE.specif_kwargs_LD("localgp_lensing", kwargs) == Dict(:N_χs => 150 , :en => 1e6)
end

@testset "test  specif_kwargs_GNC" begin
    kwargs = Dict(:N_lob => 20, :N_χs_2 => 200 , :en => 1e6, :N_χs => 150, :suit_sampling => false)
    @test_throws AssertionError GaPSE.specif_kwargs_GNC("ciao", kwargs)
    @test_throws AssertionError GaPSE.specif_kwargs_GNC("automatic_newton", kwargs)
    @test_throws AssertionError GaPSE.specif_kwargs_GNC("newtonian_localgp", kwargs)

    @test GaPSE.specif_kwargs_GNC("auto_newton", kwargs) == Dict(:N_lob => 20)
    @test GaPSE.specif_kwargs_GNC("auto_doppler", kwargs) == Dict(:N_lob => 20)
    @test GaPSE.specif_kwargs_GNC("auto_lensing", kwargs) == Dict(:N_lob => 20,:N_χs_2 => 200, :en => 1e6, :suit_sampling => false)
    @test GaPSE.specif_kwargs_GNC("auto_localgp", kwargs) == Dict(:N_lob => 20)
    @test GaPSE.specif_kwargs_GNC("auto_integratedgp", kwargs) == Dict(:N_lob => 20, :N_χs_2 => 200, :en => 1e6, :suit_sampling => false)

    @test GaPSE.specif_kwargs_GNC("newton_doppler", kwargs) == Dict(:N_lob => 20)
    @test GaPSE.specif_kwargs_GNC("newton_lensing", kwargs) == Dict(:N_lob => 20, :N_χs => 150, :en => 1e6, :suit_sampling => false)
    @test GaPSE.specif_kwargs_GNC("newton_localgp", kwargs) == Dict(:N_lob => 20)
    @test GaPSE.specif_kwargs_GNC("newton_integratedgp", kwargs) == Dict(:N_lob => 20, :N_χs => 150, :en => 1e6, :suit_sampling => false)

    @test GaPSE.specif_kwargs_GNC("doppler_newton", kwargs) == Dict(:N_lob => 20)
    @test GaPSE.specif_kwargs_GNC("doppler_lensing", kwargs) == Dict(:N_lob => 20, :N_χs => 150, :en => 1e6, :suit_sampling => false)
    @test GaPSE.specif_kwargs_GNC("doppler_localgp", kwargs) == Dict(:N_lob => 20)
    @test GaPSE.specif_kwargs_GNC("doppler_integratedgp", kwargs) == Dict(:N_lob => 20, :N_χs => 150, :en => 1e6, :suit_sampling => false)

    @test GaPSE.specif_kwargs_GNC("lensing_newton", kwargs) == Dict(:N_lob => 20, :N_χs => 150, :en => 1e6, :suit_sampling => false)
    @test GaPSE.specif_kwargs_GNC("lensing_doppler", kwargs) == Dict(:N_lob => 20, :N_χs => 150, :en => 1e6, :suit_sampling => false)
    @test GaPSE.specif_kwargs_GNC("lensing_localgp", kwargs) == Dict(:N_lob => 20, :N_χs => 150, :en => 1e6, :suit_sampling => false)
    @test GaPSE.specif_kwargs_GNC("lensing_integratedgp", kwargs) == Dict(:N_lob => 20, :N_χs_2 => 200, :en => 1e6, :suit_sampling => false)

    @test GaPSE.specif_kwargs_GNC("localgp_newton", kwargs) == Dict(:N_lob => 20)
    @test GaPSE.specif_kwargs_GNC("localgp_doppler", kwargs) == Dict(:N_lob => 20)
    @test GaPSE.specif_kwargs_GNC("localgp_lensing", kwargs) == Dict(:N_lob => 20, :N_χs => 150, :en => 1e6, :suit_sampling => false)
    @test GaPSE.specif_kwargs_GNC("localgp_integratedgp", kwargs) == Dict(:N_lob => 20, :N_χs => 150, :en => 1e6, :suit_sampling => false)

    @test GaPSE.specif_kwargs_GNC("integratedgp_newton", kwargs) == Dict(:N_lob => 20, :N_χs => 150, :en => 1e6, :suit_sampling => false)
    @test GaPSE.specif_kwargs_GNC("integratedgp_doppler", kwargs) == Dict(:N_lob => 20, :N_χs => 150, :en => 1e6, :suit_sampling => false)
    @test GaPSE.specif_kwargs_GNC("integratedgp_lensing", kwargs) == Dict(:N_lob => 20, :N_χs_2 => 200, :en => 1e6, :suit_sampling => false)
    @test GaPSE.specif_kwargs_GNC("integratedgp_localgp", kwargs) == Dict(:N_lob => 20, :N_χs => 150, :en => 1e6, :suit_sampling => false)
end

@testset "test  specif_kwargs_GNCxLD" begin
    kwargs = Dict(:N_χs_2 => 200 , :en => 1e6, :N_χs=> 150)
    @test_throws AssertionError GaPSE.specif_kwargs_GNCxLD("ciao", kwargs)
    @test_throws AssertionError GaPSE.specif_kwargs_GNCxLD("auto_newton", kwargs)
    @test_throws AssertionError GaPSE.specif_kwargs_GNCxLD("newtonian_localgp", kwargs)

    @test GaPSE.specif_kwargs_GNCxLD("doppler_doppler", kwargs) == Dict()
    @test GaPSE.specif_kwargs_GNCxLD("lensing_lensing", kwargs) == Dict(:N_χs_2 => 200 , :en => 1e6)
    @test GaPSE.specif_kwargs_GNCxLD("integratedgp_lensing", kwargs) == Dict(:N_χs_2 => 200 , :en => 1e6)
    @test GaPSE.specif_kwargs_GNCxLD("localgp_lensing", kwargs) == Dict(:N_χs => 150 , :en => 1e6)
    @test GaPSE.specif_kwargs_GNCxLD("newton_integratedgp", kwargs) == Dict(:N_χs => 150 , :en => 1e6)
end

@testset "test  specif_kwargs_LDxGNC" begin
    kwargs = Dict(:N_χs_2 => 200 , :en => 1e6, :N_χs=> 150)
    @test_throws AssertionError GaPSE.specif_kwargs_LDxGNC("ciao", kwargs)
    @test_throws AssertionError GaPSE.specif_kwargs_LDxGNC("auto_newton", kwargs)
    @test_throws AssertionError GaPSE.specif_kwargs_LDxGNC("newtonian_localgp", kwargs)

    @test GaPSE.specif_kwargs_LDxGNC("localgp_doppler", kwargs) == Dict()
    @test GaPSE.specif_kwargs_LDxGNC("lensing_lensing", kwargs) == Dict(:N_χs_2 => 200 , :en => 1e6)
    @test GaPSE.specif_kwargs_LDxGNC("integratedgp_lensing", kwargs) == Dict(:N_χs_2 => 200 , :en => 1e6)
    @test GaPSE.specif_kwargs_LDxGNC("localgp_lensing", kwargs) == Dict(:N_χs => 150 , :en => 1e6)
    @test GaPSE.specif_kwargs_LDxGNC("integratedgp_newton", kwargs) == Dict(:N_χs => 150 , :en => 1e6)
end

