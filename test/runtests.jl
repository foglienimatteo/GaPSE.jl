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


using GaPSE, Test
using Dierckx, DelimitedFiles


const TEST_FILE = "TEST_DATA.txt"
include("TEST_DATA.jl")

##########################################################################################92

@testset "test_MathUtils" begin
     include("test_MathUtils.jl")
end

@testset "test_IPSTools" begin
     include("test_IPSTools.jl")
end

@testset "test_BackgroundData" begin
     include("test_BackgroundData.jl")
end

@testset "test_PowerSpectrum" begin
     include("test_PowerSpectrum.jl")
end

@testset "test_AutoDoppler" begin
     include("test_AutoDoppler.jl")
end


