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

module FFTLog

using FFTW
using Base: @kwdef
using SpecialFunctions: gamma
import Base: *

export SingleBesselPlan, HankelPlan
export prepare_FFTLog!, evaluate_FFTLog, evaluate_FFTLog!
export prepare_Hankel!, evaluate_Hankel, evaluate_Hankel!
export mul!, get_y


include("./FFTLog_files/common.jl")
include("./FFTLog_files/SingleBessel.jl")


##########################################################################################92


function mul!(Y, Q::SingleBesselPlan, A)
    evaluate_FFTLog!(Y, Q, A)
end

function mul!(Y, Q::HankelPlan, A)
    Y[:, :] .= evaluate_Hankel!(Y, Q, A)
end

function *(Q::SingleBesselPlan, A)
    evaluate_FFTLog(Q, A)
end

function *(Q::HankelPlan, A)
    evaluate_Hankel(Q, A)
end

end 
