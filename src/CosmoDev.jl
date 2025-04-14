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

##########################################################################################92



struct DevMySpline{A,B,C}
    xs::A       #Vector{Float64}
    coeffs::B   #Vector{Tuple{Float64,Float64,Float64,Float64}}
    N::C        #Int64
end

Adapt.@adapt_structure DevMySpline



##########################################################################################92



struct DevInputPS{F,S}
    l_si::F     #Float64
    l_b::F     #Float64
    l_a::F     #Float64
    left::F     #Float64

    spline::S    #GaPSE.MySpline

    r_si::F     #Float64
    r_b::F     #Float64
    r_a::F     #Float64
    right::F     #Float64
end

Adapt.@adapt_structure DevInputPS


struct DevIntegralIPS{F,S}
    l_si::F     #Float64
    l_b::F     #Float64
    l_a::F     #Float64
    left::F     #Float64

    spline::S    #GaPSE.MySpline

    r_si::F     #Float64
    r_b::F     #Float64
    r_a::F     #Float64
    right::F     #Float64
end

Adapt.@adapt_structure DevIntegralIPS


struct DevIPSTools{F, IPS}
    I00::IPS     #IntegralIPS
    I20::IPS     #IntegralIPS
    I40::IPS     #IntegralIPS
    I02::IPS     #IntegralIPS
    I22::IPS     #IntegralIPS
    I31::IPS     #IntegralIPS
    I13::IPS     #IntegralIPS
    I11::IPS     #IntegralIPS

    I04_tilde::IPS     #IntegralIPS

    σ_0::F     #Float64
    σ_1::F     #Float64
    σ_2::F     #Float64
    σ_3::F     #Float64
    σ_4::F     #Float64

    fit_min::F     #Float64
    fit_max::F     #Float64
    k_min::F     #Float64
    k_max::F     #Float64

end

Adapt.@adapt_structure DevIPSTools



##########################################################################################92



struct DevWindowF{V,M}
    xs::V     #Vector{Float64}
    μs::V     #Vector{Float64}
    Fs::M     #Matrix{Float64}
end

Adapt.@adapt_structure DevWindowF

struct DevWindowFIntegrated{V,M}
    ss::V     #Vector{Float64}
    μs::V     #Vector{Float64}
    IFs::M     #Matrix{Float64}
end

Adapt.@adapt_structure DevWindowFIntegrated

###

struct DevCosmoParams{F}
    z_min::F     #Float64
    z_max::F     #Float64
    θ_max::F     #Float64

    Ω_b::F     #Float64
    Ω_cdm::F     #Float64
    Ω_M0::F     #Float64
    h_0::F     #Float64

    b1::F     #Float64
    b2::F     #Float64
    s_b1::F     #Float64
    s_b2::F     #Float64
    𝑓_evo1::F     #Float64
    𝑓_evo2::F     #Float64

    s_lim::F     #Float64
    z_spline_lim::F     #Float64

    #IPS::Dict{Symbol,T1} where {T1}
    #IPSTools::Dict{Symbol,T2} where {T2}
    #WFI::Dict{Symbol,T3} where {T3}
end

Adapt.@adapt_structure DevCosmoParams



##########################################################################################92



struct DevCosmology{S,F}
    IPS::InputPS
    #ξ_matter::EPLs
    params::DevCosmoParams #CosmoParams
    tools::IPSTools
    windowF::WindowF
    windowFint::WindowFIntegrated
    WFI_norm::F     #Float64

    z_of_s::S     #GaPSE.MySpline
    D_of_s::S     #GaPSE.MySpline
    f_of_s::S     #GaPSE.MySpline
    ℋ_of_s::S     #GaPSE.MySpline
    ℋ_p_of_s::S     #GaPSE.MySpline
    ℛ_LD_of_s::S     #GaPSE.MySpline
    ℛ_GNC1_of_s::S     #GaPSE.MySpline
    ℛ_GNC2_of_s::S     #GaPSE.MySpline

    s_of_z::S     #GaPSE.MySpline

    z_eff::F     #Float64
    s_min::F     #Float64
    s_max::F     #Float64
    s_eff::F     #Float64
    s_spline_lim::F     #Float64

    volume::F     #Float64

    #file_data::String
    #file_ips::String
    #file_windowF::String
    #file_IWF::Union{String,Nothing}
end

Adapt.@adapt_structure DevCosmology

