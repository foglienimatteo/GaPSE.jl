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


struct DevMySpline{V,VT,I}
    xs::V        #Vector{Float64}
    coeffs::VT   #Vector{Tuple{Float64,Float64,Float64,Float64}}
    N::I         #Int64
end

function gpu_searchsortedlast(xs, x)
    lo = 1
    hi = length(xs)
    while lo < hi
        mid = (lo + hi + 1) >>> 1
        # direct device indexing
        if xs[mid] <= x
            lo = mid
        else
            hi = mid - 1
        end
    end
    return lo
end

function (S::DevMySpline)(x)
    #@assert S.xs[1] ≤ x ≤ S.xs[end] "BC Error: $(S.xs[1]) ≤ $x ≤ $(S.xs[end]) does not hold!"
    #i = searchsortedlast(S.xs, x)
    i = gpu_searchsortedlast(S.xs, x)
    #u = x - S.xs[i]
    u, a, b, c, d = (i == length(S.xs)) ? (x - S.xs[i-1], S.coeffs[i-1]...) : (x - S.xs[i], S.coeffs[i]...)
    (u ≈ zero(u)) && (return S.coeffs[i][1])

    return a + b * u + c * u^2 + d * u^3
    #return @evalpoly(u, S.coeffs[i]...)
end


# XA: Works with XB
#function Adapt.adapt_structure(to, s::MySpline)
#    DevMySpline(
#        adapt(to, devfloat.(s.xs)), 
#        adapt(to, [devfloat.(T) for T in s.coeffs]), 
#        adapt(to, s.N)
#    )
#end

function Adapt.adapt_structure(to, s::MySpline; devfloat=DevFloat)
    if devfloat == typeof(s.xs[begin])
        return DevMySpline(
            adapt(to, s.xs),
            adapt(to, s.coeffs),
            adapt(to, s.N)
        )
    else
        return DevMySpline(
            adapt(to, devfloat.(s.xs)),
            adapt(to, [devfloat.(T) for T in s.coeffs]),
            adapt(to, s.N)
        )
    end
end


function Adapt.adapt_structure(to, s::DevMySpline)
    DevMySpline(adapt(to, s.xs), adapt(to, s.coeffs), adapt(to, s.N))
end



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

function (IPS::DevInputPS)(x)
    if x < IPS.left
        return power_law(x, IPS.l_si, IPS.l_b, IPS.l_a)
    elseif x > IPS.right
        return power_law(x, IPS.r_si, IPS.r_b, IPS.r_a)
    else
        return IPS.spline(x)
    end
end

function Adapt.adapt_structure(to, s::InputPS; devfloat=DevFloat)
    if devfloat == typeof(s.l_si)
        return DevInputPS(
            adapt(to, s.l_si), adapt(to, s.l_b), adapt(to, s.l_a), adapt(to, s.left),
            Adapt.adapt_structure(to, s.spline; devfloat=devfloat),
            adapt(to, s.r_si), adapt(to, s.r_b), adapt(to, s.r_a), adapt(to, s.right),
        )
    else
        return DevInputPS(
            adapt(to, devfloat(s.l_si)), adapt(to, devfloat(s.l_b)), adapt(to, devfloat(s.l_a)), adapt(to, devfloat(s.left)),
            Adapt.adapt_structure(to, s.spline; devfloat=devfloat),
            adapt(to, devfloat(s.r_si)), adapt(to, devfloat(s.r_b)), adapt(to, devfloat(s.r_a)), adapt(to, devfloat(s.right)),
        )
    end
end


function Adapt.adapt_structure(to, s::DevInputPS)
    return DevInputPS(
        adapt(to, s.l_si), adapt(to, s.l_b), adapt(to, s.l_a), adapt(to, s.left),
        Adapt.adapt_structure(to, s.spline),
        adapt(to, s.r_si), adapt(to, s.r_b), adapt(to, s.r_a), adapt(to, s.right),
    )
end



##########################################################################################92



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

function (IPS::DevIntegralIPS)(x)
    if x < IPS.left
        return power_law(x, IPS.l_si, IPS.l_b, IPS.l_a)
    elseif x > IPS.right
        return power_law(x, IPS.r_si, IPS.r_b, IPS.r_a)
    else
        return IPS.spline(x)
    end
end

function Adapt.adapt_structure(to, s::IntegralIPS; devfloat=DevFloat)
    if devfloat == typeof(s.l_si)
        return DevIntegralIPS(
            adapt(to, s.l_si), adapt(to, s.l_b), adapt(to, s.l_a), adapt(to, s.left),
            Adapt.adapt_structure(to, s.spline; devfloat=devfloat),
            adapt(to, s.r_si), adapt(to, s.r_b), adapt(to, s.r_a), adapt(to, s.right),
        )
    else
        return DevIntegralIPS(
            adapt(to, devfloat(s.l_si)), adapt(to, devfloat(s.l_b)), adapt(to, devfloat(s.l_a)), adapt(to, devfloat(s.left)),
            Adapt.adapt_structure(to, s.spline; devfloat=devfloat),
            adapt(to, devfloat(s.r_si)), adapt(to, devfloat(s.r_b)), adapt(to, devfloat(s.r_a)), adapt(to, devfloat(s.right)),
        )
    end
end


function Adapt.adapt_structure(to, s::DevIntegralIPS)
    return DevIntegralIPS(
        adapt(to, s.l_si), adapt(to, s.l_b), adapt(to, s.l_a), adapt(to, s.left),
        Adapt.adapt_structure(to, s.spline),
        adapt(to, s.r_si), adapt(to, s.r_b), adapt(to, s.r_a), adapt(to, s.right),
    )
end


##########################################################################################92


struct DevIPSTools{F,IntIPS}
    I00::IntIPS     #IntegralIPS
    I20::IntIPS     #IntegralIPS
    I40::IntIPS     #IntegralIPS
    I02::IntIPS     #IntegralIPS
    I22::IntIPS     #IntegralIPS
    I31::IntIPS     #IntegralIPS
    I13::IntIPS     #IntegralIPS
    I11::IntIPS     #IntegralIPS

    I04_tilde::IntIPS     #IntegralIPS

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



function Adapt.adapt_structure(to, s::IPSTools; devfloat=DevFloat)
    if devfloat == typeof(s.fit_min)
        return DevIPSTools(
            Adapt.adapt_structure(to, s.I00; devfloat=devfloat),
            Adapt.adapt_structure(to, s.I20; devfloat=devfloat),
            Adapt.adapt_structure(to, s.I40; devfloat=devfloat),
            Adapt.adapt_structure(to, s.I02; devfloat=devfloat),
            Adapt.adapt_structure(to, s.I22; devfloat=devfloat),
            Adapt.adapt_structure(to, s.I31; devfloat=devfloat),
            Adapt.adapt_structure(to, s.I13; devfloat=devfloat),
            Adapt.adapt_structure(to, s.I11; devfloat=devfloat),
            
            Adapt.adapt_structure(to, s.I04_tilde; devfloat=devfloat),

            adapt(to, s.σ_0), adapt(to, s.σ_1), adapt(to, s.σ_2), 
            adapt(to, s.σ_3), adapt(to, s.σ_4),

            adapt(to, s.fit_min), adapt(to, s.fit_max), 
            adapt(to, s.k_min), adapt(to, s.k_max),
        )
    else
        return DevIPSTools(
            Adapt.adapt_structure(to, s.I00; devfloat=devfloat),
            Adapt.adapt_structure(to, s.I20; devfloat=devfloat),
            Adapt.adapt_structure(to, s.I40; devfloat=devfloat),
            Adapt.adapt_structure(to, s.I02; devfloat=devfloat),
            Adapt.adapt_structure(to, s.I22; devfloat=devfloat),
            Adapt.adapt_structure(to, s.I31; devfloat=devfloat),
            Adapt.adapt_structure(to, s.I13; devfloat=devfloat),
            Adapt.adapt_structure(to, s.I11; devfloat=devfloat),
            
            Adapt.adapt_structure(to, s.I04_tilde; devfloat=devfloat),

            adapt(to, devfloat(s.σ_0)), adapt(to, devfloat(s.σ_1)), adapt(to, devfloat(s.σ_2)), 
            adapt(to, devfloat(s.σ_3)), adapt(to, devfloat(s.σ_4)),
            adapt(to, devfloat(s.fit_min)), adapt(to, devfloat(s.fit_max)), 
            adapt(to, devfloat(s.k_min)), adapt(to, devfloat(s.k_max)),
        )
    end
end


function Adapt.adapt_structure(to, s::DevIPSTools)
    return DevIPSTools(
        Adapt.adapt_structure(to, s.I00),
        Adapt.adapt_structure(to, s.I20),
        Adapt.adapt_structure(to, s.I40),
        Adapt.adapt_structure(to, s.I02),
        Adapt.adapt_structure(to, s.I22),
        Adapt.adapt_structure(to, s.I31),
        Adapt.adapt_structure(to, s.I13),
        Adapt.adapt_structure(to, s.I11),

        Adapt.adapt_structure(to, s.I04_tilde),

        adapt(to, s.σ_0), adapt(to, s.σ_1), adapt(to, s.σ_2), 
        adapt(to, s.σ_3), adapt(to, s.σ_4),

        adapt(to, s.fit_min), adapt(to, s.fit_max), 
        adapt(to, s.k_min), adapt(to, s.k_max),
    )
end



##########################################################################################92


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

function Adapt.adapt_structure(to, s::CosmoParams; devfloat=DevFloat)
    if devfloat == typeof(s.z_min)
        return DevCosmoParams(
            s.z_min, s.z_max, s.θ_max, s.Ω_b, s.Ω_cdm, s.Ω_M0, s.h_0,
            s.b1, s.b2, s.s_b1, s.s_b2, s.𝑓_evo1, s.𝑓_evo2,
            s.s_lim, s.z_spline_lim
        )
    else
        return DevCosmoParams(
            adapt(to, devfloat(s.z_min)),
            adapt(to, devfloat(s.z_max)),
            adapt(to, devfloat(s.θ_max)),
            adapt(to, devfloat(s.Ω_b)),
            adapt(to, devfloat(s.Ω_cdm)),
            adapt(to, devfloat(s.Ω_M0)),
            adapt(to, devfloat(s.h_0)),
            adapt(to, devfloat(s.b1)),
            adapt(to, devfloat(s.b2)),
            adapt(to, devfloat(s.s_b1)),
            adapt(to, devfloat(s.s_b2)),
            adapt(to, devfloat(s.𝑓_evo1)),
            adapt(to, devfloat(s.𝑓_evo2)),
            adapt(to, devfloat(s.s_lim)),
            adapt(to, devfloat(s.z_spline_lim)),
        )
    end
end

Adapt.@adapt_structure DevCosmoParams



##########################################################################################92



struct DevCosmology{IPS,CP,IPST,S,F}
    IPS::IPS    #DevInputPS
    #ξ_matter::EPLs
    params::CP  #CosmoParams
    tools::IPST  #IPSTools
    #windowF::WindowF
    #windowFint::WindowFIntegrated
    #WFI_norm::F     #Float64

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



function Adapt.adapt_structure(to, s::Cosmology; devfloat=DevFloat)
    if devfloat == typeof(s.z_eff)
        return DevCosmology(
            Adapt.adapt_structure(to, s.IPS),
            Adapt.adapt_structure(to, s.params),
            Adapt.adapt_structure(to, s.tools), Adapt.adapt_structure(to, s.z_of_s),
            Adapt.adapt_structure(to, s.D_of_s),
            Adapt.adapt_structure(to, s.f_of_s),
            Adapt.adapt_structure(to, s.ℋ_of_s),
            Adapt.adapt_structure(to, s.ℋ_p_of_s),
            Adapt.adapt_structure(to, s.ℛ_LD_of_s),
            Adapt.adapt_structure(to, s.ℛ_GNC1_of_s),
            Adapt.adapt_structure(to, s.ℛ_GNC2_of_s), Adapt.adapt_structure(to, s.s_of_z), s.z_eff,
            s.s_min,
            s.s_max,
            s.s_eff,
            s.s_spline_lim,
            s.volume,
        )
    else
        return DevCosmology(
            Adapt.adapt_structure(to, s.IPS; devfloat=devfloat),
            Adapt.adapt_structure(to, s.params; devfloat=devfloat),
            Adapt.adapt_structure(to, s.tools; devfloat=devfloat), Adapt.adapt_structure(to, s.z_of_s; devfloat=devfloat),
            Adapt.adapt_structure(to, s.D_of_s; devfloat=devfloat),
            Adapt.adapt_structure(to, s.f_of_s; devfloat=devfloat),
            Adapt.adapt_structure(to, s.ℋ_of_s; devfloat=devfloat),
            Adapt.adapt_structure(to, s.ℋ_p_of_s; devfloat=devfloat),
            Adapt.adapt_structure(to, s.ℛ_LD_of_s; devfloat=devfloat),
            Adapt.adapt_structure(to, s.ℛ_GNC1_of_s; devfloat=devfloat),
            Adapt.adapt_structure(to, s.ℛ_GNC2_of_s; devfloat=devfloat), Adapt.adapt_structure(to, s.s_of_z; devfloat=devfloat), adapt(to, devfloat(s.z_eff)),
            adapt(to, devfloat(s.s_min)),
            adapt(to, devfloat(s.s_max)),
            adapt(to, devfloat(s.s_eff)),
            adapt(to, devfloat(s.s_spline_lim)),
            adapt(to, devfloat(s.volume)),
        )
    end
end



Adapt.@adapt_structure DevCosmology


struct DevPoint{F}
    z::F    #::Float64
    #conftime::F    #::Float64
    comdist::F  #::Float64
    #angdist::F #::Float64
    #lumdist::F #::Float64
    D::F    #::Float64
    f::F    #::Float64
    ℋ::F    #::Float64
    ℋ_p::F  #::Float64
    ℛ_LD::F #::Float64
    ℛ_GNC1::F   #::Float64
    ℛ_GNC2::F   #::Float64
    a::F    #::Float64
end

#Point(z, comdist, D, f, ℋ, ℛ_LD) = new(z, comdist, D, f, ℋ, ℛ_LD, 1.0/(1.0+z))
function DevPoint(s, cosmo::DevCosmology)
    z = cosmo.z_of_s(s)
    DevPoint(z, s, cosmo.D_of_s(s), cosmo.f_of_s(s), cosmo.ℋ_of_s(s),
        cosmo.ℋ_p_of_s(s), cosmo.ℛ_LD_of_s(s), cosmo.ℛ_GNC1_of_s(s), cosmo.ℛ_GNC2_of_s(s),
        1.0 / (1.0 + z))
end


function Adapt.adapt_structure(to, s::Point; devfloat=DevFloat)
    if devfloat == typeof(s.z)
        return DevPoint(
            s.z, s.comdist, s.D, s.f, s.ℋ, s.ℋ_p, s.ℛ_LD, s.ℛ_GNC1, s.ℛ_GNC2, s.a
        )
    else
        return DevPoint(
            adapt(to, devfloat(s.z)),
            adapt(to, devfloat(s.comdist)),
            adapt(to, devfloat(s.D)),
            adapt(to, devfloat(s.f)),
            adapt(to, devfloat(s.ℋ)),
            adapt(to, devfloat(s.ℋ_p)),
            adapt(to, devfloat(s.ℛ_LD)),
            adapt(to, devfloat(s.ℛ_GNC1)),
            adapt(to, devfloat(s.ℛ_GNC2)),
            adapt(to, devfloat(s.a)),
        )
    end
end

Adapt.@adapt_structure DevPoint

