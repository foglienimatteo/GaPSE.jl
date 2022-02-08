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


struct Cosmology
     IPS::InputPS
     params::CosmoParams
     tools::IPSTools
     windowF::WindowF
     z_of_s::Dierckx.Spline1D
     s_of_z::Dierckx.Spline1D
     D_of_s::Dierckx.Spline1D
     f_of_s::Dierckx.Spline1D
     ℋ_of_s::Dierckx.Spline1D
     ℛ_of_s::Dierckx.Spline1D

     z_eff::Float64
     s_min::Float64
     s_max::Float64
     s_eff::Float64

     volume::Float64

     file_data::String
     file_ips::String
     file_windowF::String


     function Cosmology(
          params::CosmoParams,
          file_data::String,
          file_ips::String,
          file_windowF::String,
          file_Is::Union{String, Nothing} = nothing;
          expand::Bool=true,
          names_bg = NAMES_BACKGROUND
     )
     
          BD = BackgroundData(file_data::String, params.z_min, params.z_max;
               names = names_bg, h = params.h_0)
          IPS = InputPS(file_ips; expand = expand)
          windowF = WindowF(file_windowF)
          tools = isnothing(file_Is) ?
                  IPSTools(IPS; k_min = params.k_min, k_max = params.k_max,
                    N = params.N, fit_min = params.fit_min,
                    fit_max = params.fit_max, con = params.con) :
                  IPSTools(IPS, file_Is) 
     
          ℛs = 1.0 .- 1.0 ./ (BD.ℋ .* BD.comdist)
     
          z_of_s = Spline1D(BD.comdist, BD.z)
          s_of_z = Spline1D(BD.z, BD.comdist)
     
          s_min = s_of_z(params.z_min)
          s_max = s_of_z(params.z_max)
          z_eff = func_z_eff(s_min, s_max, z_of_s)
          s_eff = s_of_z(z_eff)
          vol = V_survey(s_min, s_max, params.θ_max)
          new(
               IPS,
               params,
               tools,
               windowF,
               z_of_s,
               s_of_z,
               Spline1D(BD.comdist, BD.D),
               Spline1D(BD.comdist, BD.f),
               Spline1D(BD.comdist, BD.ℋ),
               Spline1D(BD.comdist, ℛs),
               z_eff, s_min, s_max, s_eff,
               vol,
               file_data,
               file_ips,
               file_windowF
          )
     end
end




struct Point
     z::Float64
     #conftime::Float64
     comdist::Float64
     #angdist::Float64
     #lumdist::Float64
     D::Float64
     f::Float64
     ℋ::Float64
     #ℋ_p::Float64
     ℛ::Float64
     a::Float64

     Point(z, comdist, D, f, ℋ, ℛ) = new(z, comdist, D, f, ℋ, ℛ, 1.0/(1.0+z))
     function Point(s, cosmo::Cosmology)
          z = cosmo.z_of_s(s)
          new(z, s, cosmo.D_of_s(s), cosmo.f_of_s(s),
               cosmo.ℋ_of_s(s), cosmo.ℛ_of_s(s), 1.0/(1.0+z))
     end
end
