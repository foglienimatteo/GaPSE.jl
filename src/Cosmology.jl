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


function func_ℛ(s, ℋ; s_lim=0.01, ℋ_0 = ℋ0)
     if s > s_lim
          return 1.0 - 1.0/(s*ℋ)
     else
          return 1.0 - 1.0/(s_lim*ℋ_0)
     end
end


@doc raw"""
     Cosmology(
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
          )

Struct that contains all the information that may be used for the 
Correlation Function computations.


See also:  [`InputPS`](@ref), [`CosmoParams`](@ref), [`IPSTools`](@ref),
[`WindowF`](@ref)
"""
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
               file_Is::Union{String,Nothing} = nothing;
               names_bg = NAMES_BACKGROUND
          )

          BD = BackgroundData(file_data::String, params.z_max;
               names = names_bg, h = params.h_0)
          IPS = InputPS(file_ips)
          windowF = WindowF(file_windowF)
          tools = isnothing(file_Is) ?
               IPSTools(IPS; k_min = params.k_min, k_max = params.k_max,
               N = params.N, fit_min = params.fit_min,
               fit_max = params.fit_max, con = params.con) :
               IPSTools(IPS, file_Is)

          s_lim = isnothing(params.s_lim) ? 1.0 : params.s_lim

          #=
          z_of_s_lim = my_interpolation(BD.comdist[1], BD.z[1], BD.comdist[2], BD.z[2], s_lim)
          D_of_s_lim = my_interpolation(BD.comdist[1], BD.D[1], BD.comdist[2], BD.D[2], s_lim)
          f_of_s_lim = my_interpolation(BD.comdist[1], BD.f[1], BD.comdist[2], BD.f[2], s_lim)
          ℋ_of_s_lim = my_interpolation(BD.comdist[1], BD.ℋ[1], BD.comdist[2], BD.ℋ[2], s_lim)

          new_BD_comdist = vcat(0.0, s_lim, BD.comdist[2:end])
          new_BD_z = vcat(0.0, z_of_s_lim, BD.z[2:end])
          new_BD_D = vcat(D_of_s_lim, D_of_s_lim, BD.D[2:end])
          new_BD_f = vcat(f_of_s_lim, f_of_s_lim, BD.f[2:end])
          new_BD_ℋ = vcat(ℋ_of_s_lim, ℋ_of_s_lim, BD.ℋ[2:end])

          another_BD_comdist = vcat(s_lim, s_lim, BD.comdist[2:end])
          another_BD_z = vcat(z_of_s_lim, z_of_s_lim, BD.z[2:end])

          z_of_s = Spline1D(new_BD_comdist, another_BD_z; bc = "error")
          s_of_z = Spline1D(new_BD_z, another_BD_comdist; bc = "error")
          D_of_s = Spline1D(new_BD_comdist, new_BD_D; bc = "error")
          f_of_s = Spline1D(new_BD_comdist, new_BD_f; bc = "error")
          ℋ_of_s = Spline1D(new_BD_comdist, new_BD_ℋ; bc = "error")
          =#

          z_of_s = Spline1D(BD.comdist, BD.z; bc = "error")
          s_of_z = Spline1D(BD.z, BD.comdist; bc = "error")
          D_of_s = Spline1D(BD.comdist, BD.D; bc = "error")
          f_of_s = Spline1D(BD.comdist, BD.f; bc = "error")
          ℋ_of_s = Spline1D(BD.comdist, BD.ℋ; bc = "error")

          ss = 10 .^ range(-4, log10(BD.comdist[end]), length = 1000)
          ℛs = [func_ℛ(s, ℋ_of_s(s); s_lim = s_lim) for s in ss]
          ℛ_of_s = Spline1D(vcat(0.0, ss), vcat(ℛs[begin], ℛs); bc = "error")

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
               z_of_s, s_of_z,
               D_of_s, f_of_s, ℋ_of_s, ℛ_of_s,
               z_eff, s_min, s_max, s_eff,
               vol,
               file_data,
               file_ips,
               file_windowF
          )
     end
end



"""
     Point(
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
          a::Float64)
     
A point in the Universe, placed at redshift `z` from us.
It contains all the relevant cosmological information at that redshift.

## Constructors

`Point(s, cosmo::Cosmology)` : given a comoving distance `s`, extrapolate all 
the data from the given input `Cosmology`.

See also: [`Cosmology`](@ref)
"""
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

     #Point(z, comdist, D, f, ℋ, ℛ) = new(z, comdist, D, f, ℋ, ℛ, 1.0/(1.0+z))
     function Point(s, cosmo::Cosmology)
          z = cosmo.z_of_s(s)
          new(z, s, cosmo.D_of_s(s), cosmo.f_of_s(s),
               cosmo.ℋ_of_s(s), cosmo.ℛ_of_s(s), 1.0/(1.0+z))
     end
end
