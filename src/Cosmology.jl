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


@doc raw"""
     func_ℛ(s, ℋ; s_lim=0.01, ℋ_0 = ℋ0)

Return the following value:
```math
\mathrm{func_ℛ}(s, \scrH)=
\begin{cases}
1 - \frac{1}{\scrH \, s} \; ,
    \quad s > s_\mathrm{lim}\\
1 - \frac{1}{\scrH_0 \, s_\mathrm{lim}} \; , 
     \quad \quad 0 \leq s \leq s_\mathrm{lim}
\end{cases}
```
"""
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
          D_of_s::Dierckx.Spline1D
          f_of_s::Dierckx.Spline1D
          ℋ_of_s::Dierckx.Spline1D
          ℛ_of_s::Dierckx.Spline1D

          s_of_z::Dierckx.Spline1D

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

## Arguments 

- `IPS::InputPS` : the matter Input Power Spectrum of the Universe we are focusiong on.

- `params::CosmoParams` : options and parameters decided for this Cosmology.

- `tools::IPSTools` : all the functions and integrals depending on the Input PS.

- `windowF::WindowF` : the window function `F`, defined as:
  ```math
     \begin{split}
     F(x,\mu; \theta_\mathrm{max}) = & \;4\pi 
     \int_0^{\theta_\mathrm{max}} \mathrm{d}\theta_1 \int_0^\pi \mathrm{d} \theta \; 
     \, \Theta\left(\frac
          {x \cos \theta + \cos \theta_1}{\sqrt{x^1+2+2x\mu}} - 
          \cos(\theta_\mathrm{max}) 
          \right) 
     \, \Theta(\mu-\cos(\theta+\theta_1)) \\
     &\Theta(\cos(\theta - \theta_1)-\mu) \;
     \frac{\sin\theta\sin\theta_1}
          {\sqrt{(\sin\theta\sin\theta_1)^2-(\cos\theta\cos\theta_1-\mu)^2}}
     \end{split}
  ```

- `z_of_s, D_of_s, f_of_s, ℋ_of_s, ℛ_of_s ::Dierckx.Spline1D` : splines that returns the
  value of `z`, `D`, `f`, `ℋ` and `ℛ` (respectively) corresponding to an input comoving
  distance `s`. These splines are obtained from the data stored by `BackgroundData` applied
  to the input background data file.

- `s_of_z ::Dierckx.Spline1D` : splines that returns the value of the comoving distance `s`
  corresponding to an input redshift `z`. Also this spline is obtained from the data stored by 
  `BackgroundData` applied to the input background data file.

- `z_eff::Float64` : effective redshift of this survey; its value is obtained through
  the function `func_z_eff`, with inputs the `s_min`, `s_max` and `z_of_s` here stored.

- `s_min::Float64` and `s_max::Float64` : the minimum and maximum comoving distances of
  the survey considered; they are the corresponding comoving distance to the chosen minimum and
  maximum redshifts `z_min` and `z_max`, stored in the input `CosmoParams`.

- `s_eff::Float64` : the corresponding comoving distance to the computed effective 
  redshifts `z_eff`.

- `volume::Float64` : volume of this survey. It is computed applying the function `V_survey`
  with inputs `s_min`, `s_max` here stored and the `θ_max` in the input `CosmoParams`.

- `file_data, file_ips, file_windowF::String` : the file names used for this Cosmology.

## Constructors

`Cosmology(
     params::CosmoParams,
     file_data::String,
     file_ips::String,
     file_windowF::String,
     file_Is::Union{String,Nothing} = nothing;
     names_bg = NAMES_BACKGROUND)`

- `params::CosmoParams` : parameters to be used for this Cosmology. See the docstring
  of `CosmoParams` for more information on the possible inputs.

- `file_data::String` : file containing all the background data; it is expected that such file
  is a background output of the CLASS program (link: https://github.com/lesgourg/class_public).
  It is managed through `BackgroundData`.

- `file_ips::String` : file containing the Input Power Spectrum; it is expected that such file
  is a power spectrum output of the CLASS program (link: https://github.com/lesgourg/class_public).
  It is managed through `InputPS`.

- `file_windowF::String` : file containing a map of the window function `F`.
  This file is managed through `WindowF`, and can be produced with `F_map`; see their
  docstrings for more information.

- `file_Is::Union{String,Nothing} = nothing` : if you want to given in input manually
  all the ``I_\ell^n`` integrals, you can set as input the file containing them.
  It is expected that they are ordered in colums with the following order:
  `s  I00  I20  I40  I02  I22  I31  I11  I13  I04_tilde`.
  If nothing is passed (recommended), they are manually calculated from the Input Power Spectrum.

- `names = NAMES_BACKGROUND` : the column names of the `file_data`. If the colum order change from
  the default one `NAMES_BACKGROUND`, you must set as input the vector of string with the correct
  one, with the SAME names. They are, with the default order:\n
  $(NAMES_BACKGROUND)  

See also:  [`InputPS`](@ref), [`CosmoParams`](@ref), [`IPSTools`](@ref),
[`BackgroundData`](@ref), [`WindowF`](@ref), [`F_map`](@ref), [`func_z_eff`](@ref),
[`V_survey`](@ref)
"""
struct Cosmology
     IPS::InputPS
     params::CosmoParams
     tools::IPSTools
     windowF::WindowF
     
     z_of_s::Dierckx.Spline1D
     D_of_s::Dierckx.Spline1D
     f_of_s::Dierckx.Spline1D
     ℋ_of_s::Dierckx.Spline1D
     ℛ_of_s::Dierckx.Spline1D

     s_of_z::Dierckx.Spline1D

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
               IPSTools(IPS; k_min = params.IPS[:k_min], k_max = params.IPS[:k_max],
               N = params.IPS[:N], fit_min = params.IPS[:fit_min],
               fit_max = params.fit_max, con = params.con) :
               IPSTools(IPS, file_Is)

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
          ℛs = [func_ℛ(s, ℋ_of_s(s); s_lim =  params.s_lim) for s in ss]
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
               z_of_s, D_of_s, f_of_s, ℋ_of_s, ℛ_of_s,
               s_of_z,
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
