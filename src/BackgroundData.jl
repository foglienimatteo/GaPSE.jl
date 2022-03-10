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

"""
     const f0 :: Float64

Linear growth rate at present time. Its value is equal to:
```math
     f_0 \\simeq 0.5126998572951
```
"""
const f0 = 5.126998572951e-01


"""
     const D0 :: Float64

Linear growth factor at present time. Its value is equal to:
```math
     D_0 = 1.0
```
"""
const D0 = 1.0


"""
     const ℋ0 :: Float64

Comoving Hubble constant at present time. Its value is, in natural system
(where the speed of light c=1): 
``\\mathcal{H}_0 \\simeq 3.335641\\times10^{-4} \\; h_0^{-1}\\mathrm{Mpc}``
"""
const ℋ0 = 3.3356409519815204e-4 # h_0/Mpc



##########################################################################################92


"""
     BackgroundData(
          z::Vector{Float64}
          conftime::Vector{Float64}
          comdist::Vector{Float64}
          angdist::Vector{Float64}
          lumdist::Vector{Float64}
          D::Vector{Float64}
          f::Vector{Float64}
          ℋ::Vector{Float64}
          ℋ_p::Vector{Float64})

Struct that contains all the relevant cosmological information for
future computations.
The data are stored with increasing distance values 
(so the first ones are associated to `z=0`).
It is internally used in `Cosmology.`

## Arguments

- `z::Vector{Float64}` : redshifts (adimensionals).

- `conftime::Vector{Float64}` : conformal times, measured in [Mpc/h].

- `comdist::Vector{Float64}` : comoving distances, measured in [Mpc/h].

- `angdist::Vector{Float64}` : angular diameter distances, measured in [Mpc/h].

- `lumdist::Vector{Float64}` : luminosity distances, measured in [Mpc/h].

- `D::Vector{Float64}` : linear growth factors, normalized to 1.0 at the present day (adimensional).

- `f::Vector{Float64}` : linear growth rates (adimensional).

- `ℋ::Vector{Float64}` : comoving Hubble parameters, measured in [h/Mpc].

- `ℋ_p::Vector{Float64}` : derivatives of the comoving Hubble parameter wrt the conformal time.
  It is here manually computed with the Dierckx function `derivative`.


## Constructors

`BackgroundData(file::String, z_max; names = NAMES_BACKGROUND, h = 0.7)`

- `file::string` : input file where the data are stored; it is expected that such file
  is a background output of the CLASS program (link: https://github.com/lesgourg/class_public)

- `z_max` : the maximum redhsift we are interested in our analysis. The constructor will
  store the data necessary for a study only in `0 < z < z_max`, for optimisation purposes
  (More precisely, the maximum distance stored will be `3*z_max`).

- `names = NAMES_BACKGROUND` : the column names of the `file`. If the colum order change from
  the default one `NAMES_BACKGROUND`, you must set as input the vector of string with the correct
  one, with the SAME names. They are, with the default order:\n
  $(NAMES_BACKGROUND)

- `h = 0.7` : the adimensional hubble constant. By default, CLASS background data are measured with
  it nuymerically expressed (so distances are measured in `Mpc`, for example), while this code works
  with `h` in the unit of measure (so distances are measured in `Mpc/h`, for example).
  Change this value to `1.0` if the input data do not have this issue, or to your value of interest 
  (`0.67`, `0.5`, ...).

See also: [`CosmoParams`](@ref), [`Cosmology`](@ref)
"""
struct BackgroundData
     z::Vector{Float64}
     conftime::Vector{Float64}
     comdist::Vector{Float64}
     angdist::Vector{Float64}
     lumdist::Vector{Float64}
     D::Vector{Float64}
     f::Vector{Float64}
     ℋ::Vector{Float64}
     ℋ_p::Vector{Float64}

     function BackgroundData(file::String, z_max;
          names = NAMES_BACKGROUND, h = 0.7)
     
          I_redshift = findfirst(x -> x == "z", names)
          I_comdist = findfirst(x -> x == "comov. dist.", names)
     
          data = readdlm(file, comments = true)
     
          N_z_MAX = findfirst(z -> z <= z_max, data[:, I_redshift]) - 1
          com_dist_z_MAX = data[:, I_comdist][N_z_MAX]
          N_2_com_dist_z_MAX = findfirst(s -> s <= 3.0 * com_dist_z_MAX, data[:, I_comdist]) - 1
     
          data_dict = Dict([name => reverse(data[:, i][N_2_com_dist_z_MAX:end])
                            for (i, name) in enumerate(names)]...)
     
          com_H = data_dict["H [1/Mpc]"] ./ h ./ (1.0 .+ data_dict["z"])
          conf_time = data_dict["conf. time [Mpc]"] .* h
          spline_com_H = Spline1D(reverse(conf_time), reverse(com_H); bc = "nearest")
          com_H_p = [Dierckx.derivative(spline_com_H, t) for t in conf_time]
     
          new(
               data_dict["z"],
               conf_time,
               data_dict["comov. dist."] .* h,
               data_dict["ang.diam.dist."] .* h,
               data_dict["lum. dist."] .* h,
               data_dict["gr.fac. D"],
               data_dict["gr.fac. f"],
               com_H,
               com_H_p,
          )
     end
end



DEFAULT_IPS_OPTS = Dict(
     :fit_left_min => 1e-6::Float64, 
     :fit_left_max => 3e-6::Float64,
     :fit_right_min => 1e1::Float64, 
     :fit_right_max => 2e1::Float64,
     )

DEFAULT_IPSTOOLS_OPTS = Dict(
     :N => 1024::Integer, 
     :fit_min => 0.05::Float64, 
     :fit_max => 0.5::Float64, 
     :con => true::Bool,
     :k_min => 1e-6::Float64,
     :k_max => 10.0::Float64,
     )

"""
     CosmoParams(
          z_min::Float64
          z_max::Float64
          θ_max::Float64

          k_min::Float64
          k_max::Float64

          Ω_b::Float64
          Ω_cdm::Float64
          Ω_M0::Float64
          h_0::Float64

          N::Integer
          fit_min::Float64
          fit_max::Float64
          con::Bool
          s_lim::Float64)


Struct that contains all the parameters and options that are 
matter of concerns for the `Cosmology` we are interested in.

## Arguments

- `z_min::Float64` and `z_max::Float64` : the minimum and maximum redshifts of the
  survey we want to study.

- `θ_max::Float64` : Angular maximum value of the survey. It is
  implicitly assumed an azimutal simmetry of the survey.

- `k_min::Float64` and `k_max::Float64` : extremes of integration for the `σ_i`
  integrals in `IPSTools`.

- `Ω_b::Float64`, `Ω_cdm::Float64` and `Ω_M0::Float64` : barionic, cold-dark-matter and
  total matter density parameters.

- `h_0::Float64` : today's Hubble adimensional parameter (`H_0 = h_0 * 100 km/(s * Mpc)`).

- `N::Integer` : number of points to be used in the Sperical Bessel Fourier Transform made
  by `xicalc` in `IPSTools`.

- `fit_min::Float64` and `fit_max::Float64` : the limits (min and max) where the integral ``I_\\ell^n``
  in `Cosmology` must be fitted with a power law, for small distances. This operation is necessary, 
  because `xicalc`, in this context, gives wrong results for too small input distance `s`; nevertheless, all
  these ``I_\\ell^n`` integrals have fixed power-law trends for ``s \\rightarrow 0``, so this approach gives
  good results.

- `con::Bool` : do you want that the fit of all the ``I_\\ell^n`` in `IPSTools` for the LEFT edge
  is not a simple power-law ``y = f(x) = b \\, x^s``, but also consider a constant ``a``,
  such that ``y = f(x) = a + b \\, x^s``?

- `s_lim::Float64` : the lower-bound value for the function `func_ℛ`; it is necessary, because
  `ℛ` blows up for ``s \\rightarrow 0^{+}``. Consequently, if the `func_ℛ` input value is 
  `0 ≤ s < s_lim`, the returned value is always `func_ℛ(s_lim)`.

## Constructors

`CosmoParams(z_min, z_max, θ_max; k_min = 1e-8, k_max = 10.0,
     Ω_b = 0.0489, Ω_cdm = 0.251020, h_0 = 0.70, N::Integer = 1024, 
     fit_min = 0.05, fit_max = 0.5, con::Bool = true, 
     s_lim = 1e-2)`
The associations are trivials, with `Ω_M0 = Ω_cdm + Ω_b`.

See also: [`Cosmology`](@ref), [`IPSTools`](@ref), [`func_ℛ`](@ref)
"""
struct CosmoParams
     z_min::Float64
     z_max::Float64
     θ_max::Float64

     Ω_b::Float64
     Ω_cdm::Float64
     Ω_M0::Float64
     h_0::Float64

     s_lim::Float64

     IPS::Dict{Symbol, T} where T
     IPSTools::Dict{Symbol, T} where T

     function CosmoParams(z_min, z_max, θ_max;
               Ω_b = 0.0489, Ω_cdm = 0.251020, h_0 = 0.70, s_lim = 1e-2,
               IPS_opts::Dict = Dict{Symbol, Any}(),
               IPSTools_opts::Dict = Dict{Symbol, Any}(),
               ) 

          @assert typeof(IPS_opts) <: Dict{Symbol, T} where T "the keys of "*
          "the IPS_opts dict have to be Symbols (like :k_min, :N, ...)"

          @assert typeof(IPSTools_opts) <: Dict{Symbol, T} where T "the keys of "*
          "the IPSTools_opts dict have to be Symbols (like :k_min, :N, ...)"

          check_compatible_dicts(DEFAULT_IPS_OPTS, IPS_opts, "IPS_opts")
          check_compatible_dicts(DEFAULT_IPSTOOLS_OPTS, IPSTools_opts, "IPSTools_opts")

          IPS = merge(DEFAULT_IPS_OPTS, IPS_opts)
          IPSTools = merge(DEFAULT_IPSTOOLS_OPTS, IPSTools_opts)
     
          @assert 0.0 < z_min < z_max " 0.0 < z_min < z_max must hold!"
          @assert 0.0 ≤ θ_max ≤ π / 2.0 " 0.0 ≤ θ_max ≤ π/2.0 must hold!"
          @assert 0.0 ≤ Ω_b ≤ 1.0 " 0.0 ≤ Ω_b ≤ 1.0 must hold!"
          @assert 0.0 ≤ Ω_cdm ≤ 1.0 " 0.0 ≤ Ω_cdm ≤ 1.0 must hold!"
          @assert 0.0 < h_0 ≤ 1.0 " 0.0 < h_0 ≤ 1.0 must hold!"
          @assert 0.0 < s_lim < 10.0 "0.0 < s_lim < 10.0 must hold!"

          @assert 0.0 < IPS[:fit_left_min] < IPS[:fit_left_max] < 1e-1 
               " 0 < fit_left_min < fit_left_max < 0.1 must hold!"
          @assert 0.5 < IPS[:fit_right_min] < IPS[:fit_right_max] < 1e6
               " 0.5 < fit_right_min < fit_right_max < 1e6 must hold!"

          @assert 0.0 ≤ IPSTools[:k_min] < IPSTools[:k_max] " 0.0 ≤ k_min < k_max must hold!"
          @assert IPSTools[:N] > 7 " N > 7 must hold!"
          @assert 1e-2 ≤ IPSTools[:fit_min] < IPSTools[:fit_max] < 10.0 " 1e-2 "*
               "≤ fit_min < fit_max < 10.0 must hold!"
     
          new(z_min, z_max, θ_max, Ω_b, Ω_cdm, Ω_cdm + Ω_b, h_0, s_lim,
               IPS, IPSTools)
     end
end



##########################################################################################92



"""
     func_z_eff(s_min, s_max, z_of_s) :: Float64

Return the effective redshift ``z_\\mathrm{eff}``, calcuated as follows:
```math
\\begin{split}
z_\\mathrm{eff} := 
    \\frac{
        \\int \\mathrm{d}^3\\mathbf{s} \\, \\phi^2(\\mathbf{s}) \\, z(s)
     }{
         \\int \\mathrm{d}^3\\mathbf{s}\\, \\phi^2(\\mathbf{s}) 
      } &= \\frac{
          \\int_0^\\infty \\mathrm{d}s  \\, s^2 \\, \\phi^2(s) \\, z(s) \\times
          \\int_{4\\pi}\\mathrm{d}^2\\hat{\\mathbf{s}} \\, W^2(\\hat{\\mathbf{s}})
      }{
          \\int_0^\\infty \\mathrm{d}s \\, s^2 \\, \\phi^2(s)\\times
          \\int_{4\\pi}\\mathrm{d}^2\\hat{\\mathbf{s}} \\, W^2(\\hat{\\mathbf{s}})
      } \\\\[5pt]
      &= \\frac{
          \\int_0^\\infty \\mathrm{d}s  \\, s^2 \\, \\phi^2(s) \\, z(s)
      }{
          \\int_0^\\infty \\mathrm{d}s \\, s^2 \\, \\phi^2(s)
      } \\\\[4pt]
      &= \\frac{3}{s_\\mathrm{max}^3 - s_\\mathrm{min}^3} \\,
          \\int_{s_\\mathrm{min}}^{s_\\mathrm{max}} \\mathrm{d}s  \\, s^2 \\, z(s)
\\end{split}
```
where we have used our assuption on separability of the window function
```math
     \\phi(\\mathbf{s}) = \\phi(s) \\, W(\\hat{s})
```
and their definitions.


See also: [`ϕ`](@ref), [`W`](@ref)
"""
function func_z_eff(s_min, s_max, z_of_s)
     3.0 / (s_max^3 - s_min^3) * quadgk(s -> s^2 * z_of_s(s), s_min, s_max)[1]
end


"""
     s(s1, s2, y) :: Float64

Return the value ``s = \\sqrt{s_1^2 + s_2^2 - 2 \\, s_1 \\, s_2 \\, y}``

See also: [`μ`](@ref), [`s2`](@ref), [`y`](@ref)
"""
s(s1, s2, y) = √(s1^2 + s2^2 - 2 * s1 * s2 * y)


"""
     μ(s1, s2, y) :: Float64

Return the value ``\\mu=\\hat{\\mathbf{s}}_1\\cdot\\hat{\\mathbf{s}}``, defined as:
```math
\\mu = \\mu(s_1, s_2, y) = \\frac{y \\, s_2 - s_1}{s(s_1, s_2, y)} \\;,
\\quad s(s_1, s_2, y) = \\sqrt{s_1^2 + s^2 - 2 \\, s_1 \\, s_2 \\, y}
```
with ``y=\\cos\\theta=\\hat{\\mathbf{s}}_1\\cdot\\hat{\\mathbf{s}}`` and where ``s`` is 
obtained from the function `s`

See also: [`s`](@ref), [`s2`](@ref), [`y`](@ref)
"""
μ(s1, s2, y) = (y * s2 - s1) / s(s1, s2, y)



"""
     s2(s1, s, μ) :: Float64

Return the value ``s_2 = \\sqrt{s_1^2 + s^2 + 2 \\, s_1 \\, s \\, \\mu}``

See also: [`s`](@ref), [`μ`](@ref), [`y`](@ref)
"""
s2(s1, s, μ) = √(s1^2 + s^2 + 2 * s1 * s * μ)



"""
     y(s1, s, μ) :: Float64

Return the value ``y=\\cos\\theta``, defined as:
```math
y = y(s_1, s, \\mu) = \\frac{\\mu \\, s + s_1}{s2(s_1, s, \\mu)} \\;,
\\quad s_2 = \\sqrt{s_1^2 + s^2 + 2 \\, s_1 \\, s \\, \\mu}
```
with ``\\mu=\\hat{\\mathbf{s}}_1\\cdot\\hat{\\mathbf{s}}_2`` and 
where ``s_2`` is btained from the function `s2`

See also: [`s`](@ref), [`μ`](@ref), [`s2`](@ref)
"""
y(s1, s, μ) = (μ * s + s1) / s2(s, s1, μ)




