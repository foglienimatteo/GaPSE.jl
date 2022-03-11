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
     const DEFAULT_IPS_OPTS = Dict(
          :fit_left_min => 1e-6::Float64, 
          :fit_left_max => 3e-6::Float64,
          :fit_right_min => 1e1::Float64, 
          :fit_right_max => 2e1::Float64,
          )

The default values to be stored in `CosmoParams` concerning the 
Input Power Spectrum. In the `Cosmology` that will have such `CosmoParams` as input,
they will be used in its `InputPS`.

See also: [`CosmoParams`](@ref), [`Cosmology`](@ref), [`InputPS`](@ref)
"""
const DEFAULT_IPS_OPTS = Dict(
     :fit_left_min => 1e-6::Float64,
     :fit_left_max => 3e-6::Float64,
     :fit_right_min => 1e1::Float64,
     :fit_right_max => 2e1::Float64,
)


"""
     const DEFAULT_IPSTOOLS_OPTS = Dict(
          :N => 1024::Integer,
          :fit_min => 0.05::Float64,
          :fit_max => 0.5::Float64,
          :con => true::Bool,
          :k_min => 1e-6::Float64,
          :k_max => 10.0::Float64,
     )

The default values to be stored in `CosmoParams` concerning the 
Input Power Spectrum Tools. In the `Cosmology` that will have such `CosmoParams` as input,
they will be used in its `IPSTools`.

See also: [`CosmoParams`](@ref), [`Cosmology`](@ref), [`IPSTools`](@ref)
"""
const DEFAULT_IPSTOOLS_OPTS = Dict(
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

     IPS::Dict{Symbol,T} where {T}
     IPSTools::Dict{Symbol,T} where {T}

     function CosmoParams(z_min, z_max, θ_max;
          Ω_b = 0.0489, Ω_cdm = 0.251020, h_0 = 0.70, s_lim = 1e-2,
          IPS_opts::Dict = Dict{Symbol,Any}(),
          IPSTools_opts::Dict = Dict{Symbol,Any}()
     )

          @assert typeof(IPS_opts) <: Dict{Symbol,T} where {T} "the keys of " *
                                                               "the IPS_opts dict have to be Symbols (like :k_min, :N, ...)"

          @assert typeof(IPSTools_opts) <: Dict{Symbol,T} where {T} "the keys of " *
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
          @assert 1e-2 ≤ IPSTools[:fit_min] < IPSTools[:fit_max] < 10.0 " 1e-2 " *
                                                                        "≤ fit_min < fit_max < 10.0 must hold!"

          new(z_min, z_max, θ_max, Ω_b, Ω_cdm, Ω_cdm + Ω_b, h_0, s_lim,
               IPS, IPSTools)
     end
end

