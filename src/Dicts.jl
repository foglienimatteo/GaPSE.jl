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
     const IMPLEMENTED_GR_EFFECTS = [
          "auto_doppler", "auto_lensing",
          "auto_localgp", "auto_integratedgp", 
          
          "lensing_doppler", "doppler_lensing",
          "doppler_localgp", "localgp_doppler",
          "doppler_integratedgp", "integratedgp_doppler",
          "lensing_localgp", "localgp_lensing",
          "lensing_integratedgp", "integratedgp_lensing",
          "localgp_integratedgp", "integratedgp_localgp",
     ]

The names of the GR effects implemented.
Their order is associated with the one in `IMPLEMENTED_ξs`, so be careful
to change it.

See also: [`IMPLEMENTED_ξs`](@ref)
"""
const IMPLEMENTED_GR_EFFECTS = [
     "auto_doppler", "auto_lensing",
     "auto_localgp", "auto_integratedgp", 
     
     "lensing_doppler", "doppler_lensing",
     "doppler_localgp", "localgp_doppler",
     "doppler_integratedgp", "integratedgp_doppler",
     "lensing_localgp", "localgp_lensing",
     "lensing_integratedgp", "integratedgp_lensing",
     "localgp_integratedgp", "integratedgp_localgp",
];


"""
     const IMPLEMENTED_ξs = [
          ξ_Doppler, ξ_Lensing, ξ_LocalGP, ξ_IntegratedGP, 
          ξ_Lensing_Doppler, ξ_Doppler_Lensing,
          ξ_Doppler_LocalGP, ξ_LocalGP_Doppler,
          ξ_Doppler_IntegratedGP, ξ_IntegratedGP_Doppler,
          ξ_Lensing_LocalGP, ξ_LocalGP_Lensing,
          ξ_Lensing_IntegratedGP, ξ_IntegratedGP_Lensing,
          ξ_LocalGP_IntegratedGP, ξ_IntegratedGP_LocalGP
     ]

The names of the GR effect TPCFs implemented.
Their order is associated with the one in `IMPLEMENTED_GR_EFFECTS`, so be careful
to change it.

See also: [`IMPLEMENTED_GR_EFFECTS`](@ref)
"""
const IMPLEMENTED_ξs = [
     ξ_Doppler, ξ_Lensing, ξ_LocalGP, ξ_IntegratedGP, 
     ξ_Lensing_Doppler, ξ_Doppler_Lensing,
     ξ_Doppler_LocalGP, ξ_LocalGP_Doppler,
     ξ_Doppler_IntegratedGP, ξ_IntegratedGP_Doppler,
     ξ_Lensing_LocalGP, ξ_LocalGP_Lensing,
     ξ_Lensing_IntegratedGP, ξ_IntegratedGP_Lensing,
     ξ_LocalGP_IntegratedGP, ξ_IntegratedGP_LocalGP
];


"""
     const DICT_GR_ξs::Dict{String,Function}

For an input key string `effect` from `IMPLEMENTED_GR_EFFECTS`, return the 
associated TPCF `DICT_GR_ξs[effect]` from `IMPLEMENTED_ξs`.

## Example

```jldoctest
julia> GaPSE.DICT_GR_ξs["auto_doppler"]
ξ_Doppler
```

See also: [`IMPLEMENTED_GR_EFFECTS`](@ref), [`IMPLEMENTED_ξs`](@ref)
"""
const DICT_GR_ξs = Dict([k => v for (k, v) in zip(IMPLEMENTED_GR_EFFECTS, IMPLEMENTED_ξs)]...)


"""
     const INDEX_GR_EFFECT::Dict{String,Integer}

For an input key string `effect` from `IMPLEMENTED_GR_EFFECTS`, return the 
associated index position in that vector.

## Example

```jldoctest
julia> GaPSE.INDEX_GR_EFFECT["auto_doppler"]
1
```

See also: [`IMPLEMENTED_GR_EFFECTS`](@ref)
"""
const INDEX_GR_EFFECT = Dict([name => i for (i, name) in
                        enumerate(IMPLEMENTED_GR_EFFECTS)]...)


"""
     const GR_EFFECT_INDEXED::Dict{Integer,String}

For an input index position `i` of `IMPLEMENTED_GR_EFFECTS`, return the 
associated key string `effect`.

## Example

```jldoctest
julia> GaPSE.DICT_GR_ξs[1]
"auto_doppler"
```

See also: [`IMPLEMENTED_GR_EFFECTS`](@ref)
"""
const GR_EFFECT_INDEXED = Dict([i => name for (i, name) in
                          enumerate(IMPLEMENTED_GR_EFFECTS)]...)

