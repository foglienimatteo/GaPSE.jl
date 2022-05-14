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


##################### PERTURBED LUMINOSITY DISTANCE #############################



"""
     const GR_EFFECTS_LD = [
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
Their order is associated with the one in `VEC_ξs_LD`, so be careful
to change it.

See also: [`VEC_ξs_LD`](@ref)
"""
const GR_EFFECTS_LD = [
     "auto_doppler", "auto_lensing",
     "auto_localgp", "auto_integratedgp", "lensing_doppler", "doppler_lensing",
     "doppler_localgp", "localgp_doppler",
     "doppler_integratedgp", "integratedgp_doppler",
     "lensing_localgp", "localgp_lensing",
     "lensing_integratedgp", "integratedgp_lensing",
     "localgp_integratedgp", "integratedgp_localgp",
];


"""
     const VEC_ξs_LD = [
          ξ_LD_Doppler, ξ_LD_Lensing, ξ_LD_LocalGP, ξ_LD_IntegratedGP, 
          ξ_LD_Lensing_Doppler, ξ_LD_Doppler_Lensing,
          ξ_LD_Doppler_LocalGP, ξ_LD_LocalGP_Doppler,
          ξ_LD_Doppler_IntegratedGP, ξ_LD_IntegratedGP_Doppler,
          ξ_LD_Lensing_LocalGP, ξ_LD_LocalGP_Lensing,
          ξ_LD_Lensing_IntegratedGP, ξ_LD_IntegratedGP_Lensing,
          ξ_LD_LocalGP_IntegratedGP, ξ_LD_IntegratedGP_LocalGP
     ]

The names of the GR effect TPCFs implemented.
Their order is associated with the one in `GR_EFFECTS_LD`, so be careful
to change it.

See also: [`GR_EFFECTS_LD`](@ref)
"""
const VEC_ξs_LD = [
     ξ_LD_Doppler, ξ_LD_Lensing, ξ_LD_LocalGP, ξ_LD_IntegratedGP,
     ξ_LD_Lensing_Doppler, ξ_LD_Doppler_Lensing,
     ξ_LD_Doppler_LocalGP, ξ_LD_LocalGP_Doppler,
     ξ_LD_Doppler_IntegratedGP, ξ_LD_IntegratedGP_Doppler,
     ξ_LD_Lensing_LocalGP, ξ_LD_LocalGP_Lensing,
     ξ_LD_Lensing_IntegratedGP, ξ_LD_IntegratedGP_Lensing,
     ξ_LD_LocalGP_IntegratedGP, ξ_LD_IntegratedGP_LocalGP
];



"""
     const DICT_GR_ξs_LD::Dict{String,Function}

For an input key string `effect` from `GR_EFFECTS_LD`, return the 
associated TPCF `DICT_GR_ξs_LD[effect]` from `VEC_ξs_LD`.

## Example

```jldoctest
julia> GaPSE.DICT_GR_ξs_LD["auto_doppler"]
ξ_LD_Doppler
```

See also: [`GR_EFFECTS_LD`](@ref), [`VEC_ξs_LD`](@ref)
"""
const DICT_GR_ξs_LD = Dict([k => v for (k, v) in zip(GR_EFFECTS_LD, VEC_ξs_LD)]...)


"""
     const INDEX_GR_EFFECT_LD::Dict{String,Integer}

For an input key string `effect` from `GR_EFFECTS_LD`, return the 
associated index position in that vector.

## Example

```jldoctest
julia> GaPSE.INDEX_GR_EFFECT_LD["auto_doppler"]
1
```

See also: [`GR_EFFECTS_LD`](@ref)
"""
const INDEX_GR_EFFECT_LD = Dict([name => i for (i, name) in
                                 enumerate(GR_EFFECTS_LD)]...)


"""
     const GR_EFFECT_INDEX_LD::Dict{Integer,String}

For an input index position `i` of `GR_EFFECTS_LD`, return the 
associated key string `effect`.

## Example

```jldoctest
julia> GaPSE.DICT_GR_ξs_LD[1]
"auto_doppler"
```

See also: [`GR_EFFECTS_LD`](@ref)
"""
const GR_EFFECT_INDEX_LD = Dict([i => name for (i, name) in
                                 enumerate(GR_EFFECTS_LD)]...)



##################### GALAXY NUMBER COUNTS #############################




"""
     const GR_EFFECTS_GNC = [
          "auto_newton", "auto_doppler", "auto_lensing",
          "auto_localgp", "auto_integratedgp", 
          
          "newton_doppler", "doppler_newton",
          "newton_lensing", "lensing_newton",
          "newton_localgp", "localgp_newton",
          "newton_integratedgp", "integratedgp_newton",
          "lensing_doppler", "doppler_lensing",
          "doppler_localgp", "localgp_doppler",
          "doppler_integratedgp", "integratedgp_doppler",
          "lensing_localgp", "localgp_lensing",
          "lensing_integratedgp", "integratedgp_lensing",
          "localgp_integratedgp", "integratedgp_localgp",
     ]

The names of the GR effects implemented.
Their order is associated with the one in `VEC_ξs_GNC`, so be careful
to change it.

See also: [`VEC_ξs_GNC`](@ref)
"""
const GR_EFFECTS_GNC = [
     "auto_newton", "auto_doppler", "auto_lensing",
     "auto_localgp", "auto_integratedgp", "newton_doppler", "doppler_newton",
     "newton_lensing", "lensing_newton",
     "newton_localgp", "localgp_newton",
     "newton_integratedgp", "integratedgp_newton",
     "lensing_doppler", "doppler_lensing",
     "doppler_localgp", "localgp_doppler",
     "doppler_integratedgp", "integratedgp_doppler",
     "lensing_localgp", "localgp_lensing",
     "lensing_integratedgp", "integratedgp_lensing",
     "localgp_integratedgp", "integratedgp_localgp",
];


"""
     const VEC_ξs_GNC = [
          ξ_GNC_Newtonian, ξ_GNC_Doppler, ξ_GNC_Lensing,
          ξ_GNC_LocalGP, ξ_GNC_IntegratedGP, 
          
          ξ_GNC_Newtonian_Doppler, ξ_GNC_Doppler_Newtonian,
          ξ_GNC_Newtonian_Lensing, ξ_GNC_Lensing_Newtonian,
          ξ_GNC_Newtonian_LocalGP, ξ_GNC_LocalGP_Newtonian,
          ξ_GNC_Newtonian_IntegratedGP, ξ_GNC_IntegratedGP_Newtonian,
          ξ_GNC_Lensing_Doppler, ξ_GNC_Doppler_Lensing,
          ξ_GNC_Doppler_LocalGP, ξ_GNC_LocalGP_Doppler,
          ξ_GNC_Doppler_IntegratedGP, ξ_GNC_IntegratedGP_Doppler,
          ξ_GNC_Lensing_LocalGP, ξ_GNC_LocalGP_Lensing,
          ξ_GNC_Lensing_IntegratedGP, ξ_GNC_IntegratedGP_Lensing,
          ξ_GNC_LocalGP_IntegratedGP, ξ_GNC_IntegratedGP_LocalGP
     ]

The names of the GR effect TPCFs implemented.
Their order is associated with the one in `GR_EFFECTS_GNC`, so be careful
to change it.

See also: [`GR_EFFECTS_GNC`](@ref)
"""
const VEC_ξs_GNC = [
     ξ_GNC_Newtonian, ξ_GNC_Doppler, ξ_GNC_Lensing,
     ξ_GNC_LocalGP, ξ_GNC_IntegratedGP, ξ_GNC_Newtonian_Doppler, ξ_GNC_Doppler_Newtonian,
     ξ_GNC_Newtonian_Lensing, ξ_GNC_Lensing_Newtonian,
     ξ_GNC_Newtonian_LocalGP, ξ_GNC_LocalGP_Newtonian,
     ξ_GNC_Newtonian_IntegratedGP, ξ_GNC_IntegratedGP_Newtonian,
     ξ_GNC_Lensing_Doppler, ξ_GNC_Doppler_Lensing,
     ξ_GNC_Doppler_LocalGP, ξ_GNC_LocalGP_Doppler,
     ξ_GNC_Doppler_IntegratedGP, ξ_GNC_IntegratedGP_Doppler,
     ξ_GNC_Lensing_LocalGP, ξ_GNC_LocalGP_Lensing,
     ξ_GNC_Lensing_IntegratedGP, ξ_GNC_IntegratedGP_Lensing,
     ξ_GNC_LocalGP_IntegratedGP, ξ_GNC_IntegratedGP_LocalGP
];



"""
     const DICT_GR_ξs_GNC::Dict{String,Function}

For an input key string `effect` from `GR_EFFECTS_GNC`, return the 
associated TPCF `DICT_GR_ξs_GNC[effect]` from `VEC_ξs_GNC`.

## Example

```jldoctest
julia> GaPSE.DICT_GR_ξs_GNC["auto_doppler"]
ξ_GNC_Doppler
```

See also: [`GR_EFFECTS_GNC`](@ref), [`VEC_ξs_GNC`](@ref)
"""
const DICT_GR_ξs_GNC = Dict([k => v for (k, v) in zip(GR_EFFECTS_GNC, VEC_ξs_GNC)]...)


"""
     const INDEX_GR_EFFECT_GNC::Dict{String,Integer}

For an input key string `effect` from `GR_EFFECTS_GNC`, return the 
associated index position in that vector.

## Example

```jldoctest
julia> GaPSE.INDEX_GR_EFFECT_GNC["auto_doppler"]
1
```

See also: [`GR_EFFECTS_GNC`](@ref)
"""
const INDEX_GR_EFFECT_GNC = Dict([name => i for (i, name) in
                                  enumerate(GR_EFFECTS_GNC)]...)


"""
     const GR_EFFECT_INDEX_GNC::Dict{Integer,String}

For an input index position `i` of `GR_EFFECTS_GNC`, return the 
associated key string `effect`.

## Example

```jldoctest
julia> GaPSE.DICT_GR_ξs_GNC[1]
"auto_doppler"
```

See also: [`GR_EFFECTS_GNC`](@ref)
"""
const GR_EFFECT_INDEX_GNC = Dict([i => name for (i, name) in
                                  enumerate(GR_EFFECTS_GNC)]...)




##################### GALAXY NUMBER COUNTS x LUMINOSITY DISTANCE #############################




"""
     const GR_EFFECTS_GNCxLD = [
          "newton_doppler", 
          "newton_lensing", 
          "newton_localgp", 
          "newton_integratedgp",

          "doppler_doppler",
          "doppler_lensing",
          "doppler_localgp", 
          "doppler_integratedgp",

          "lensing_doppler",
          "lensing_lensing",
          "lensing_localgp",
          "lensing_integratedgp",
          
          "localgp_doppler",
          "localgp_lensing",
          "localgp_localgp",
          "localgp_integratedgp",

          "integratedgp_doppler",
          "integratedgp_lensing",
          "integratedgp_localgp",
          "integratedgp_integratedgp",
     ]

The names of the GR effects implemented.
Their order is associated with the one in `VEC_ξs_GNCxLD`, so be careful
to change it.

See also: [`VEC_ξs_GNCxLD`](@ref)
"""
const GR_EFFECTS_GNCxLD = [
     "newton_doppler", 
     "newton_lensing", 
     "newton_localgp", 
     "newton_integratedgp",

     "doppler_doppler",
     "doppler_lensing",
     "doppler_localgp", 
     "doppler_integratedgp",

     "lensing_doppler",
     "lensing_lensing",
     "lensing_localgp",
     "lensing_integratedgp",
     
     "localgp_doppler",
     "localgp_lensing",
     "localgp_localgp",
     "localgp_integratedgp",

     "integratedgp_doppler",
     "integratedgp_lensing",
     "integratedgp_localgp",
     "integratedgp_integratedgp",
];


"""
     const GR_EFFECTS_LDxGNC = [
          "doppler_newton", 
          "lensing_newton", 
          "localgp_newton", 
          "integratedgp_newton",

          "doppler_doppler",
          "lensing_doppler",
          "localgp_doppler", 
          "integratedgp_doppler",

          "doppler_lensing",
          "lensing_lensing",
          "localgp_lensing",
          "integratedgp_lensing",
          
          "doppler_localgp",
          "lensing_localgp",
          "localgp_localgp",
          "integratedgp_localgp",

          "doppler_integratedgp",
          "lensing_integratedgp",
          "localgp_integratedgp",
          "integratedgp_integratedgp",
     ]

The names of the GR effects implemented.
Their order is associated with the one in `VEC_ξs_LDxGNC`, so be careful
to change it.

See also: [`VEC_ξs_LDxGNC`](@ref)
"""
const GR_EFFECTS_LDxGNC = [
     "doppler_newton", 
     "lensing_newton", 
     "localgp_newton", 
     "integratedgp_newton",

     "doppler_doppler",
     "lensing_doppler",
     "localgp_doppler", 
     "integratedgp_doppler",

     "doppler_lensing",
     "lensing_lensing",
     "localgp_lensing",
     "integratedgp_lensing",
     
     "doppler_localgp",
     "lensing_localgp",
     "localgp_localgp",
     "integratedgp_localgp",

     "doppler_integratedgp",
     "lensing_integratedgp",
     "localgp_integratedgp",
     "integratedgp_integratedgp",
];


##########


"""
     const VEC_ξs_GNCxLD = [
          ξ_GNCxLD_Newtonian_Doppler,
          ξ_GNCxLD_Newtonian_Lensing, 
          ξ_GNCxLD_Newtonian_LocalGP, 
          ξ_GNCxLD_Newtonian_IntegratedGP, 

          ξ_GNCxLD_Doppler_Doppler,
          ξ_GNCxLD_Doppler_Lensing,
          ξ_GNCxLD_Doppler_LocalGP, 
          ξ_GNCxLD_Doppler_IntegratedGP, 
          
          ξ_GNCxLD_Lensing_Doppler,
          ξ_GNCxLD_Lensing_Lensing,
          ξ_GNCxLD_Lensing_LocalGP,
          ξ_GNCxLD_Lensing_IntegratedGP,

          ξ_GNCxLD_LocalGP_Doppler,
          ξ_GNCxLD_LocalGP_Lensing,
          ξ_GNCxLD_LocalGP_LocalGP,
          ξ_GNCxLD_LocalGP_IntegratedGP,

          ξ_GNCxLD_IntegratedGP_Doppler,
          ξ_GNCxLD_IntegratedGP_Lensing,
          ξ_GNCxLD_IntegratedGP_LocalGP,
          ξ_GNCxLD_IntegratedGP_IntegratedGP,
     ]

The names of the GR effect TPCFs implemented.
Their order is associated with the one in `GR_EFFECTS_GNCxLD`, so be careful
to change it.

See also: [`GR_EFFECTS_GNCxLD`](@ref)
"""
const VEC_ξs_GNCxLD = [
     ξ_GNCxLD_Newtonian_Doppler,
     ξ_GNCxLD_Newtonian_Lensing, 
     ξ_GNCxLD_Newtonian_LocalGP, 
     ξ_GNCxLD_Newtonian_IntegratedGP, 

     ξ_GNCxLD_Doppler_Doppler,
     ξ_GNCxLD_Doppler_Lensing,
     ξ_GNCxLD_Doppler_LocalGP, 
     ξ_GNCxLD_Doppler_IntegratedGP, 
     
     ξ_GNCxLD_Lensing_Doppler,
     ξ_GNCxLD_Lensing_Lensing,
     ξ_GNCxLD_Lensing_LocalGP,
     ξ_GNCxLD_Lensing_IntegratedGP,

     ξ_GNCxLD_LocalGP_Doppler,
     ξ_GNCxLD_LocalGP_Lensing,
     ξ_GNCxLD_LocalGP_LocalGP,
     ξ_GNCxLD_LocalGP_IntegratedGP,

     ξ_GNCxLD_IntegratedGP_Doppler,
     ξ_GNCxLD_IntegratedGP_Lensing,
     ξ_GNCxLD_IntegratedGP_LocalGP,
     ξ_GNCxLD_IntegratedGP_IntegratedGP,
];


"""
     const VEC_ξs_LDxGNC = [
          ξ_LDxGNC_Doppler_Newtonian,
          ξ_LDxGNC_Lensing_Newtonian,
          ξ_LDxGNC_LocalGP_Newtonian,
          ξ_LDxGNC_IntegratedGP_Newtonian,

          ξ_LDxGNC_Doppler_Doppler,
          ξ_LDxGNC_Lensing_Doppler,
          ξ_LDxGNC_LocalGP_Doppler,
          ξ_LDxGNC_IntegratedGP_Doppler,

          ξ_LDxGNC_Doppler_Lensing,
          ξ_LDxGNC_Lensing_Lensing,
          ξ_LDxGNC_LocalGP_Lensing,
          ξ_LDxGNC_IntegratedGP_Lensing,

          ξ_LDxGNC_Doppler_LocalGP,
          ξ_LDxGNC_Lensing_LocalGP,
          ξ_LDxGNC_LocalGP_LocalGP,
          ξ_LDxGNC_IntegratedGP_LocalGP,

          ξ_LDxGNC_Doppler_IntegratedGP,
          ξ_LDxGNC_Lensing_IntegratedGP,
          ξ_LDxGNC_LocalGP_IntegratedGP,
          ξ_LDxGNC_IntegratedGP_IntegratedGP,
     ]

The names of the GR effect TPCFs implemented.
Their order is associated with the one in `GR_EFFECTS_LDxGNC`, so be careful
to change it.

See also: [`GR_EFFECTS_LDxGNC`](@ref)
"""
const VEC_ξs_LDxGNC = [
     ξ_LDxGNC_Doppler_Newtonian,
     ξ_LDxGNC_Lensing_Newtonian,
     ξ_LDxGNC_LocalGP_Newtonian,
     ξ_LDxGNC_IntegratedGP_Newtonian,

     ξ_LDxGNC_Doppler_Doppler,
     ξ_LDxGNC_Lensing_Doppler,
     ξ_LDxGNC_LocalGP_Doppler,
     ξ_LDxGNC_IntegratedGP_Doppler,

     ξ_LDxGNC_Doppler_Lensing,
     ξ_LDxGNC_Lensing_Lensing,
     ξ_LDxGNC_LocalGP_Lensing,
     ξ_LDxGNC_IntegratedGP_Lensing,

     ξ_LDxGNC_Doppler_LocalGP,
     ξ_LDxGNC_Lensing_LocalGP,
     ξ_LDxGNC_LocalGP_LocalGP,
     ξ_LDxGNC_IntegratedGP_LocalGP,

     ξ_LDxGNC_Doppler_IntegratedGP,
     ξ_LDxGNC_Lensing_IntegratedGP,
     ξ_LDxGNC_LocalGP_IntegratedGP,
     ξ_LDxGNC_IntegratedGP_IntegratedGP,
];


##########


"""
     const DICT_GR_ξs_GNCxLD::Dict{String,Function}

For an input key string `effect` from `GR_EFFECTS_GNCxLD`, return the 
associated TPCF `DICT_GR_ξs_GNCxLD[effect]` from `VEC_ξs_GNCxLD`.

## Example

```jldoctest
julia> GaPSE.DICT_GR_ξs_GNCxLD["lensing_doppler"]
ξ_GNCxLD_Lensing_Doppler
```

See also: [`GR_EFFECTS_GNCxLD`](@ref), [`VEC_ξs_GNCxLD`](@ref)
"""
const DICT_GR_ξs_GNCxLD = Dict([k => v for (k, v) in zip(GR_EFFECTS_GNCxLD, VEC_ξs_GNCxLD)]...)


"""
     const DICT_GR_ξs_LDxGNC::Dict{String,Function}

For an input key string `effect` from `GR_EFFECTS_LDxGNC`, return the 
associated TPCF `DICT_GR_ξs_LDxGNC[effect]` from `VEC_ξs_v`.

## Example

```jldoctest
julia> GaPSE.DICT_GR_ξs_LDxGNC["lensing_doppler"]
ξ_LDxGNC_Lensing_Doppler
```

See also: [`GR_EFFECTS_LDxGNC`](@ref), [`VEC_ξs_LDxGNC`](@ref)
"""
const DICT_GR_ξs_LDxGNC = Dict([k => v for (k, v) in zip(GR_EFFECTS_LDxGNC, VEC_ξs_LDxGNC)]...)


##########


"""
     const INDEX_GR_EFFECT_GNCxLD::Dict{String,Integer}

For an input key string `effect` from `GR_EFFECTS_GNCxLD`, return the 
associated index position in that vector.

## Example

```jldoctest
julia> GaPSE.INDEX_GR_EFFECT_GNCxLD["newton_lensing"]
2
```

See also: [`GR_EFFECTS_GNCxLD`](@ref)
"""
const INDEX_GR_EFFECT_GNCxLD = Dict([name => i for (i, name) in
                                  enumerate(GR_EFFECTS_GNCxLD)]...)


"""
     const INDEX_GR_EFFECT_LDxGNC::Dict{String,Integer}

For an input key string `effect` from `GR_EFFECTS_LDxGNC`, return the 
associated index position in that vector.

## Example

```jldoctest
julia> GaPSE.INDEX_GR_EFFECT_LDxGNC["lensing_newton"]
2
```

See also: [`GR_EFFECTS_LDxGNC`](@ref)
"""
const INDEX_GR_EFFECT_LDxGNC = Dict([name => i for (i, name) in
                                  enumerate(GR_EFFECTS_LDxGNC)]...)


##########


"""
     const GR_EFFECT_INDEX_GNCxLD::Dict{Integer,String}

For an input index position `i` of `GR_EFFECTS_GNCxLD`, return the 
associated key string `effect`.

## Example

```jldoctest
julia> GaPSE.DICT_GR_ξs_GNCxLD[2]
"newton_doppler"
```

See also: [`GR_EFFECTS_GNCxLD`](@ref)
"""
const GR_EFFECT_INDEX_GNCxLD = Dict([i => name for (i, name) in
                                  enumerate(GR_EFFECTS_GNCxLD)]...)


"""
     const GR_EFFECT_INDEX_LDxGNC::Dict{Integer,String}

For an input index position `i` of `GR_EFFECTS_LDxGNC`, return the 
associated key string `effect`.

## Example

```jldoctest
julia> GaPSE.DICT_GR_ξs_LDxGNC[2]
"doppler_newton"
```

See also: [`GR_EFFECTS_LDxGNC`](@ref)
"""
const GR_EFFECT_INDEX_LDxGNC = Dict([i => name for (i, name) in
                                  enumerate(GR_EFFECTS_LDxGNC)]...)

