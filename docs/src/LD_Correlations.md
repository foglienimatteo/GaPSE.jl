```@meta
DocTestSetup = quote
    using GaPSE
end
```

# LD TPCFs

## Two-Point Auto-Correlation Functions

```@docs
GaPSE.ξ_LD_Doppler
GaPSE.ξ_LD_Lensing
GaPSE.ξ_LD_LocalGP
GaPSE.ξ_LD_IntegratedGP
```

## Two-Point Cross-Correlation Functions

```@docs
GaPSE.ξ_LD_Doppler_Lensing
GaPSE.ξ_LD_Lensing_Doppler
GaPSE.ξ_LD_Doppler_LocalGP
GaPSE.ξ_LD_LocalGP_Doppler
GaPSE.ξ_LD_Doppler_IntegratedGP
GaPSE.ξ_LD_IntegratedGP_Doppler
GaPSE.ξ_LD_Lensing_LocalGP
GaPSE.ξ_LD_LocalGP_Lensing
GaPSE.ξ_LD_Lensing_IntegratedGP
GaPSE.ξ_LD_IntegratedGP_Lensing
GaPSE.ξ_LD_LocalGP_IntegratedGP
GaPSE.ξ_LD_IntegratedGP_LocalGP
```


## Two-Point Cross-Correlation Function multipoles

```@docs
GaPSE.integrand_ξ_LD_multipole
GaPSE.ξ_LD_multipole
GaPSE.map_ξ_LD_multipole
GaPSE.print_map_ξ_LD_multipole
```

## Two-Point Cross-Correlation Function Sum multipoles

```@docs
GaPSE.sum_ξ_LD_multipole
GaPSE.map_sum_ξ_LD_multipole
GaPSE.print_map_sum_ξ_LD_multipole
```
