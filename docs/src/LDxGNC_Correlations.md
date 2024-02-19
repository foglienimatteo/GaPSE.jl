```@meta
DocTestSetup = quote
    using GaPSE
end
```

# LDxGNC TPCFs

## Two-Point Cross-Correlation Functions

```@docs
GaPSE.ξ_LDxGNC_Doppler_Newtonian
GaPSE.ξ_LDxGNC_Lensing_Newtonian
GaPSE.ξ_LDxGNC_LocalGP_Newtonian
GaPSE.ξ_LDxGNC_IntegratedGP_Newtonian
GaPSE.ξ_LDxGNC_Doppler_Doppler
GaPSE.ξ_LDxGNC_Lensing_Doppler
GaPSE.ξ_LDxGNC_LocalGP_Doppler
GaPSE.ξ_LDxGNC_IntegratedGP_Doppler
GaPSE.ξ_LDxGNC_Doppler_Lensing
GaPSE.ξ_LDxGNC_Lensing_Lensing
GaPSE.ξ_LDxGNC_LocalGP_Lensing
GaPSE.ξ_LDxGNC_IntegratedGP_Lensing
GaPSE.ξ_LDxGNC_Doppler_LocalGP
GaPSE.ξ_LDxGNC_Lensing_LocalGP
GaPSE.ξ_LDxGNC_LocalGP_LocalGP
GaPSE.ξ_LDxGNC_IntegratedGP_LocalGP
GaPSE.ξ_LDxGNC_Doppler_IntegratedGP
GaPSE.ξ_LDxGNC_Lensing_IntegratedGP
GaPSE.ξ_LDxGNC_LocalGP_IntegratedGP
GaPSE.ξ_LDxGNC_IntegratedGP_IntegratedGP
```


## Two-Point Cross-Correlation Function multipoles

```@docs
GaPSE.integrand_ξ_LDxGNC_multipole
GaPSE.ξ_LDxGNC_multipole
GaPSE.map_ξ_LDxGNC_multipole
GaPSE.print_map_ξ_LDxGNC_multipole
```

## Two-Point Cross-Correlation Function Sum multipoles

```@docs
GaPSE.sum_ξ_LDxGNC_multipole
GaPSE.map_sum_ξ_LDxGNC_multipole
GaPSE.print_map_sum_ξ_LDxGNC_multipole
```
