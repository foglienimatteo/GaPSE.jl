```@meta
DocTestSetup = quote
    using GaPSE
end
```

# GNC TPCFs

## Two-Point Cross-Correlation Functions integrands

```@docs
GaPSE.integrand_ξ_GNC_Newtonian_Lensing
GaPSE.integrand_ξ_GNC_Newtonian_IntegratedGP
GaPSE.integrand_ξ_GNC_Doppler_IntegratedGP
GaPSE.integrand_ξ_GNC_Lensing_Doppler
GaPSE.integrand_ξ_GNC_Lensing_LocalGP
GaPSE.integrand_ξ_GNC_Lensing_IntegratedGP
GaPSE.integrand_ξ_GNC_LocalGP_IntegratedGP
```


## Two-Point Cross-Correlation Function multipoles

```@docs
GaPSE.integrand_ξ_GNC_multipole
GaPSE.ξ_GNC_multipole
GaPSE.map_ξ_GNC_multipole
GaPSE.print_map_ξ_GNC_multipole
```

## Two-Point Cross-Correlation Function Sum multipoles

```@docs
GaPSE.sum_ξ_GNC_multipole
GaPSE.map_sum_ξ_GNC_multipole
GaPSE.print_map_sum_ξ_GNC_multipole
```