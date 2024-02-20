```@meta
DocTestSetup = quote
    using GaPSE
end
```

# LD TPCFs

## Two-Point Auto-Correlation Functions integrands

```@docs
GaPSE.integrand_ξ_LD_Lensing
GaPSE.integrand_ξ_LD_IntegratedGP
```

## Two-Point Cross-Correlation Functions integrands

```@docs
GaPSE.integrand_ξ_LD_Lensing_LocalGP
GaPSE.integrand_ξ_LD_Doppler_IntegratedGP
GaPSE.integrand_ξ_LD_LocalGP_IntegratedGP
GaPSE.integrand_ξ_LD_Lensing_IntegratedGP
GaPSE.integrand_ξ_LD_Lensing_Doppler
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
