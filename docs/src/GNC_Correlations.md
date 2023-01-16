```@meta
DocTestSetup = quote
    using GaPSE
end
```

# Two-Point Auto-Correlation Functions

```@docs
GaPSE.ξ_GNC_Newtonian
GaPSE.ξ_GNC_Doppler
GaPSE.ξ_GNC_Lensing
GaPSE.ξ_GNC_LocalGP
GaPSE.ξ_GNC_IntegratedGP

GaPSE.integrand_ξ_GNC_Lensing
GaPSE.integrand_ξ_GNC_IntegratedGP
```

# Two-Point Cross-Correlation Functions

```@docs
GaPSE.ξ_GNC_Newtonian_Doppler
GaPSE.ξ_GNC_Doppler_Newtonian
GaPSE.ξ_GNC_Newtonian_Lensing
GaPSE.ξ_GNC_Lensing_Newtonian
GaPSE.ξ_GNC_Newtonian_LocalGP
GaPSE.ξ_GNC_LocalGP_Newtonian
GaPSE.ξ_GNC_Newtonian_IntegratedGP
GaPSE.ξ_GNC_IntegratedGP_Newtonian
GaPSE.ξ_GNC_Doppler_Lensing
GaPSE.ξ_GNC_Lensing_Doppler
GaPSE.ξ_GNC_Doppler_LocalGP
GaPSE.ξ_GNC_LocalGP_Doppler
GaPSE.ξ_GNC_Doppler_IntegratedGP
GaPSE.ξ_GNC_IntegratedGP_Doppler
GaPSE.ξ_GNC_Lensing_LocalGP
GaPSE.ξ_GNC_LocalGP_Lensing
GaPSE.ξ_GNC_Lensing_IntegratedGP
GaPSE.ξ_GNC_IntegratedGP_Lensing
GaPSE.ξ_GNC_LocalGP_IntegratedGP
GaPSE.ξ_GNC_IntegratedGP_LocalGP

GaPSE.integrand_ξ_GNC_Newtonian_Lensing
GaPSE.integrand_ξ_GNC_Newtonian_IntegratedGP
GaPSE.integrand_ξ_GNC_Doppler_IntegratedGP
GaPSE.integrand_ξ_GNC_Lensing_Doppler
GaPSE.integrand_ξ_GNC_Lensing_LocalGP
GaPSE.integrand_ξ_GNC_Lensing_IntegratedGP
GaPSE.integrand_ξ_GNC_LocalGP_IntegratedGP
```


# Two-Point Cross-Correlation Function multipoles

```@docs
GaPSE.integrand_ξ_GNC_multipole
GaPSE.ξ_GNC_multipole
GaPSE.map_ξ_GNC_multipole
GaPSE.print_map_ξ_GNC_multipole
```

# Two-Point Cross-Correlation Function Sum multipoles

```@docs
GaPSE.sum_ξ_GNC_multipole
GaPSE.map_sum_ξ_GNC_multipole
GaPSE.print_map_sum_ξ_GNC_multipole
```
