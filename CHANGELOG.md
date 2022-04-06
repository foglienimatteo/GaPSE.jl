# HEAD

# VERSION 0.3.0

- added the functions for the evaluation of the Doppler auto-CF with the Plane-Parallel approximation and their tests

- modified default range 10 .^ range(-1, 3, length= N_log) to  10 .^ range(0, 3, length= N_log) 

- added tests for PS evaluation


# VERSION 0.2.0

- Modified the `CosmoParams` struct: now you should pass dictionaries for `InputPS` and `IPSTools` options

- Added the `ipynbs/TUTORIAL.ipynb` notebook, which is a small tour of how the code should be used

- Changed the keyword argument `s_1` into `s1` (in functions such as `map_ξ_multipole`, `map_sum_ξ_multipole`, ...)

- Removed the functions `integrand_on_mu`, `integral_on_mu`, `map_integral_on_mu` and `print_map_integral_on_mu`; added `integrand_ξ_multipole`. This restructuration let the code more flexible.

- Optimized `F_map` and re-written its keyword argument structure

- Improved the documentation of the code, and the `README.md`



# VERSION 0.1.0

Created the basic structure of the program, with the structs `WindowF`, `InputPS`, `CosmoParams` and `Cosmology`.
It's possible to:
- read an arbitrary input matter power spectrum, produced by CLASS
- read an arbitrary input background data file, produced by CLASS
- write and read a map of the window F function
- calculate all the General Relativistic (GR) effects Two-Point Correlation Functions (TPCFs) multipole functions (both auto-correlation and cross-correlation) arising from the perturbed luminosity distance [[1]](#1) for an arbitrary multipole order 
- calculate the sum of all these TPCFs

Some notebook are already provided, inside the directory `ipynb`, such as:
- `ALL_CF_L.ipynb`
- `ALL_Multipoles.ipynb`
- `PS_Multipoles.ipynb`
- `SUM_CF.ipynb`
  


## References

<a id="1">[1]</a> 
See the equation (C.23) with ``s`` Dalal, Doré et al., _Imprints of primordial non-Gaussianities on large-scale structure_ (2008), American Physical Society, DOI: 10.1103/PhysRevD.77.123514, 
url: https://journals.aps.org/prd/abstract/10.1103/PhysRevD.77.123514
