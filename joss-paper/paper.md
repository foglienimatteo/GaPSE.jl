---
title: 'GaPSE.jl: Full Relativistic two-point statistics of Galaxy Number Counts and Luminosity-Distance fluctuations'
tags:
  - Julia
  - cosmology
  - large-scale structure
  - galaxy clustering
  - general relativity
authors:
  - name: Matteo Foglieni
    orcid: 0000-0001-8736-4517
    affiliation: 1
  - name: Mattia Pantiri
    orcid: 0009-0007-2259-2368
    affiliation: 2
affiliations:
  - name: Leibniz Supercomputing Centre (LRZ), Boltzmannstr. 1, 85748 Garching bei Muenchen, Germany
    index: 1
  - name: Leiden Institute of Physics (LION), Niels Bohrweg 2, 2333 CA Leiden, The Netherlands
    index: 2
date: 14 October 2026
bibliography: paper.bib
---


# Summary

`GaPSE.jl` (Galaxy Power Spectrum Estimator) is an open-source Julia code written for computing 
the two-point statistics of relativistic effects in galaxy clustering at first order in perturbation theory. 
The software computes two-point correlation functions (TPCFs) of relativistic galaxy number counts (GNC) 
and luminosity-distance (LD) perturbations, including their cross-correlations, as well as their corresponding power spectra multipoles, 
with or without a survey window function. Its modular design provides individual contributions from the different physical effects, 
allowing users to investigate their relative importance across redshifts, multipoles, and survey geometries. 
The code supports multithreaded calculations, and an extensive documentation and many Jupyter notebooks are provided.


# Statement of need


The current standard paradigm for the origin of the structure formation in the Universe
is the inflation, i.e. an exponentially rapidly growth of the Universe scale factor $a$ of about
28 orders of magnitude. Almost all the inflationary models predict a flat universe and a nearly scale-invariant 
spectrum of initial fluctuations. Consequently, studying constraining variables that might distinguish between
these models is of utmost importance. Local Primordial Non-Gaussianities (PNG) are one of the most promising 
variable of this kind [@Achucarro:2022].

The standard parameterization of the local PNG is the so-called scale-independent $f_NL$ parameterization,
which would imprint a scale-dependent bias at low wavenumbers on the observable Power Spectrum [@Dalal:2008; @Slosar:2008].
However, on those scale, the pure Newtonian and plane-parallel approximations that underpin most clustering analyses cease to be adequate,
because the separation between a pair of galaxies becomes comparable to the distance from us to them.
Consequenly, a full General Relativistic (GR) approach must be employed.

In such analysis, the GR effects manifest themself into two observable groups, that can be derived analytically [@Castorina:2022]: 
the Luminosity Distance (LD) fluctuations and the Galaxy Number Counts (GNC).
Due to the presence of the gravitational potential in many GR effects, some of these terms might _mimic_ the PNG detection and measurements,
leading to wrong results. For this reason, we need a deeper understanding of these GR effects,
especially for future survey that aim to map the Universe on the largest scales, like Euclid.
This is why we developed `GaPSE.jl`.

`GaPSE.jl` (Galaxy Power Spectrum Estimator) is an open-source code, written in the Julia programming language, that computes the
two-point statistics of these GR effects at first order in perturbation theory.
This enables users to assess the importance of relativistic corrections, compare them with the PNG signal, 
and investigate how survey geometry affects the resulting observables.


# State of the field and research impact statement

Publicly available Boltzmann codes, such as CAMB [@CAMB:2026], compute the observed number counts as angular power
spectra $C_\ell(z_i, z_j)$, a representation that is exact but expensive to connect to the
Fourier-space quantities that clustering analyses actually fit.
`GaPSE.jl` implements instead the Fourier-space formalism of [@CastorinaDiDio:2022], 
in which the observed Power Spectrum multipoles are obtained from a configuration-space calculation that keeps the two
lines of sight exact. 
The code was written for, and used to produce the results of, [@Foglieni:2023] and [@Pantiri:2024], 
which respectively derived the large-scale limit of the observed galaxy Power Spectrum and extended the formalism to luminosity-distance
fluctuations.

The practical need it fills is the ability to isolate individual contributions. 
Every relativistic term is computed and returned separately, a user can ask which effect
dominates at a given redshift, at which multipole order, and under which survey geometry.
This term-by-term decomposition facilitates physical interpretation and systematic assessments of potential degeneracies with PNG. 
Such analyses would be considerably less straightforward with software that provides only the total observable.


# Software design

Giving in input the linear matter Power Spectrum at $z = 0$ and the background quantities of the 
chosen cosmology - both read directly from the output of `CLASS` [@CLASS:2011] - `GaPSE.jl` computes:

- the 25 Two-Point Correlation Functions (TPCFs) of the relativistic Galaxy Number Counts (GNC), at arbitrary multipole
  order, covering the auto- and cross-correlations of the density and redshift-space
  distortion ("Newtonian"), Doppler, lensing, local gravitational potential and integrated
  gravitational potential terms; observer terms can be optionally included;
- the 16 TPCFs of the Luminosity-Distance (LD) perturbations;
- the 20 GNC $\times$ LD cross-correlations and their 20 LD $\times$ GNC counterparts;
- the Power Spectrum multipoles of all of the above, based on the Yamamoto estimator;
- the Doppler and matter TPCFs in the plane-parallel approximation, for comparison;
- the Local PNG contribution, modelled via the scale-dependent bias

All of these can be evaluated with or without a survey window function; a toy survey with
azimuthal symmetry is provided, and an arbitrary window can be supplied as a tabulated map.

The numerical core of the code rests on an FFTLog evaluation of the $I_\ell^n$ integrals of the input
Power Spectrum, through `TwoFAST.jl` [@Gebhardt:2018] and `FFTLog.jl` [@Bonici:2021], and on three 
interchangeable quadrature schemes: trapezoidal, Gauss–Lobatto and adaptive Gauss–Kronrod.

The loop over input separation points is parallelized using Julia's multithreading capabilities. 
In the tested configuration, the multithreaded implementation achieves a speedup of 
approximately 3.3× using  threads. The numerical results are independent of the number of threads used.
The code repo ships a test suite that compares every correlation against stored reference values, 
run in continuous integration on Linux and macOS, single and multi-threaded.
The extensive documentation covers both the API and the theory behind each expression, and a set of
Jupyter notebooks reproduces all the figures and the numerical claims.
Furthermore, a Docker container with all the aforementioned features and files can be pulled from `matteofoglieni:gapse`.


# AI usage disclosure

No AI tools were used in the development of `GaPSE.jl` before and up to version `v0.8.0`.
Beginning with version `v0.8.0`, we took advantage of AI models developed by Anthropic
to investigate numerical issues and inconsistencies arising in specific corner cases and to improve the software documentation.
We take full responsibility for reviewing and validating the AI-assisted work and for the resulting code and documentation.


# References
