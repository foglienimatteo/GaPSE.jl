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


function ξ_GNC_Newtonian_Doppler(P1::Point, P2::Point, y, cosmo::Cosmology; 
    b1=nothing, b2=nothing, s_b1=nothing, s_b2=nothing, 𝑓_evo1=nothing, 𝑓_evo2=nothing,
    s_lim=nothing, obs::Union{Bool, Symbol} = :noobsvel)

    s1, D1, f1 = P1.comdist, P1.D, P1.f
    s2, D2, f2, ℋ2 = P2.comdist, P2.D, P2.f, P2.ℋ

    b1 = isnothing(b1) ? cosmo.params.b1 : b1
    s_b2 = isnothing(s_b2) ? cosmo.params.s_b2 : s_b2
    𝑓_evo2 = isnothing(𝑓_evo2) ? cosmo.params.𝑓_evo2 : 𝑓_evo2

    s_lim = isnothing(s_lim) ? cosmo.params.s_lim : s_lim
    ℛ2 = func_ℛ_GNC(s2, P2.ℋ, P2.ℋ_p; s_b=s_b2, 𝑓_evo=𝑓_evo2, s_lim=s_lim)

    Δs = s(s1, s2, y)

    common = D1 * D2 * f2 * ℋ2 * ℛ2

    J00 = 1 / 15 * (5 * b1 * (s2 - y * s1) + f1 * (2 * y^2 * s2 - 3 * y * s1 + s2))
    J02 = 1 / (21 * Δs^2) * (
        7 * b1 * (y * s1 - s2) * (2 * y * s1 * s2 - s1^2 - s2^2) +
        f1 * (
            (10 * y^2 - 1) * s1^2 * s2
            -
            y * (5 * y^2 + 4) * s1 * s2^2
            +
            (y^2 + 2) * s2^3 - 3 * y * s1^3
        )
    )
    J04 = 1 / (35 * Δs^2) * f1 * (
            -2 * (y^2 + 2) * s1^2 * s2
            + y * (y^2 + 5) * s1 * s2^2
            + (1 - 3 * y^2) * s2^3 + 2 * y * s1^3
        )

    I00 = cosmo.tools.I00(Δs)
    I20 = cosmo.tools.I20(Δs)
    I40 = cosmo.tools.I40(Δs)

    if obs == false || obs == :no || obs == :noobsvel
        return common * (J00 * I00 + J02 * I20 + J04 * I40)
        
    elseif obs == true || obs == :yes
        #### New observer terms #########

        I31_s1 = cosmo.tools.I31(s1)
        I11_s1 = cosmo.tools.I11(s1)

        obs_common = 1 / 5 * y * f0 * ℋ0 * s1 * (ℛ2 - 5 * s_b2 + 2)

        obs_terms = D1 * obs_common * ((5 * b1 + 3 * f1) * I11_s1 - 2 * f1 * I31_s1)

        #################################

        return common * (J00 * I00 + J02 * I20 + J04 * I40) + obs_terms
    else
        throw(AssertionError(":$obs is not a valid Symbol for \"obs\"; they are: \n\t"*
            "$(":".*string.(VALID_OBS_VALUES) .* vcat([" , " for i in 1:length(VALID_OBS_VALUES)-1], " .")... )" 
            ))
    end
end


function ξ_GNC_Newtonian_Doppler(s1, s2, y, cosmo::Cosmology; kwargs...)
    P1, P2 = Point(s1, cosmo), Point(s2, cosmo)
    return ξ_GNC_Newtonian_Doppler(P1, P2, y, cosmo; kwargs...)
end



"""
    ξ_GNC_Newtonian_Doppler(
        P1::Point, P2::Point, y, cosmo::Cosmology; 
        b1=nothing, b2=nothing, s_b1=nothing, s_b2=nothing, 
        𝑓_evo1=nothing, 𝑓_evo2=nothing,
        s_lim=nothing, obs::Union{Bool, Symbol} = :noobsvel
        ) ::Float64

    ξ_GNC_Newtonian_Doppler(s1, s2, y, cosmo::Cosmology;
        kwargs...) ::Float64

Return the Two-Point Correlation Function (TPCF) given by the cross correlation between the 
Newtonian and the Doppler effects arising from the Galaxy Number Counts (GNC).

In the first method, you should pass the two `Point` (`P1` and `P2`) where to 
evaluate the function, while in the second method (that internally recalls the first) 
you must provide the two corresponding comoving distances `s1` and `s2`.
We remember that all the distances are measured in ``h_0^{-1}\\mathrm{Mpc}``.

The analytical expression of this TPCF is the following:

```math
\\begin{split}
    \\xi^{\\delta v_{\\parallel}}( s_1 , s_2, y ) &= 
    D_1 D_2 J^{\\delta v_{\\parallel}}_{\\alpha} \\left[ 
        J^{\\delta v_{\\parallel}}_{00} I^0_{0} (s) + 
        J^{\\delta v_{\\parallel}}_{02} I^0_2 (s) +
        J^{\\delta v_{\\parallel}}_{04} I^0_4 (s) 
        \\right]\\nonumber \\\\
    &+ D_1 J^{\\delta v_{\\parallel}}_{\\beta} \\left[ 
        J^{\\delta v_{\\parallel}}_{11} I^1_{1} (s_1) + 
        J^{\\delta v_{\\parallel}}_{13} I^1_{3} (s_1) 
    \\right]
\\end{split}
```

with

```math
\\begin{split}
    J^{\\delta v_{\\parallel}}_{\\alpha} &= f_2 \\mathcal{H}_2 \\mathcal{R}_2
    \\, , \\\\
    %%%%%%%%%%%%%%%%%%%%%
    J^{\\delta v_{\\parallel}}_{\\beta} &= 
    y f_0 \\mathcal{H}_0 s_1 (\\mathcal{R}_2 - 5s_{\\mathrm{b}, 2}  + 2)
    \\, , \\\\
    %%%%%%%%%%%%%%%%%%%%%
    J^{\\delta v_{\\parallel}}_{00} &=
    \\frac{1}{15} \\left\\{
        s_2 \\left[ 5 b_1 + (2 y^2 + 1) f_1 \\right] -
        y s_1 \\left[ 3 f_1 + 5 b_1 \\right]
    \\right\\} 
    \\, , \\\\
    %%%%%%%%%%%%%%%%%%%%%
    J^{\\delta v_{\\parallel}}_{02} &=
    \\frac{1}{21 s^2} 
    \\left\\{ 
        \\left[
            (y^2 + 1) f_1 + 7 b_1
        \\right] s_2^3 -
        y \\left[
            21 b_1 + (5 y^2 + 4) f_1 
        \\right] s_1 s_2^2 +
        \\right.\\nonumber \\\\
        &\\left. \\qquad \\qquad
        \\left[
            7 (2 y^2 + 1) b_1 + (10 y^2 - 1) f_1
        \\right] s_1^2 s_2 -
        y \\left[
            7 b_1 + 3 f_1
        \\right] s_1^3
    \\right\\}
    \\, , \\\\
    %%%%%%%%%%%%%%%%%%%%%
    J^{\\delta v_{\\parallel}}_{04} &=
    \\frac{f_1}{35 s^2} 
    \\left[
        2 y s_1 ^3
        - 2 (y^2 + 2) s_1 ^2 s_2 +
        y (y^2 + 5) s_1 s_2^2 +
        (1 - 3 y^2) s_2 ^3
    \\right] 
    \\, , \\\\
    %%%%%%%%%%%%%%%%%%%%%
    J^{\\delta v_{\\parallel}}_{11} &= \\frac{1}{5} (5b_1+3 f_1 )  
    \\, , \\\\
    %%%%%%%%%%%%%%%%%%%%%
    J^{\\delta v_{\\parallel}}_{13} &= -\\frac{2}{5} f_1 
    \\, ,
\\end{split}
```

where:

- ``s_1`` and ``s_2`` are comoving distances;

- ``D_1 = D(s_1)``, ... is the linear growth factor (evaluated in ``s_1``);

- ``a_1 = a(s_1)``, ... is the scale factor (evaluated in ``s_1``);

- ``f_1 = f(s_1)``, ... is the linear growth rate (evaluated in ``s_1``);

- ``\\mathcal{H}_1 = \\mathcal{H}(s_1)``, ... is the comoving 
  Hubble parameter (evaluated in ``s_1``, ...);

- ``y = \\cos{\\theta} = \\hat{\\mathbf{s}}_1 \\cdot \\hat{\\mathbf{s}}_2``;

- ``\\mathcal{R}_1 = \\mathcal{R}(s_1)``, ... is 
  computed by `func_ℛ_GNC` in `cosmo::Cosmology` (and evaluated in ``s_1`` );
  the definition of ``\\mathcal{R}(s)`` is the following:
  ```math
  \\mathcal{R}(s) = 5 s_{\\mathrm{b}}(s) + \\frac{2 - 5 s_{\\mathrm{b}}(s)}{\\mathcal{H}(s) \\, s} +  
  \\frac{\\dot{\\mathcal{H}}(s)}{\\mathcal{H}(s)^2} - \\mathit{f}_{\\mathrm{evo}} \\quad ;
  ```

- ``b_1``, ``s_{\\mathrm{b}, 1}``, ``\\mathit{f}_{\\mathrm{evo}, 1}`` 
  (and ``b_2``, ``s_{\\mathrm{b}, 2}``, ``\\mathit{f}_{\\mathrm{evo}, 2}``) : 
  galaxy bias, magnification bias (i.e. the slope of the luminosity function at the luminosity threshold), 
  and evolution bias for the first (second) effect;

- ``\\Omega_{\\mathrm{M}0} = \\Omega_{\\mathrm{cdm}} + \\Omega_{\\mathrm{b}}`` is the sum of 
  cold-dark-matter and barionic density parameters (again, stored in `cosmo`);

- ``I_\\ell^n`` and ``\\sigma_i`` are defined as
  ```math
  I_\\ell^n(s) = \\int_0^{+\\infty} \\frac{\\mathrm{d}q}{2\\pi^2} 
  \\, q^2 \\, P(q) \\, \\frac{j_\\ell(qs)}{(qs)^n} \\quad , 
  \\quad \\sigma_i = \\int_0^{+\\infty} \\frac{\\mathrm{d}q}{2\\pi^2} 
  \\, q^{2-i} \\, P(q)
  ```
  with ``P(q)`` as the matter Power Spectrum at ``z=0`` (stored in `cosmo`) 
  and ``j_\\ell`` as spherical Bessel function of order ``\\ell``;

- ``\\tilde{I}_0^4`` is defined as
  ```math
  \\tilde{I}_0^4 = \\int_0^{+\\infty} \\frac{\\mathrm{d}q}{2\\pi^2} 
  \\, q^2 \\, P(q) \\, \\frac{j_0(qs) - 1}{(qs)^4}
  ``` 
  with ``P(q)`` as the matter Power Spectrum at ``z=0`` (stored in `cosmo`) 
  and ``j_\\ell`` as spherical Bessel function of order ``\\ell``;

- ``\\mathcal{H}_0``, ``f_0`` and so on are evaluated at the observer position (i.e. at present day);

- ``\\Delta\\chi_1 := \\sqrt{\\chi_1^2 + s_2^2-2\\,\\chi_1\\,s_2\\,y}`` and 
  ``\\Delta\\chi_2 := \\sqrt{s_1^2 + \\chi_2^2-2\\,s_1\\,\\chi_2\\,y}``;

- ``s=\\sqrt{s_1^2 + s_2^2 - 2 \\, s_1 \\, s_2 \\, y}`` and 
  ``\\Delta\\chi := \\sqrt{\\chi_1^2 + \\chi_2^2-2\\,\\chi_1\\,\\chi_2\\,y}``.

In this TPCF, the term proportional to only``D(s_1)`` is an observer term 
(while the term proportional to ``D(s_1) \\, D(s_2)`` is not) and it does depend 
on the observer velocity. Consequently, if you set `obs = :yes` or `obs = true`
both of them will computed, while for `obs = :no`, `obs = false` or
even `obs = :noobsvel` only the first one will be taken into account.

## Inputs

- `P1::Point` and `P2::Point`, or `s1` and `s2`: `Point`/comoving distances where the 
  TPCF has to be calculated; they contain all the 
  data of interest needed for this calculus (comoving distance, growth factor and so on);
  
- `y`: the cosine of the angle between the two points `P1` and `P2` wrt the observer;

- `cosmo::Cosmology`: cosmology to be used in this computation; it contains all the splines
  used for the conversion `s` -> `Point`, and all the cosmological parameters ``b``, ...

## Keyword Arguments

- `b1=nothing`, `s_b1=nothing`, `𝑓_evo1=nothing` and `b2=nothing`, `s_b2=nothing`, `𝑓_evo2=nothing`:
  galaxy, magnification and evolutionary biases respectively for the first and the second effect 
  computed in this TPCF:
  - if not set (i.e. if you leave the default value `nothing`) the values stored in the input `cosmo`
    will be considered;
  - if you set one or more values, they will override the `cosmo` ones in this computation;
  - the two sets of values should be different only if you are interested in studing two galaxy species;
  - only the required parameters for the chosen TPCF will be used, depending on its analytical expression;
    all the others will have no effect, we still inserted them for pragmatical code compatibilities. 

- `s_lim=nothing` : parameter used in order to avoid the divergence of the ``\\mathcal{R}`` and 
  ``\\mathfrak{R}`` denominators: when ``0 \\leq s \\leq s_\\mathrm{lim}`` the returned values are
  ```math
  \\forall \\, s \\in [ 0, s_\\mathrm{lim} ] \\; : \\quad 
      \\mathfrak{R}(s) = 1 - \\frac{1}{\\mathcal{H}_0 \\, s_\\mathrm{lim}} \\; , \\quad
      \\mathcal{R}(s) = 5 s_{\\mathrm{b}} + 
          \\frac{2 - 5 s_{\\mathrm{b}}}{\\mathcal{H}_0 \\, s_\\mathrm{lim}} +  
          \\frac{\\dot{\\mathcal{H}}}{\\mathcal{H}_0^2} - \\mathit{f}_{\\mathrm{evo}} \\; .
  ```
  If `nothing`, the default value stored in `cosmo` will be considered.

- `obs::Union{Bool,Symbol} = :noobsvel` : do you want to consider the observer terms in the computation of the 
  chosen GNC TPCF effect?
  - `:yes` or `true` -> all the observer terms will be considered;
  - `:no` or `false` -> no observer term will be taken into account;
  - `:noobsvel` -> the observer terms related to the observer velocity (that you can find in the CF concerning Doppler)
    will be neglected, the other ones will be taken into account.

See also: [`Point`](@ref), [`Cosmology`](@ref), [`ξ_GNC_multipole`](@ref), 
[`map_ξ_GNC_multipole`](@ref), [`print_map_ξ_GNC_multipole`](@ref)
"""
ξ_GNC_Newtonian_Doppler

##########################################################################################92

##########################################################################################92

##########################################################################################92


"""
    ξ_GNC_Doppler_Newtonian(s1, s2, y, cosmo::Cosmology; 
        b1=nothing, b2=nothing, s_b1=nothing, s_b2=nothing, 
        𝑓_evo1=nothing, 𝑓_evo2=nothing, s_lim=nothing, 
        obs::Union{Bool, Symbol} = :noobsvel
        ) ::Float64

Return the Two-Point Correlation Function (TPCF) given by the cross correlation between the 
Doppler and the Newtonian effects arising from the Galaxy Number Counts (GNC).

It's computed through the symmetric function `ξ_GNC_Newtonian_Doppler`; check its documentation for
more details about the analytical expression and the keyword arguments.
We remember that all the distances are measured in ``h_0^{-1}\\mathrm{Mpc}``.


## Inputs

- `s1` and `s2`: comoving distances where the TPCF has to be calculated;
  
- `y`: the cosine of the angle between the two points `P1` and `P2` wrt the observer;

- `cosmo::Cosmology`: cosmology to be used in this computation; it contains all the splines
  used for the conversion `s` -> `Point`, and all the cosmological parameters ``b``, ...

## Keyword Arguments

- `kwargs...` : Keyword arguments to be passed to the symmetric TPCF.

See also: [`Point`](@ref), [`Cosmology`](@ref), [`ξ_GNC_multipole`](@ref), 
[`map_ξ_GNC_multipole`](@ref), [`print_map_ξ_GNC_multipole`](@ref),
[`ξ_GNC_Newtonian_Doppler`](@ref)
"""
function ξ_GNC_Doppler_Newtonian(s1, s2, y, cosmo::Cosmology; 
    b1=nothing, b2=nothing, s_b1=nothing, s_b2=nothing, 𝑓_evo1=nothing, 𝑓_evo2=nothing,
    s_lim=nothing, kwargs...)
    
    b1 = isnothing(b1) ? cosmo.params.b1 : b1
    b2 = isnothing(b2) ? cosmo.params.b2 : b2
    s_b1 = isnothing(s_b1) ? cosmo.params.s_b1 : s_b1
    s_b2 = isnothing(s_b2) ? cosmo.params.s_b2 : s_b2
    𝑓_evo1 = isnothing(𝑓_evo1) ? cosmo.params.𝑓_evo1 : 𝑓_evo1
    𝑓_evo2 = isnothing(𝑓_evo2) ? cosmo.params.𝑓_evo2 : 𝑓_evo2

    ξ_GNC_Newtonian_Doppler(s2, s1, y, cosmo; b1=b2, b2=b1, s_b1=s_b2, s_b2=s_b1,
        𝑓_evo1=𝑓_evo2, 𝑓_evo2=𝑓_evo1, s_lim=s_lim, kwargs...)
end

