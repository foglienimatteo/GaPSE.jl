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
# along with  If not, see <http://www.gnu.org/licenses/>.
#



"""
     TF(
          left_value::Float64
          left::Float64

          spline::Dierckx.Spline1D

          r_si::Float64
          r_b::Float64
          r_a::Float64
          right::Float64
          )

Contains all the information useful in order to return the Transfer Function value from:
- a spline inside the interval `left ≤ x ≤ right`
- the associated power law for `x > right` (with "right" coefficients `r_si`, `r_b` and `r_a`)
- the associated constant left value `left_value` for `x < left`



## Arguments 

- `left_value::Float64` : the constant value that must be returned in case `x < left`.

- `left::Float64` : the break between the left power-law (for `x <left`) and the 
  spline (for `x ≥ left`); its value is the `xs[begin]` one.

- `spline::Dierckx.Spline1D` : spline that interpolates between the real values of the 
  integral calculated inside the range `left ≤ x ≤ right`

- `right::Float64` : the break between the right power-law (for `x ≥ left`) and the 
  spline (for `x ≤ right`); its value is the `xs[end]` one.

- `r_si, r_b, r_a :: Float64` : coefficient for the spurious power-law 
  ``y = f(x) = a + b \\, x^s`` for the RIGHT edge; when an input value `x > right` is
  given, the returned one is obtained from `power_law` with this coefficients (
  where, of course, `r_si` is the exponent, `r_b`` the coefficient and `r_a` the 
  spurious adding constant). 
  NOTE: for numerical issues, the "pure" power-law ``y = f(x) = b + x^s`` should be used. 

## Constructors

`TF(ks, Tks)` : from a set of `(ks, Tks)` pairs, take the mean of the first 10 as left
value and fit with `power_law_from_data` the last 15.


See also: [`power_law_from_data`](@ref)
"""
struct TF
     left_value::Float64
     left::Float64

     spline::Dierckx.Spline1D

     r_si::Float64
     r_b::Float64
     r_a::Float64
     right::Float64


     function TF(ks, Tks)
          left_val = sum(Tks[1:10]) / 10
          r_si, r_b, r_a = power_law_from_data(
               ks, Tks, [-2.0, 1.0], ks[end-15], ks[end]; con=false)
          spline = Spline1D(ks, Tks; bc="error")

          new(left_val, ks[2], spline, r_si, r_b, r_a, ks[end])
     end

end;


"""
     (f::TF)(x)

Return the value of the `f::TF` as follows:
```math
f(x)=
\\begin{cases}
y_{\\mathrm left} \\; ,
    \\quad \\quad x < \\mathrm{left}\\\\
\\mathrm{spline}(x) \\; , \\quad \\mathrm{left} \\leq x \\leq \\mathrm{right} \\\\
a_\\mathrm{R} + b_\\mathrm{R} \\, x ^ {s_\\mathrm{R}} \\; , 
\\quad x > \\mathrm{right}
\\end{cases}
```

where ``y_\\mathrm{left}``,``\\mathrm{left}``,
``\\mathrm{spline}``, ``a_\\mathrm{R}``, ``b_\\mathrm{R}``, ``s_\\mathrm{R}`` and 
``\\mathrm{right}`` are all stored inside the `TF` considered.

See also: [`TF`](@ref)
"""
function (tf::TF)(k)
     if k < tf.left
          return tf.left_value
     elseif k > tf.right
          return power_law(k, tf.r_si, tf.r_b, tf.r_a)
     else
          return tf.spline(k)
     end
end;


##########################################################################################92


"""
     α_bias(k, tf::TF; bf=1.0, D=1.0, Ω_M0=0.29992)

Return the coefficient ``\\alpha_{\\rm bias}`` that relates the Non-Gaussian density fluctiations 
``\\delta_{\\rm NG}`` and the Non-Gaussian gravitational potential ``\\Phi_{\\rm NG}``
in Fourier space:

```math
\\delta_{\\rm NG}(k) = \\alpha(k, z) \\,  \\Phi_{\\rm NG}(k) \\; \\;  , \\quad
\\alpha(k, z) = \\frac{2}{3} \\frac{k^2 T_m(k) D(z)}{\\Omega_{\\mathrm{M}0}} \\left(\\frac{c}{H_0}\\right)^2 
\\; \\; ,  \\quad \\; \\alpha_{\\rm bias} = \\frac{b_{\\phi} f_{\\rm NL}}{\\alpha(k, z)} \\; .
```

- `bf=1.0` : value of the degenerate product ``b_{\\phi} f_{\\rm NL}``.
- `D = 1.0` : value of the linear growth factor ``D`` at present day; inside this function, 
  it is multiplied for a constant ``q`` in order to get ``q \\, D = ``  


See also: [`TF`](@ref)
"""
function α_bias(k, tf::TF; bf=1.0, D=1.0, Ω_M0=0.29992, q=0.779017)
     return 1.5 * bf * Ω_M0 * (100 / 299792.458)^2 / (0.779017 * D * k^2 * tf(k))
end



"""
     IntegralIPSalpha(
          l_si::Float64
          l_b::Float64
          l_a::Float64
          left::Float64

          spline::Dierckx.Spline1D

          r_si::Float64
          r_b::Float64
          r_a::Float64
          right::Float64
     )

Contains all the information useful in order to return the value of the integral
of the Input Power Spectrum weighted with the `α_bias` function. In other words,
return this expression:

```math
\\int_0^\\infty \\frac{\\mathrm{d} q}{2 \\pi^2} \\, q^2 \\,
    \\frac{j_\\ell(qs)}{(qs)^n} \\, P(q) \\, \\alpha_{\\mathrm{bias}} \\; ,
```

where ``P(q)`` is the Input Power Spectrum and

```math
\\delta_{\\rm NG}(k) = \\alpha(k, z) \\,  \\Phi_{\\rm NG}(k) \\; \\;  , \\quad
\\alpha(k, z) = \\frac{2}{3} \\frac{k^2 T_m(k) D(z)}{\\Omega_{\\mathrm{M}0}} \\left(\\frac{c}{H_0}\\right)^2 
\\; \\; ,  \\quad \\; \\alpha_{\\rm bias} = \\frac{b_{\\phi} f_{\\rm NL}}{\\alpha(k, z)} \\; .
```

## Arguments 

- `l_si, l_b, l_a ::Float64` : coefficient for the spurious power-law 
  ``y = f(x) = a + b \\, x^s`` for the LEFT edge; when an input value `x < left` is
  given, the returned one is obtained from `power_law` with this coefficients (
  where, of course, `l_si` is the exponent, `l_b` the coefficient and `l_a` the 
  spurious adding constant). 

- `left::Float64` : the break between the left power-law (for `x < left`) and the 
  spline (for `x ≥ left`); its value is the `fit_min` of the used constructor.

- `spline::Dierckx.Spline1D` : spline that interpolates between the real values of the 
  integral calculated inside the range `left ≤ x ≤ right`

- `right::Float64` : the break between the right power-law (for `x > right`) and the 
  spline (for `x ≤ right`); its value is the `fit_max` of the used constructor.

- `r_si, r_b, r_a ::Float64` : coefficient for the spurious power-law 
  ``y = f(x) = a + b \\, x^s`` for the RIGHT edge; when an input value `x > right` is
  given, the returned one is obtained from `power_law` with this coefficients (
  where, of course, `r_si` is the exponent, `r_b` the coefficient and `r_a` the 
  spurious adding constant). 
  NOTE: for numerical issues, only the "pure" power-law ``y = f(x) = b + x^s`` can be used. 
  In other words, it always set `r_a = 0.0`.

## Constructors

    IntegralIPSalpha(tf::TF, cosmo::Cosmology, l, n=0; D=nothing, bf=1.0,
        N::Int=1024, kmin=1e-6, kmax=1e4, s0=1e-4,
        fit_left_min=nothing, fit_left_max=nothing, p0_left=nothing,
        fit_right_min=nothing, fit_right_max=nothing, p0_right=nothing)

The integral obtained with this constructor is calculated through `xicalc`, and
expanded with power-laws at the edges.

- `tf::TF`: the struct that contains all the data concerning the Transfer Function.

- `cosmo::Cosmology` : cosmology to be used in this computation

- `l` : degree of the spherical Bessel function to be used.

- `n=0` : degree of the exponent for the denominator. The interesting case is
  only the default value ``0``.

- `D = nothing` : value of the linear growth factor ``D`` to be used. If `nothing`,
  it will be internally set as ``D(z_{\\mathrm{eff}})``, where ``z_{\\mathrm{eff}}`` is
  the effective redshift for the input cosmology.

- `bf = 1.0` : value of the degenerate product ``b_{\\phi} f_{\\rm NL}``.

- `kmin = 1e-6, kmax = 1e4, s0 = 1e-4` : values to be passed to `xicalc` for the
  integration

- `fit_left_min = 2.0, fit_left_max = 10.0` : the limits (min and max) where the integral
  must be fitted with a power law, for small distances. This operation is necessary, because `xicalc`,
  in this context, gives wrong results for too small input distance `s`; nevertheless,
  this integral has fixed power-law trends for ``s \\rightarrow 0``, so this approach gives
  good results.

- `p0_left = nothing` : vector with the initial values for the left power-law fitting; its length must
  be 2 (if you want to fit with a pure power-law ``y = f(x) = b x^s``, so only `l_si` and `l_b` 
  are matter of concern) or 3 (if you want to fit with a spurious power-law ``y = f(x) = a + b x^s``,
  so you are also interested in `l_a`); if `nothing`, it will be
  automatically set `p0 = [-1.0, 1.0]`.

- `fit_right_min = nothing, fit_right_max = nothing` : the limits (min and max) where the integral
  must be fitted with a power law, for high distances. 
  This integral has fixed power-law trends for ``s \\rightarrow \\infty``, so this approach gives
  good results. If `nothing`, the last 15 points returned from `xicalc` are used for
  this fitting.
  NOTE: for numerical issues, only the "pure" power-law ``y = f(x) = b + x^s`` can be used. 

- `p0_right = nothing` : vector with the initial values for the left power-law fitting; its length must
  be 2 (to fit with a pure power-law ``y = f(x) = b x^s``, so only `r_si` and `r_b` 
  are matter of concern); if `nothing`, it will be
  automatically set `p0 = [-4.0, 1.0]`.

All the power-law fitting (both "pure" and spurious) are made through the 
local function `power_law_from_data`.

See also: [`power_law_from_data`](@ref), [`power_law`](@ref), 
[`Cosmology`](@ref), [`α_bias`](@ref)
"""
struct IntegralIPSalpha
     l_si::Float64
     l_b::Float64
     l_a::Float64
     left::Float64

     spline::Dierckx.Spline1D

     r_si::Float64
     r_b::Float64
     r_a::Float64
     right::Float64

     function IntegralIPSalpha(tf::TF, cosmo::Cosmology, l, n=0;
          D=nothing, bf=1.0,
          N::Int=1024, kmin=1e-6, kmax=1e4, s0=1e-4,
          fit_left_min=nothing, fit_left_max=nothing, p0_left=nothing,
          fit_right_min=nothing, fit_right_max=nothing, p0_right=nothing)

          DD = isnothing(D) ? cosmo.D_of_s(cosmo.s_eff) : D
          Ω_M00 = cosmo.params.Ω_M0

          rs, xis = xicalc(k -> cosmo.IPS(k) * α_bias(k, tf; bf=bf, D=DD, Ω_M0=Ω_M00), l, n;
               N=N, kmin=kmin, kmax=kmax, r0=s0)

          fit_left_MIN = !isnothing(fit_left_min) ? fit_left_min : begin
               l ≈ 0.0 ? 5e-2 : l ≈ 2.0 ? 5e-1 : rs[2]
          end
          fit_left_MAX = !isnothing(fit_left_max) ? fit_left_max : begin
               l ≈ 0.0 ? 1e-1 : l ≈ 2.0 ? 1e0 : rs[16]
          end
          p_0_left = isnothing(p0_left) ? [-1.0, 1.0] : p0_left
          l_si, l_b, l_a = power_law_from_data(
               rs, xis, p_0_left, fit_left_MIN, fit_left_MAX; con=false)

          fit_right_MIN = isnothing(fit_right_min) ? rs[length(rs)-16] : fit_right_min
          fit_right_MAX = isnothing(fit_right_max) ? rs[length(rs)-1] : fit_right_max
          p_0_right = isnothing(p0_right) ? [-4.0, 1.0] : p0_right
          r_si, r_b, r_a = power_law_from_data(
               rs, xis, p_0_right, fit_right_MIN, fit_right_MAX; con=false)

          ind_left = findfirst(x -> x > fit_left_MIN, rs) - 1
          ind_right = findfirst(x -> x >= fit_right_MAX, rs)
          new_rs = vcat(rs[ind_left:ind_right])
          new_Js = vcat(xis[ind_left:ind_right])
          spline = Spline1D(new_rs, new_Js; bc="error")

          #println("\nleft = $l_si , $l_b , $l_a, $fit_left_min")
          #println("right = $r_si , $r_b , $r_a, $fit_right_MAX\n")

          new(l_si, l_b, l_a, fit_left_MIN, spline, r_si, r_b, r_a, fit_right_MAX)
     end
end;


"""
     (f::IntegralIPSalpha)(x)

Return the value of the `f::IntegralIPSalpha` as follows:
```math
f(x)=
\\begin{cases}
a_\\mathrm{L} + b_\\mathrm{L} \\, x ^ {s_\\mathrm{L}} \\; ,
    \\quad x < \\mathrm{left}\\\\
\\mathrm{spline}(x) \\; , \\quad \\mathrm{left} \\leq x \\leq \\mathrm{right} \\\\
a_\\mathrm{R} + b_\\mathrm{R} \\, x ^ {s_\\mathrm{R}} \\; , 
\\quad x > \\mathrm{right}
\\end{cases}
```

where ``a_\\mathrm{L}``, ``b_\\mathrm{L}``, ``s_\\mathrm{L}``, ``\\mathrm{left}``,
``\\mathrm{spline}``, ``a_\\mathrm{R}``, ``b_\\mathrm{R}``, ``s_\\mathrm{R}`` and 
``\\mathrm{right}`` are all stored inside the `IntegralIPSalpha` considered.

See also: [`IntegralIPSalpha`](@ref)
"""
function (Jl::IntegralIPSalpha)(x)
     if x < Jl.left
          return power_law(x, Jl.l_si, Jl.l_b, Jl.l_a)
     elseif x > Jl.right
          #warning("i am going too right! ")
          return power_law(x, Jl.r_si, Jl.r_b, Jl.r_a)
     else
          return Jl.spline(x)
     end
end;



##########################################################################################92



struct CosmoPNG
     tf::TF
     file_TF::String

     J0::IntegralIPSalpha
     J2::IntegralIPSalpha

     params::Dict{Symbol,T1} where {T1}

     function CosmoPNG(
          cosmo::Cosmology,
          file_TF::String;
          comments::Bool=true, D=nothing, bf=1.0,
          flm0=5e-2, flM0=1e-1, flm2=5e-1, flM2=1e0,
          kmin=1e-6, kmax=1e4, N::Int=1024
     )

          table = readdlm(file_TF; comments=comments)
          ks_tf = convert(Vector{Float64}, table[:, 1])
          pks_tf = convert(Vector{Float64}, table[:, 2])

          DD = isnothing(D) ? cosmo.D_of_s(cosmo.s_eff) : D
          tf = TF(ks_tf, pks_tf)

          J0 = IntegralIPSalpha(tf, cosmo, 0, 0; D=DD, bf=bf, kmin=kmin, kmax=kmax,
               fit_left_min=flm0, fit_left_max=flM0, N=N)
          J2 = IntegralIPSalpha(tf, cosmo, 2, 0; D=DD, bf=bf, kmin=kmin, kmax=kmax,
               fit_left_min=flm2, fit_left_max=flM2, N=N)

          params = Dict(:D => D, :bf => bf, :flm0 => flm0, :flM0 => flM0,
               :flm2 => flm2, :flM2 => flM2, :kmin => kmin, :kmax => kmax, :N => N)

          new(tf, file_TF, J0, J2, params)
     end
end


##########################################################################################92



function ξ_S_L0(P::Point, cosmo::Cosmology, cosmopng::CosmoPNG)
     b = cosmo.params.b
     s = P.comdist

     Peff = Point(cosmo.s_eff, cosmo)
     D, f = Peff.D, Peff.f

     2.0 * (b + f / 3.0) * D^2 * cosmopng.J0(s)
end

function ξ_S_L0(s1, cosmo::Cosmology, cosmopng::CosmoPNG)
     P1 = Point(s1, cosmo)
     return ξ_S_L0(P1, cosmo, cosmopng)
end

function ξ_S_L2(P::Point, cosmo::Cosmology, cosmopng::CosmoPNG)
     s = P.comdist

     Peff = Point(cosmo.s_eff, cosmo)
     D, f = Peff.D, Peff.f

     -4.0 / 3.0 * f * D^2 * cosmopng.J2(s)
end

function ξ_S_L2(s1, cosmo::Cosmology, cosmopng::CosmoPNG)
     P1 = Point(s1, cosmo)
     return ξ_S_L2(P1, cosmo, cosmopng)
end


function ξ_S(s1, y, cosmo::Cosmology, cosmopng::CosmoPNG)
     ξ_S_L0(s1, cosmo, cosmopng) + ξ_S_L2(s1, cosmo, cosmopng) * Pl(y, 2)
end



##########################################################################################92



function integrand_ξ_S_multipole(s, μ, cosmo::Cosmology, cosmopng::CosmoPNG;
     L::Int=0, use_windows::Bool=true)

     res = if use_windows == true
          ξ_S(s, μ, cosmo, cosmopng) .* (
               spline_integrF(s, μ, cosmo.windowFint) / cosmo.WFI_norm * Pl(μ, L)
          )
     else
          ξ_S(s, μ, cosmo, cosmopng) .* Pl(μ, L)
     end

     return (2.0 * L + 1.0) / 2.0 * res
end

function ξ_S_multipole(
     s, cosmo::Cosmology, cosmopng::CosmoPNG;
     L::Int=0,
     use_windows::Bool=true,
     enhancer::Float64=1e6,
     μ_atol::Float64=0.0,
     μ_rtol::Float64=1e-2)

     orig_f(μ) = enhancer * integrand_ξ_S_multipole(s, μ, cosmo, cosmopng;
          L=L, use_windows=use_windows)

     int = quadgk(μ -> orig_f(μ), -1.0, 1.0; atol=μ_atol, rtol=μ_rtol)[1]

     return int / enhancer
end

function map_ξ_S_multipole(cosmo::Cosmology, cosmopng::CosmoPNG,
     v_ss=nothing;
     pr::Bool=true,
     N_log::Int=1000,
     L::Int=0,
     kwargs...)

     t1 = time()
     ss = isnothing(v_ss) ? 10 .^ range(0, 3, length=N_log) : v_ss
     xis = pr ? begin
          @showprogress "ξ_S, L=$L: " [
               ξ_S_multipole(s, cosmo, cosmopng; L=L, kwargs...) for s in ss
          ]
     end : [
          ξ_S_multipole(s, cosmo, cosmopng; L=L, kwargs...) for s in ss
     ]

     t2 = time()
     pr && println("\ntime needed for map_ξ_S_multipole " *
                   "[in s] = $(@sprintf("%.5f", t2-t1)) ")
     return (ss, xis)
end


function print_map_ξ_S_multipole(
     cosmo::Cosmology, cosmopng::CosmoPNG,
     out::String,
     v_ss=nothing;
     kwargs...)

     check_parent_directory(out)
     check_namefile(out)

     t1 = time()
     vec = map_ξ_S_multipole(cosmo, cosmopng, v_ss; kwargs...)
     t2 = time()

     isfile(out) && run(`rm $out`)
     open(out, "w") do io
          println(io, BRAND)

          println(io, "#\n# This is an integration map on mu of the ξ multipole of S,\n" *
                      "# which is the difference between the Power Spectrum with PNG and with f_NL = 0.")
          println(io, "# Transfer function read from the file: $(cosmopng.file_TF)")
          print(io, "# Parameters used for the considered CosmoPNG: ")
          print(io, "\n")
          for key in keys(cosmopng.params)
               println(io, "# \t\t$(key) = $(cosmopng.params[key])")
          end
          println(io, "#")

          parameters_used(io, cosmo; logo=false)
          println(io, "# computational time needed (in s) : $(@sprintf("%.4f", t2-t1))")
          print(io, "# kwards passed: ")

          if isempty(kwargs)
               println(io, "none")
          else
               print(io, "\n")
               for key in keys(kwargs)
                    println(io, "# \t\t$(key) = $(kwargs[key])")
               end
          end

          println(io, "# ")
          println(io, "# s [Mpc/h_0] \t \t xi")
          for (s, xi) in zip(vec[1], vec[2])
               println(io, "$s \t $xi")
          end
     end
end


