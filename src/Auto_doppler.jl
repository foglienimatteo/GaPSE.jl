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

function ξ_doppler(P1::Point, P2::Point, y; enhancer = 1.0)
     s1, D1, f1, ℋ1, ℛ1 = P1.comdist, P1.D, P1.f, P1.ℋ, P1.ℛ
     s2, D2, f2, ℋ2, ℛ2 = P2.comdist, P2.D, P2.f, P2.ℋ, P2.ℛ

     delta_s = s(P1.comdist, P2.comdist, y)
     prefac = D1 * D2 * f1 * f2 * ℛ1 * ℛ2 * ℋ1 * ℋ2
     c1 = 3 * s1 * s2 - 2 * y * (s1^2 + s2^2) + s1 * s2 * y^2

     parenth = I00(delta_s) / 45.0 + I20(delta_s) / 31.5 + I40(delta_s) / 105.0

     first = prefac * (c1 * parenth + I02(delta_s) * y * delta_s^2 / 3.0)

     return enhancer * first
end


##########################################################################################92


function integrand_on_mu_doppler(s1, s, μ,
     BS::BackgroundSplines; L::Integer = 0, enhancer = 1.0,
     use_windows::Bool = true)

     s2_value = s2(s1, s, μ)
     y_value = y(s1, s, μ)

     if use_windows == true
          ϕ_s2 = ϕ(s2_value)
          (ϕ_s2 > 0.0) || (return 0.0)
          P1, P2 = Point(s1, BS), Point(s2_value, BS)
          val = ξ_doppler(P1, P2, y_value, enhancer = enhancer)
          return val * ϕ_s2 * spline_F(s / s1, μ) * Pl(μ, L)
     else
          P1, P2 = Point(s1, BS), Point(s2_value, BS)
          val = ξ_doppler(P1, P2, y_value, enhancer = enhancer)
          return val * Pl(μ, L)
     end
end


function integrand_on_mu_doppler(P1::Point, P2::Point, y;
     L::Integer = 0, enhancer = 1.0, use_windows::Bool = true)

     if use_windows == true
          ϕ_s2 = ϕ(P2.comdist)
          (ϕ_s2 > 0.0) || (return 0.0)
          val = ξ_doppler(P1, P2, y, enhancer = enhancer)
          return val * ϕ_s2 * spline_F(s / s1, μ) * Pl(μ, L)
     else
          val = ξ_doppler(P1, P2, y, enhancer = enhancer)
          return val * Pl(μ, L)
     end
end

#=
@doc raw"""
     ξ_doppler(s1, s2, y; enhancer = 1, tol = 1)

Return the Doppler auto-correlation function, defined as follows:
```math
\xi^{v_{\parallel}v_{\parallel}} (s_1, s_2, \cos{\theta}) 
= D_1 D_2 f_1 f_2 \mathcal{H}_1 \mathcal{H}_2 \mathcal{R}_1 \mathcal{R}_2 
(J_{00}\, I^0_0(s) + J_{02}\,I^0_2(s) + J_{04}\,I^0_4(s) + J_{20}\,I^2_0(s))
```
where ``D_1 = D(s_1)``, ``D_2 = D(s_2)`` and so on, ``\mathcal{H} = a H`` and 
the J coefficients are given by (with ``y = \cos{\theta}``):
```math
\begin{align*}
    J_{00} (s_1, s_2, y) & = \frac{1}{45} (y^2 s_1 s_2 - 2y(s_1^2 + s_2^2) + 3s_1 s_2) \\
    J_{02} (s_1, s_2, y) & = \frac{2}{63} (y^2 s_1 s_2 - 2y(s_1^2 + s_2^2) + 3s_1 s_2) \\
    J_{04} (s_1, s_2, y) & = \frac{1}{105} (y^2 s_1 s_2 - 2y(s_1^2 + s_2^2) + 3s_1 s_2) \\
    J_{20} (s_1, s_2, y) & = \frac{1}{3} y s^2
\end{align*}
```

`enhancer` is just a float number used in order to deal better with small numbers; the returned
value is actually `enhancer * xi_doppler`, where `xi_doppler` is the true value calculated
as shown before.

`tol` is a flaot number that set the minimum distance of interest between the two galaxies
placed at ``\mathbf{s_1}`` and ``\mathbf{s_2}``; the output is set to `0.0` if it is true
that:
```math
begin{split}
     s = s(s_1, s_2, y) &= \sqrt{s_1^2 + s_2^2 - 2 \, s_1 \, s_2 \, y } \ge \mathrm{tol} \\

     &\;\mathbf{where} \; \, y=\cos\theta = \hat{\mathbf{s}_1}\dot\hat{\mathbf{s}_2}
\end{split}
```
"""
function ξ_doppler(s1, s2, y; enhancer = 1.0)
     delta_s = s(s1, s2, y)
     D1, D2 = D(s1), D(s2)

     f1, ℋ1 = f(s1), ℋ(s1)
     f2, ℋ2 = f(s2), ℋ(s2)
     #f1, ℋ1, ℋ1_p, s_b1 = f(s1), ℋ(s1), ℋ_p(s1), s_b(s1)
     #f2, ℋ2, ℋ2_p, s_b2 = f(s2), ℋ(s2), ℋ_p(s2), s_b(s2)
     #ℛ1, ℛ2 = ℛ(s1, ℋ1, ℋ1_p, s_b1), ℛ(s2, ℋ2, ℋ2_p, s_b2)
     ℛ1, ℛ2 = ℛ(s1, ℋ1), ℛ(s2, ℋ2)

     prefac = D1 * D2 * f1 * f2 * ℛ1 * ℛ2 * ℋ1 * ℋ2
     c1 = 3 * s1 * s2 - 2 * y * (s1^2 + s2^2) + s1 * s2 * y^2

     parenth = I00(delta_s) / 45.0 + I20(delta_s) / 31.5 + I40(delta_s) / 105.0

     first = prefac * (c1 * parenth + I02(delta_s) * y * delta_s^2 / 3.0)

     return enhancer * first
end
=#

#=
function integrand_on_mu_doppler(s1, s, μ; L::Integer = 0, enhancer = 1.0,
     use_windows::Bool = true)
     if use_windows == true
          ϕ_s2 = ϕ(s2(s1, s, μ))
          (ϕ_s2 > 0.0) || (return 0.0)
          val = ξ_doppler(s1, s2(s1, s, μ), y(s1, s, μ), enhancer = enhancer)
          return val * ϕ_s2 * spline_F(s / s1, μ) * Pl(μ, L)
     else
          val = ξ_doppler(s1, s2(s1, s, μ), y(s1, s, μ), enhancer = enhancer)
          return val * Pl(μ, L)
     end
end
=#

#=
function integral_on_mu_doppler(s1, s; L::Integer = 0, enhancer = 1, tol = 1, kwargs...)
     f(μ) = integrand_on_mu_doppler(s1, s, μ; L = L, enhancer = enhancer, tol = tol)
     return quadgk(f, -1, 1; kwargs...)
end

function map_integral_on_mu_doppler(s1 = s_eff;
     L::Integer = 0, pr::Bool = true, enhancer = 1e6, tol = 1, kwargs...)

     t1 = time()
     ss = 10 .^ range(log10(tol), 3, length = 100)
     f(s) = integral_on_mu_doppler(s1, s; L = L, enhancer = enhancer, tol = tol, kwargs...)
     vec = @showprogress [f(s) ./ enhancer for s in ss]
     xis, xis_err = [x[1] for x in vec], [x[2] for x in vec]
     t2 = time()
     pr && println("\ntime needed for map_integral_on_mu_doppler [in s] = $(t2-t1)\n")
     return (ss, xis, xis_err)
end




function print_map_int_on_mu_doppler(out::String; L::Integer = 0,
     s1 = s_eff, pr::Bool = true, kwargs...)

     t1 = time()
     vec = map_integral_on_mu_doppler(s1; L = L, pr = pr, kwargs...)
     t2 = time()

     isfile(out) && run(`rm $out`)
     open(out, "w") do io
          println(io, "# This is an integration map on mu of xi_doppler.")
          parameters_used(io)
          println(io, "# computational time needed (in s) : $(@sprintf("%.4f", t2-t1))")
          print(io, "# kwards passed: ")

          if isempty(kwargs)
               println(io, "none")
          else
               print(io, "\n")
               for (i, key) in enumerate(keys(kwargs))
                    println(io, "# \t\t$(key) = $(kwargs[key])")
               end
          end

          println(io, "\ns \t xi \t xi_error")
          for (s, xi, xi_err) in zip(vec[1], vec[2], vec[3])
               println(io, "$s \t $xi \t $(xi_err)")
          end
     end
end

=#

##########################################################################################92


#=
# mean time of evaluation: 141 seconds!
# too long to be used
function PS_doppler_exact(; int_s_min = 1e-2, int_s_max = 2 * s_max, N = 128,
     L::Integer = 0, pr::Bool = true, kwargs...)

     if ϕ(s_eff) > 0
          t1 = time()
          ks, pks = xicalc(s -> 2 * π^2 * integral_on_mu_doppler(s_eff, s; L = L, kwargs...)[1], L, 0;
               N = N, kmin = int_s_min, kmax = int_s_max, r0 = 1 / int_s_max)
          t2 = time()
          pr && println("\ntime needed for PS_doppler [in s] = $(t2-t1)\n")

          if iseven(L)
               return ks, ((2 * L + 1) / A(s_min, s_max, θ_MAX) * ϕ(s_eff) * (-1)^(L / 2)) .* pks
          else
               return ks, ((2 * L + 1) / A(s_min, s_max, θ_MAX) * ϕ(s_eff) * (-im)^L) .* pks
          end
     else
          throw(ErrorException(
               "ϕ(s_eff) should be >0, but in s_eff=$s_eff " *
               "is ϕ(s_eff) = $(ϕ(s_eff))! \n"
          ))
     end
end
=#

