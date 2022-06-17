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

function α_bias(k, cosmo::GaPSE.Cosmology)
     D = cosmo.D_of_s(cosmo.s_eff)
     Ω_M0 = cosmo.params.Ω_M0
     #1.68647
     #return 1.0
     return 3.0 / 2.0 * Ω_M0 * (100 / 299792.458)^2 / (0.779017 * D * k^2 * TF_casto(k))
end;


struct my_TF
     left_value::Float64
     left::Float64

     spline::Dierckx.Spline1D

     r_si::Float64
     r_b::Float64
     r_a::Float64
     right::Float64


     function my_TF(ks, Tks)
          left_val = sum(Tks[1:10]) / 10
          r_si, r_b, r_a = GaPSE.power_law_from_data(
               ks, Tks, [-2.0, 1.0], ks[end-15], ks[end]; con=false)
          spline = Spline1D(ks, Tks; bc="error")

          new(left_val, ks[2], spline, r_si, r_b, r_a, ks[end])
     end

end;


function (TF::my_TF)(k)
     if k < TF.left
          return TF.left_value
     elseif k > TF.right
          return GaPSE.power_law(k, TF.r_si, TF.r_b, TF.r_a)
     else
          return TF.spline(k)
     end
end;


struct f_NL_IntegralIPS
     l_si::Float64
     l_b::Float64
     l_a::Float64
     left::Float64

     spline::Dierckx.Spline1D

     r_si::Float64
     r_b::Float64
     r_a::Float64
     right::Float64

     function f_NL_IntegralIPS(cosmo::GaPSE.Cosmology, l, n;
          N::Int=1024, kmin=1e-4, kmax=1e3, s0=1e-3,
          fit_left_min=nothing, fit_left_max=nothing, p0_left=nothing, con=false,
          fit_right_min=nothing, fit_right_max=nothing, p0_right=nothing)


          rs, xis = xicalc(k -> cosmo.IPS(k) * α_bias(k, cosmo), l, n;
               N=N, kmin=kmin, kmax=kmax, r0=s0)

          fit_left_MIN = isnothing(fit_left_min) ? rs[2] : fit_left_min
          fit_left_MAX = isnothing(fit_left_max) ? rs[16] : fit_left_max
          p_0_left = isnothing(p0_left) ? (con == true ? [-1.0, 1.0, 0.0] : [-1.0, 1.0]) : p0_left
          l_si, l_b, l_a = GaPSE.power_law_from_data(
               rs, xis, p_0_left, fit_left_MIN, fit_left_MAX; con=con)

          fit_right_MIN = isnothing(fit_right_min) ? rs[length(rs)-16] : fit_right_min
          fit_right_MAX = isnothing(fit_right_max) ? rs[length(rs)-1] : fit_right_max
          p_0_right = isnothing(p0_right) ? [-4.0, 1.0] : p0_right
          r_si, r_b, r_a = GaPSE.power_law_from_data(
               rs, xis, p_0_right, fit_right_MIN, fit_right_MAX; con=false)

          ind_left = findfirst(x -> x > fit_left_MIN, rs) - 1
          ind_right = findfirst(x -> x >= fit_right_MAX, rs)
          new_rs = vcat(rs[ind_left:ind_right])
          new_Is = vcat(xis[ind_left:ind_right])
          spline = Spline1D(new_rs, new_Is; bc="error")

          #println("\nleft = $l_si , $l_b , $l_a, $fit_left_min")
          #println("right = $r_si , $r_b , $r_a, $fit_right_MAX\n")

          new(l_si, l_b, l_a, fit_left_MIN, spline, r_si, r_b, r_a, fit_right_MAX)
     end
end;


function (Iln::f_NL_IntegralIPS)(x)
     if x < Iln.left
          return GaPSE.power_law(x, Iln.l_si, Iln.l_b, Iln.l_a)
     elseif x > Iln.right
          #warning("i am going too right! ")
          return GaPSE.power_law(x, Iln.r_si, Iln.r_b, Iln.r_a)
     else
          return Iln.spline(x)
     end
end;



##########################################################################################92



function ξ_S_L0(P::GaPSE.Point, cosmo::GaPSE.Cosmology)
     b = cosmo.params.b
     s = P.comdist

     Peff = GaPSE.Point(cosmo.s_eff, cosmo)
     D, f = Peff.D, Peff.f

     2.0 * (b + f / 3.0) * D^2 * I00(s)
end

function ξ_S_L0(s1, cosmo::GaPSE.Cosmology)
     P1 = GaPSE.Point(s1, cosmo)
     return ξ_S_L0(P1, cosmo)
end

function ξ_S_L2(P::GaPSE.Point, cosmo::GaPSE.Cosmology)
     s = P.comdist

     Peff = GaPSE.Point(cosmo.s_eff, cosmo)
     D, f = Peff.D, Peff.f

     -4.0 / 3.0 * f * D^2 * I20(s)
end

function ξ_S_L2(s1, cosmo::GaPSE.Cosmology)
     P1 = GaPSE.Point(s1, cosmo)
     return ξ_S_L2(P1, cosmo)
end

function ξ_S(s1, y, cosmo::GaPSE.Cosmology)
     return ξ_S_L0(s1, cosmo) + ξ_S_L2(s1, cosmo) * Pl(y, 2)
end;



##########################################################################################92



function integrand_ξ_S_multipole(s, μ, cosmo::GaPSE.Cosmology;
     L::Int=0, use_windows::Bool=true)

     res = if use_windows == true
          ξ_S(s, μ, cosmo) .* (GaPSE.spline_integrF(s, μ, cosmo.windowFint) / cosmo.WFI_norm * Pl(μ, L))
     else
          ξ_S(s, μ, cosmo) .* Pl(μ, L)
     end

     return (2.0 * L + 1.0) / 2.0 * res
end

function ξ_S_multipole(
     s, cosmo::GaPSE.Cosmology;
     L::Int=0,
     use_windows::Bool=true,
     enhancer::Float64=1e6,
     μ_atol::Float64=0.0,
     μ_rtol::Float64=1e-2)

     orig_f(μ) = enhancer * integrand_ξ_S_multipole(s, μ, cosmo;
          L=L, use_windows=use_windows)

     int = quadgk(μ -> orig_f(μ), -1.0, 1.0; atol=μ_atol, rtol=μ_rtol)[1]

     return int / enhancer
end

function map_ξ_S_multipole(cosmo::GaPSE.Cosmology,
     v_ss=nothing;
     pr::Bool=true,
     N_log::Int=1000,
     L::Int=0,
     kwargs...)

     t1 = time()
     ss = isnothing(v_ss) ? 10 .^ range(0, 3, length=N_log) : v_ss
     xis = pr ? begin
          @showprogress "ξ_S, L=$L: " [
               ξ_S_multipole(s, cosmo; L=L, kwargs...) for s in ss
          ]
     end : [
          ξ_S_multipole(s, cosmo; L=L, kwargs...) for s in ss
     ]

     t2 = time()
     pr && println("\ntime needed for map_ξ_S_multipole " *
                   "[in s] = $(@sprintf("%.5f", t2-t1)) ")
     return (ss, xis)
end;