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


"""
     derivate_point(xp, yp, x1, y1, x2, y2)

Return the derivative in `(xp, yp)`, given the neighboor points
`(x1,y1)` and `(x2,y2)`, with `x1 < xp < x2`.
It is not assumed that `x2 - xp = xp - x1`.
"""
function derivate_point(xp, yp, x1, y1, x2, y2)
     l2, l1 = (x2 - xp), (xp - x1)
     m2, m1 = (y2 - yp) / l2, (yp - y1) / l1
     res = (m1 * l2 + m2 * l1) / (l1 + l2)
     #println(res)
     return res
end

function derivate_vector(xs, ys; N::Integer = 1)
     if N == 1
          real_vec = [derivate_point(xs[i], ys[i], xs[i-1], ys[i-1], xs[i+1], ys[i+1])
                      for i in (N+1):(length(xs)-N)]
          return vcat(real_vec[begin], real_vec, real_vec[end])
     elseif N > 1
          vec = [derivate_point(xs[i], ys[i], xs[i-j], ys[i-j], xs[i+j], ys[i+j])
                 for i in (N+1):(length(xs)-N), j in 1:N]
          real_vec = [sum(row) / N for row in eachrow(vec)]
          return vcat([real_vec[begin] for i in 1:N], real_vec,
               [real_vec[end] for i in 1:N])
     else
          throw(ErrorException(" N must be an integer >1, not $N!"))
     end
end


function spectral_index(xs, ys; N::Integer = 1, con = false)
     derivs = derivate_vector(xs, ys; N = N)

     if con == false
          return [x * d / y for (x, y, d) in zip(xs, ys, derivs)]
     else
          sec_derivs = derivate_vector(xs, derivs; N = N)
          vec = [x * d2 / d for (x, d, d2) in zip(xs, derivs, sec_derivs)]
          return vec .+ 1.0
     end
end



@doc raw"""
     mean_spectral_index(xs, ys; N::Integer = 1, con = false)

Assuming that the input `ys` follow a power law distribution, 
return the mean spectral index ``\langle S \rangle`` of them.

The spectral index ``S`` of a generic function ``f = f(x)`` is
defined as:
```math
     S = \frac{\partial \log f(x)}{\partial \log x} 
          = \frac{x}{f(x)} \frac{\partial f(x)}{\partial x} 
```
"""
function mean_spectral_index(xs, ys; N::Integer = 1, con = false)
     vec = spectral_index(xs, ys; N = N, con = con)[begin+2*N:end-2*N]
     return sum(vec) / length(vec)
end


power_law(x, si, b, a) = a .+ b .* (x .^ si)

function power_law_b_a(xs, ys, si, p0; con = false)
     if con == true
          @assert length(p0) == 2 " b, a to be fitted, so length(p0) must be 2!"
          fit = curve_fit((x, p) -> power_law(x, si, p[1], p[2]), xs, ys, p0)
          return coef(fit)
     else
          @assert length(p0) == 1 " b to be fitted, so length(p0) must be 1!"
          fit = curve_fit((x, p) -> power_law(x, si, p[1], 0.0), xs, ys, p0)
          return vcat(coef(fit), 0.0)
     end
end

#=
power_law_b(x1, y1, x2, y2, si) = (y2 - y1)/(x2^si - x1^si)
function power_law_b(ixs, ys, sis; logscale=false)
    xs = !logscale ? ixs : begin
            fac = ixs[begin+1] / ixs[begin];
            [x/(fac^i) for (i,x) in enumerate(ixs)]
        end
    bs = [power_law_b(xs[i], ys[i], xs[i+1], ys[i+1], sis[i]) for i in 1:length(xs)-1]
    return vcat(bs, bs[end])
end

power_law_a(x, y, b, si) = y - b*(x^si)

function power_law_a(ixs, ys, bs, sis; logscale=false)
   xs = !logscale ? ixs : begin
            fac = ixs[begin+1] / ixs[begin];
            [x/(fac^i) for (i,x) in enumerate(ixs)]
        end
    [y - b*(x^si) for (x,y,b,si) in zip(xs,ys,bs,sis)]
end
=#

##########################################################################################92


#=
function my_power_law_from_data(xs, ys, p0, x1::Number, x2::Number; N = 3, con = false)
     @assert length(xs) == length(ys) "xs and ys must have same length"
     #Num = length(xs)
     new_xs = xs[x1.<xs.<x2]
     new_ys = ys[x1.<xs.<x2]

     if con == false
          @assert length(p0) == 2 " si,b to be fitted, so length(p0) must be 2!"
          my_si = mean_spectral_index(new_xs, new_ys; N = N, con = con)
          my_b, my_a = power_law_b_a(new_xs, new_ys, my_si, [p0[2]]; con = con)
          return my_si, my_b, my_a
     else
          @assert length(p0) == 3 " si,b,a to be fitted, so length(p0) must be 3!"
          my_si = mean_spectral_index(new_xs, new_ys; N = N, con = con)
          my_b, my_a = power_law_b_a(new_xs, new_ys, my_si, [p0[2], p0[3]]; con = con)
          return my_si, my_b, my_a
     end
end

function my_power_law_from_data(xs, ys, p0; con = false)
     my_power_law_from_data(xs, ys, p0, xs[begin], xs[end]; con = con)
end
=#


function power_law_from_data(xs, ys, p0, x1::Number, x2::Number; con = false)
     @assert length(xs) == length(ys) "xs and ys must have same length"
     new_xs = xs[x1.<xs.<x2]
     mean_exp = sum([log10(abs(y)) for y in ys[x1.<xs.<x2]]) / length(ys[x1.<xs.<x2])
     enhancer = 10.0^(-mean_exp)
     new_ys = ys[x1.<xs.<x2] .* enhancer

     #si = mean_spectral_index(xs, ys; N=N, con=con)
     si, b, a =
          if con == false
               @assert length(p0) == 2 " si,b to be fitted, so length(p0) must be 2!"
               vec = coef(curve_fit((x, p) -> power_law(x, p[1], p[2], 0.0),
                    new_xs, new_ys, p0))
               vcat(vec, 0.0)
          else
               @assert length(p0) == 3 " si,b,a to be fitted, so length(p0) must be 3!"
               coef(curve_fit((x, p) -> power_law(x, p[1], p[2], p[3]),
                    new_xs, new_ys, p0))
          end

     return si, b / enhancer, a / enhancer
end



function power_law_from_data(xs, ys, p0; con = false)
     power_law_from_data(xs, ys, p0, xs[begin], xs[end]; con = con)
end

#=
function power_law_from_data(xs, ys, x1::Number, x2::Number; N=3, 
        con=false, logscale=false)
    @assert length(xs) == length(ys) "xs and ys must have same length"
    Num = length(xs)
    new_xs = xs[x1 .< xs .< x2]
    new_ys = ys[x1 .< xs .< x2]

    sis = spectral_index(new_xs, new_ys; N = N, con = con, logscale = logscale)

    bs = power_law_b(new_xs, new_ys, sis, logscale = logscale)
    println(bs)

    as = con ? power_law_a(new_xs, new_ys, bs, sis; logscale = logscale) : 
        [0.0 for i in 1:length(new_xs)]
    r = 3*N
    a = sum(as[begin+r:end-r])/ (Num - 2*r)
    b = sum(bs[begin+r:end-r])/ (Num - 2*r)
    si = sum(sis[begin+r:end-r])/ (Num - 2*r)
    return si, b, a
end
=#

function expand_left_log(xs, ys; lim = 1e-8, fit_min = 2.0,
     fit_max = 10.0, p0 = [-1.0, 1.0, 0.0], con = true)

     si, b, a = power_law_from_data(
          xs[fit_min.<xs.<fit_max], ys[fit_min.<xs.<fit_max],
          p0, fit_min, fit_max; con = con)

     i = findfirst(x -> x > fit_min, xs) - 1
     f = xs[begin] / xs[begin+1]

     new_left_xs = unique(10 .^ range(log10(lim), log10(xs[i]), step = -log10(f)))
     new_left_ys = [power_law(x, si, b, a) for x in new_left_xs]

     return new_left_xs, new_left_ys
end

function expand_right_log(xs, ys; lim = 3e3, fit_min = 5.0,
     fit_max = 10.0, p0 = [-3.0, 1.0, 0.0], con = true)

     si, b, a = power_law_from_data(
          xs[fit_min.<xs.<fit_max], ys[fit_min.<xs.<fit_max],
          p0, fit_min, fit_max; con = con)

     i = findfirst(x -> x > fit_max, xs)
     f = xs[end] / xs[end-1]

     new_right_xs = unique(10 .^ range(log10(xs[i]), log10(lim), step = log10(f)))
     new_right_ys = [power_law(x, si, b, a) for x in new_right_xs]

     return new_right_xs, new_right_ys
end



##########################################################################################92


"""
     expanded_IPS(ks, pks; k_in = 1e-8, k_end = 3e3, con = false)

Given the `ks` and `pks` of a chosen Power Spectrum, returns the same PS
with "longer tails", i.e. it is prolonged for higher and lower `ks` than 
the input ones.
"""
function expanded_IPS(ks, pks; k_in = 1e-8, k_end = 3e3, con = false)
     k1, k2 = 1e-6, 1e-4
     k3, k4 = 1e1, 2e1

     p0_beg = con ? [1.0, 1.0, 1.0] : [1.0, 1.0]
     p0_end = con ? [-3.0, 1.0, 1.0] : [-3.0, 1.0]

     #=
     si_beg, b_beg, a_beg, step_beg, fac_beg =  ks[begin]>k_in ? 
         (power_law_from_data(ks, pks, p0_beg, k1, k2; con=con)..., 
         ks[begin+1] - ks[begin] ,  ks[begin] / ks[begin+1]) :
         (nothing, nothing, nothing, nothing, nothing)


     si_end, b_end, a_end, step_end, fac_end = ks[end]<k_end ?  
         (power_law_from_data(ks, pks, p0_end, k3, k4; con=con)... ,
         ks[end] - ks[end-1], ks[end] / ks[end-1]) :
         (nothing, nothing, nothing, nothing, nothing)

     println("$si_beg , $b_beg , $a_beg , $step_beg , $fac_beg")
     println("$si_end , $b_end , $a_end , $step_end , $fac_end")
     =#

     new_left_ks, new_left_pks =
          ks[begin] > k_in ?
          expand_left_log(ks, pks; lim = k_in, fit_min = k1,
               fit_max = k2, p0 = p0_beg, con = con) : (nothing, nothing)

     new_right_ks, new_right_pks =
          ks[end] < k_end ?
          expand_right_log(ks, pks; lim = k_end, fit_min = k3,
               fit_max = k4, p0 = p0_end, con = con) : (nothing, nothing)


     new_ks, new_pks =
          if !isnothing(new_left_ks) && !isnothing(new_right_ks)
               (vcat(new_left_ks, ks[k1.<ks.<k4], new_right_ks),
                    vcat(new_left_pks, pks[k1.<ks.<k4], new_right_pks))
          elseif isnothing(new_left_ks) && !isnothing(new_right_ks)
               (vcat(ks[ks.<k4], new_right_ks),
                    vcat(pks[ks.<k4], new_right_pks))
          elseif !isnothing(new_left_ks) && isnothing(new_right_ks)
               (vcat(new_left_ks, ks[k1.<ks]),
                    vcat(new_left_pks, pks[k1.<ks]))
          else
               (ks, pks)
          end

     return new_ks, new_pks
end


"""
     expanded_Iln(PK, l, n; N = 1024, kmin = 1e-4, kmax = 1e3, s0 = 1e-3,
          fit_min = 2.0, fit_max = 10.0, p0 = [-1.0, 1.0, 0.0], con = true)


"""
function expanded_Iln(PK, l, n; lim = 1e-4, N = 1024, kmin = 1e-4, kmax = 1e3, s0 = 1e-3,
     fit_min = 2.0, fit_max = 10.0, p0 = [-1.0, 1.0, 0.0], con = true)

     rs, xis = xicalc(PK, l, n; N = N, kmin = kmin, kmax = kmax, r0 = s0)

     new_left_rs, new_left_Is = expand_left_log(rs, xis; lim = lim, fit_min = fit_min,
          fit_max = fit_max, p0 = p0, con = con)

     new_rs = vcat(new_left_rs, rs[rs.>fit_min])
     new_Is = vcat(new_left_Is, xis[rs.>fit_min])

     return new_rs, new_Is
end



##########################################################################################92



@doc raw"""
     func_I04_tilde(PK, s, kmin, kmax; kwargs...)

Return the following integral:
```math
\tilde{I}^4_0 (s) = \int_0^\infty \frac{\mathrm{d}q}{2\pi^2} 
     q^2 \, P(q) \, \frac{j_0(q s) - 1}{(q s)^4}
```

"""
function func_I04_tilde(PK, s, kmin, kmax; kwargs...)
     res = quadgk(lq -> (sphericalbesselj(0, s * exp(lq)) - 1.0) * PK(exp(lq)) / (2.0 * π^2 * exp(lq)),
          log(kmin), log(kmax); kwargs...)[1]

     return res / (s^4)
end


function expanded_I04_tilde(PK, ss;
     kmin = 1e-6, kmax = 1e3, kwargs...)

     fit_1, fit_2 = 0.1, 1.0

     if all(ss .> fit_1)
          return [func_I04_tilde(PK, s, kmin, kmax; kwargs...) for s in ss]
     else
          cutted_ss = ss[ss.>fit_1]
          cutted_I04_tildes = [func_I04_tilde(PK, s, kmin, kmax; kwargs...) for s in cutted_ss]
          l_si, l_b, l_a = GaPSE.power_law_from_data(cutted_ss, cutted_I04_tildes,
               [-2.0, -1.0, 0.0], fit_1, fit_2; con = true)
          #println("l_si, l_b, l_a = $l_si , $l_b , $l_a")
          left_I04_tildes = [GaPSE.power_law(s, l_si, l_b, l_a) for s in ss[ss.<=fit_1]]

          return vcat(left_I04_tildes, cutted_I04_tildes)
     end
end

#=
function expanded_I04_tilde(PK, ss;
     kmin = 1e-6, kmax = 1e3, kwargs...)

     fit_1, fit_2 = 0.1, 1.0
     fit_3, fit_4 = 1e3, 1e4

     if all(fit_1 .< ss .< fit_4)
          return [func_I04_tilde(PK, s, kmin, kmax; kwargs...) for s in ss]

     elseif all(ss .> fit_1)
          cutted_ss = ss[ss.<fit_4]
          cutted_I04_tildes = [func_I04_tilde(PK, s, kmin, kmax; kwargs...) for s in cutted_ss]
          r_si, r_b, r_a = GaPSE.power_law_from_data(cutted_ss, cutted_I04_tildes,
               [-4.0, -1.0, 0.0], fit_3, fit_4; con = true)
          #println("r_si, r_b, r_a = $r_si , $r_b , $r_a")
          right_I04_tildes = [GaPSE.power_law(s, r_si, r_b, r_a) for s in ss[ss.>=fit_4]]

          return vcat(cutted_I04_tildes, right_I04_tildes)

     elseif all(ss .< fit_4)
          cutted_ss = ss[ss.>fit_1]
          cutted_I04_tildes = [func_I04_tilde(PK, s, kmin, kmax; kwargs...) for s in cutted_ss]
          l_si, l_b, l_a = GaPSE.power_law_from_data(cutted_ss, cutted_I04_tildes,
               [-2.0, -1.0, 0.0], fit_1, fit_2; con = true)
          #println("l_si, l_b, l_a = $l_si , $l_b , $l_a")
          left_I04_tildes = [GaPSE.power_law(s, l_si, l_b, l_a) for s in ss[ss.<=fit_1]]

          return vcat(left_I04_tildes, cutted_I04_tildes)

     else
          cutted_ss = ss[fit_1.<ss.<fit_4]
          cutted_I04_tildes = [func_I04_tilde(PK, s, kmin, kmax; kwargs...) for s in cutted_ss]

          l_si, l_b, l_a = GaPSE.power_law_from_data(cutted_ss, cutted_I04_tildes,
               [-2.0, -1.0, 0.0], fit_1, fit_2; con = true)
          #println("l_si, l_b, l_a = $l_si , $l_b , $l_a")
          r_si, r_b, r_a = GaPSE.power_law_from_data(cutted_ss, cutted_I04_tildes,
               [-4.0, -1.0, 0.0], fit_3, fit_4; con = true)
          #println("r_si, r_b, r_a = $r_si , $r_b , $r_a")
          left_I04_tildes = [GaPSE.power_law(s, l_si, l_b, l_a) for s in ss[ss.<=fit_1]]
          right_I04_tildes = [GaPSE.power_law(s, r_si, r_b, r_a) for s in ss[ss.>=fit_4]]

          I04_tildes = vcat(left_I04_tildes, cutted_I04_tildes, right_I04_tildes)

          return I04_tildes
     end
end

=#


function my_interpolation(x1, y1, x2, y2, x)
     @assert x1 ≤ x ≤ x2 "x1 ≤ x ≤ x2 must hold!"
     x > x1 || (return y1)
     x < x2 || (return y2)
     m = (y2 - y1) / (x2 - x1)
     q = y1 - m * x1
     y = m * x + q
     return y
end
