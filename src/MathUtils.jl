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

function derivate_point(xp, yp, x1, y1, x2, y2)
     m2 = (y2 - yp) / (x2 - xp)
     m1 = (yp - y1) / (xp - x1)
     res = (m1 + m2) / 2.0
     #println(res)
     return res
end

function derivate_vector(ixs, ys; N::Integer = 1, logscale = false)
     xs = !logscale ? ixs : begin
          fac = ixs[begin+1] / ixs[begin]
          [x / (fac^i) for (i, x) in enumerate(ixs)]
     end

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


function spectral_index(ixs, ys; N::Integer = 1, con = false, logscale = false)
     xs = !logscale ? ixs : begin
          fac = ixs[begin+1] / ixs[begin]
          [x / (fac^i) for (i, x) in enumerate(ixs)]
     end
     derivs = derivate_vector(xs, ys; N = N)
     if con == false
          return [x * d / y for (x, y, d) in zip(xs, ys, derivs)]
     else
          sec_derivs = derivate_vector(xs, derivs; N = N)
          vec = [x * d2 / d for (x, d, d2) in zip(xs, derivs, sec_derivs)]
          return vec .+ 1.0
     end
end


function mean_spectral_index(ixs, ys; N::Integer = 1, con = false, logscale = false)
     vec = spectral_index(ixs, ys; N = N, con = con, logscale = logscale)[begin+2*N:end-2*N]
     return sum(vec) / length(vec)
end

power_law(x, si, b, a) = a .+ b .* (x .^ si)

function power_law_b_a(ixs, ys, si; con = false, logscale = false)
     xs = !logscale ? ixs : begin
          fac = ixs[begin+1] / ixs[begin]
          [x / (fac^i) for (i, x) in enumerate(ixs)]
     end
     if con == true
          func_fitted(x, p) = power_law(x, si, p[1], p[2])
          fit = curve_fit(func_fitted, xs, ys, [0.5, 0.5])
          return coef(fit)
     else
          func_fitted_2(x, p) = power_law(x, si, p[1], 0.0)
          fit = curve_fit(func_fitted_2, xs, ys, [1.0])
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

power_law(x, si, b, a=0.0) = a + b *(x^si)
=#

function power_law_from_data(xs, ys, p0, x1::Number, x2::Number; con=false)
    @assert length(xs) == length(ys) "xs and ys must have same length"
    #Num = length(xs)
    new_xs = xs[x1 .< xs .< x2]
    new_ys = ys[x1 .< xs .< x2]
    
    #si = mean_spectral_index(xs, ys; N=N, con=con)
    si, b, a =
        if con ==false
            @assert length(p0) == 2 " si,b to be fitted, so length(p0) must be 2!"
            vec = coef(curve_fit((x,p)-> p[2] .* x .^p[1], 
                        new_xs, new_ys, p0 ))
            vcat(vec, 0.0)
        else
            @assert length(p0) == 3 " si,b,a to be fitted, so length(p0) must be 3!"
             coef(curve_fit((x,p) -> p[3] .+ p[2] .* x .^p[1], 
                    new_xs, new_ys, p0))
        end
 
    return si, b, a
end



function power_law_from_data(xs, ys, p0; con=false)
    power_law_from_data(xs, ys, p0, xs[begin], xs[end]; con=con)
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


function expanded_IPS(ks, pks; con = false)
     k1, k2 = 1e-6, 1e-4
     k3, k4 = 1e1, 2e1
     k_in, k_end = 1e-10, 3e3

     N_k1 = findfirst(x -> x > k1, ks)
     N_k4 = findfirst(x -> x > k4, ks) - 1

     p0_beg = con ? [1.0, 1.0, 1.0] : [1.0, 1.0]
     p0_end = con ? [-3.0, 1.0, 1.0] : [-3.0, 1.0]

     si_beg, b_beg, a_beg, step_beg, fac_beg = ks[begin] > k_in ?
          (power_law_from_data(ks, pks, p0_beg, k1, k2; con = con)...,
          ks[begin+1] - ks[begin], ks[begin] / ks[begin+1]) :
          (nothing, nothing, nothing, nothing, nothing)


     si_end, b_end, a_end, step_end, fac_end = ks[end] < k_end ?
          (power_law_from_data(ks, pks, p0_end, k3, k4; con = con)...,
          ks[end] - ks[end-1], ks[end] / ks[end-1]) :
          (nothing, nothing, nothing, nothing, nothing)

     #println("$si_beg , $b_beg , $a_beg , $step_beg , $fac_beg")
     #println("$si_end , $b_end , $a_end , $step_end , $fac_end")

     new_left_ks = ks[begin] > k_in ?
          #(new_left_temp = ks[N_k1]  .- step_beg .* cumsum([fac_beg^i for i in 1:1000])
          (new_left_temp = 10 .^ range(log10(k_in), log10(ks[N_k1]), step = -log10(fac_beg));
          unique(new_left_temp[new_left_temp.>k_in])) : nothing

     new_left_pks = isnothing(new_left_ks) ? 
          nothing :
          [power_law(k, si_beg, b_beg, a_beg) for k in new_left_ks]


     new_right_ks = ks[end] < k_end ?
          (new_right_temp = 10 .^ range(log10(ks[N_k4+1]), log10(k_end), step = log10(fac_end));
          new_right_temp[new_right_temp.<k_end]) : nothing


     new_right_pks = isnothing(new_left_ks) ? nothing :
                         [power_law(k, si_end, b_end, a_end) for k in new_right_ks]


     new_ks, new_pks =
          if !isnothing(new_left_ks) && !isnothing(new_right_ks)
               (vcat(new_left_ks, ks[N_k1:N_k4], new_right_ks),
                    vcat(new_left_pks, pks[N_k1:N_k4], new_right_pks))
          elseif isnothing(new_left_ks) && !isnothing(new_right_ks)
               (vcat(ks[begin:N_k4], new_right_ks),
                    vcat(pks[begin:N_k4], new_right_pks))
          elseif !isnothing(new_left_ks) && isnothing(new_right_ks)
               (vcat(new_left_ks, ks[N_k1:end]),
                    vcat(new_left_pks, pks[N_k1:end]))
          else
               (ks, pks)
          end

     return new_ks, new_pks
end



function expanded_Iln(PK, l, n; N = 1024, kmin = 1e-4, kmax = 1e3, s0 = 1e-3,
     fit_min = 2.0, fit_max = 10.0, p0 = [-1.0, 1.0, 0.0], con = true)

     r_in = 1e-4
     rs, xis = xicalc(PK, l, n; N = N, kmin = kmin, kmax = kmax, r0 = s0)
     index = findfirst(x -> x > fit_min, rs) - 1
     fac = rs[begin] / rs[begin+1]

     si, b, a = power_law_from_data(
          rs[fit_min.<rs.<fit_max], xis[fit_min.<rs.<fit_max],
          p0, fit_min, fit_max; con = con)

     new_left_rs = (
          new_left_temp = 10 .^ range(log10(r_in), log10(rs[index]), step = -log10(fac));
          unique(new_left_temp[new_left_temp.>r_in]))
     new_left_Is = [power_law(x, si, b, a) for x in new_left_rs]

     new_rs = vcat(new_left_rs, rs[index+1:end])
     new_Is = vcat(new_left_Is, xis[index+1:end])

     return new_rs, new_Is
end

