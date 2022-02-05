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

#=
function power_law_a_b(x1, y1, x2, y2, si)
    b = (y2 - y1)/(x2^si - x1^si)
    a = y1 - b*(x1^si)
    return a,b
end

function power_law_a_b(xs, ys, sis)
    bs = [power_law_b(xs[i], ys[i], xs[i+1], ys[i+1], sis[i]) for i in 1:length(xs)-1]
    as = [y - b*(x^si) for (x,y,b,si) in zip(xs, ys, bs, sis)]
    return vcat(bs, bs[end]), vcat(as, as[end])
end
=#

power_law_b(x1, y1, x2, y2, si) = (y2 - y1) / (x2^si - x1^si)
function power_law_b(xs, ys, sis)
     bs = [power_law_b(xs[i], ys[i], xs[i+1], ys[i+1], sis[i]) for i in 1:length(xs)-1]
     return vcat(bs, bs[end])
end

power_law_a(x, y, b, si) = y - b * (x^si)

power_law(x, si, b, a = 0.0) = a + b * (x^si)


function power_law_from_data(xs, ys, x1, x2; N = 3, con = false)
     @assert length(xs) == length(ys) "xs and ys must have same length"
     Num = length(xs)

     sis = spectral_index(xs[x1.<xs.<x2], ys[x1.<xs.<x2];
          N = N, con = con)
     bs = power_law_b(xs[x1.<xs.<x2], ys[x1.<xs.<x2], sis)

     as = con ? power_law_a.(xs, ys, bs, sis) : [0.0 for i in 1:length(xs)]

     r = N + 1
     a = sum(as[begin+r:end-r]) / (Num - 2 * r)
     b = sum(bs[begin+r:end-r]) / (Num - 2 * r)
     si = sum(sis[begin+r:end-r]) / (Num - 2 * r)
     return si, b, a
end


function spline_IPS(ks, pks; N = 3, con = false)
     k1, k2 = 1e-6, 1e-4
     k3, k4 = 1e1, 2e1
     k_in, k_end = 1e-8, 3e3

     si_beg, b_beg, a_beg = power_law_from_data(ks, pks, k1, k2; N = N, con = con)
     si_end, b_end, a_end = power_law_from_data(ks, pks, k3, k4; N = N, con = con)

     l_beg, l_end = ks[begin+1] - ks[begin], ks[end] - ks[end-1]
     c_beg, c_end = ks[begin+1] / ks[begin], ks[end] / ks[end-1]


end
