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

function PhiTimesWindowF(s1, s, μ, windowf::WindowF; s_min, s_max)
     return ϕ(√(s1^2 + s^2 + 2 * s1 * s * μ), s_min, s_max) * spline_F(s / s1, μ, windowf)
end

function PhiTimesWindowF_multipole(
     s1, s, windowf::WindowF;
     s_min, s_max,
     L::Int=0, alg::Symbol=:lobatto,
     N_lob::Int=100, N_trap::Int=200,
     atol_quad::Float64=0.0, rtol_quad::Float64=1e-2,
     enhancer::Float64=1e6,
     kwargs...)

     @assert alg ∈ VALID_INTEGRATION_ALGORITHM ":$alg is not a valid Symbol for \"alg\"; they are: \n\t" *
                                                     "$(":".*string.(VALID_INTEGRATION_ALGORITHM) .* vcat([" , " for i in 1:length(VALID_INTEGRATION_ALGORITHM)-1], " .")... )"

     @assert N_trap > 2 "N_trap must be >2,  N_trap = $N_trap is not!"
     @assert N_lob > 2 "N_lob must be >2,  N_lob = $N_lob is not!"
     @assert atol_quad ≥ 0.0 "atol_quad must be ≥ 0.0,  atol_quad = $atol_quad is not!"
     @assert rtol_quad ≥ 0.0 "rtol_trap must be ≥ 0.0,  rtol_quad = $rtol_quad is not!"
     @assert L ≥ 0 "L must be ≥ 0, L = $L is not!"


     orig_f(μ) = enhancer * PhiTimesWindowF(s1, s, μ, windowf; s_min=s_min, s_max=s_max) * Pl(μ, L)

     int = if alg == :lobatto
          xs, ws = gausslobatto(N_lob)
          dot(ws, orig_f.(xs))

     elseif alg == :quad
          quadgk(μ -> orig_f(μ), -1.0, 1.0; atol=atol_quad, rtol=rtol_quad)[1]

     elseif alg == :trap
          μs = union(
               range(-1.0, -0.98, length=Int(ceil(N_trap / 3) + 1)),
               range(-0.98, 0.98, length=Int(ceil(N_trap / 3) + 1)),
               range(0.98, 1.0, length=Int(ceil(N_trap / 3) + 1))
          )
          #μs = range(-1.0 + 1e-6, 1.0 - 1e-6, length=N_trap)
          orig_fs = orig_f.(μs)
          trapz(μs, orig_fs)

     else
          throw(AssertionError("how the hell did you arrive here?"))
     end

     return int / enhancer
end



##########################################################################################92



function print_map_PhiTimesWindowF_multipole(
     s1_ss::Vector{Float64}, s_ss::Vector{Float64},
     windowF::Union{String,WindowF}, out::String;
     s_min, s_max,
     pr::Bool=true, L_max::Int=4, kwargs...)

     check_parent_directory(out)
     check_namefile(out)

     @assert L_max ≥ 0 "L_max must be ≥ 0!"
     @assert 0.0 < s_min < s_max " 0.0 < s_min < s_max must hold!"
     @assert all(s_ss .≥ 0.0) "All s_ss must be ≥ 0.0!"
     @assert s_ss[begin] ≈ 0.0 "Why don't you start sampling from s=0 instead from s=$(s_ss[begin])?"
     @assert all([s_ss[i+1] > s_ss[i] for i in 1:(length(s_ss)-1)]) "s_ss must be a float vector of increasing values!"
     @assert all(s1_ss .≥ 0.0) "All s1_ss must be ≥ 0.0!"
     @assert s1_ss[begin] > 0.0 "The vector s1_ss must start from a value > 0, not $(s1_ss[begin])!"
     @assert all([s1_ss[i+1] > s1_ss[i] for i in 1:(length(s1_ss)-1)]) "s1_ss must be a float vector of increasing values!"

     WINDOWF = typeof(windowF) == String ? WindowF(windowF) : windowF

     s1_ss_grid = [s1 for s1 in s1_ss for s in s_ss]
     s_ss_grid = [s for s in s1_ss for s in s_ss]

     t1 = time()
     PTFs = [zeros(length(s1_ss_grid)) for L in 0:L_max]

     if pr == true
          for L in 0:L_max
               PTFs[L+1] = @showprogress "calculating PhiTimesF L=$L: " [
                    begin
                         res = PhiTimesWindowF_multipole(
                              s1, s, WINDOWF; s_min=s_min, s_max=s_max,
                              L=L, kwargs...)

                         #println("s1 = $s1, s=$s, res = $res")
                         res
                    end for (s1, s) in zip(s1_ss_grid, s_ss_grid)]
          end
     else
          for L in 0:L_max
               PTFs[L+1] = [
                    PhiTimesWindowF_multipole(
                         s1, s, WINDOWF; s_min=s_min, s_max=s_max,
                         L=L, kwargs...)
                    for (s1, s) in zip(s1_ss_grid, s_ss_grid)]
          end
     end


     t2 = time()

     (pr) && println("\ntime needed for print_map_PhiTimesWindowF_multipole " *
                     "[in s] = $(@sprintf("%.5f", t2-t1)) \n")


     open(out, "w") do io

          println(io, BRAND)
          println(io, "# This is an integration map of the F_{l_1} multipoles, defined as:")
          println(io, "#      F_{l_1}(s_1, s \\mu) = \\int_{-1}^{+1} \\mathrm{d}\\mu \\mathcal{L}_{l_1}(\\mu) F(s_1, s, \\mu)")
          println(io, "#      F(s_1, s \\mu) =  \\phi(\\sqrt(s_1^2 + s^2 + 2 s_1 s \\mu)) F(s/s_1, \\mu)")
          println(io, "#                    =  \\sum_{l_1=0}^{\\infty} (2 l_1 + 1) \\mathcal{L}_{l_1}(\\mu) F_{l_1}(s_1, s) / 2 ")
          println(io, "# where \\mathcal{L}_{l_1}(\\mu) is thre Legendre polynomial if order l1 and")
          println(io, "# F(x, \\mu) is stored in a WindowF struct (for its analytical definition, check the code).\n#")

          println(io, "#\n# Time needed for this computation [in s]: $(t2-t1)")
          println(io, "# Range of interest:")
          println(io, "# \t s_min = $s_min h_0^{-1} Mpc")
          println(io, "# \t s_max = $s_max h_0^{-1} Mpc")
          println(io, "# The keyword arguments were:")

          if !isempty(kwargs)
               for key in keys(kwargs)
                    println(io, "# \t\t$(key) = $(kwargs[key])")
               end
          end

          println(io, "#\n# s1 [h_0^{-1} Mpc] \t s [h_0^{-1} Mpc] \t " *
                      join(["F_{l_1=$L} \t " for L in 0:L_max]))
          for (i, (s1, s)) in enumerate(zip(s1_ss_grid, s_ss_grid))
               println(io, "$s1 \t $s \t " *
                           join(["$(PTFs[L+1][i]) \t " for L in 0:L_max]))
          end

     end

end

function print_map_PhiTimesWindowF_multipole(
     s1_zs::Vector{Float64}, s_zs::Vector{Float64},
     windowF::Union{String,WindowF}, out::String,
     file_data::String; z_min, z_max,
     names_bg=NAMES_BACKGROUND, h_0=0.7, kwargs...)

     @assert 0.0 ≤ z_min < z_max "0.0 ≤ z_min < z_max must hold!"
     @assert all(s1_zs .≥ 0.0) "All s1_zs must be ≥ 0.0!"
     @assert s1_zs[begin] > 0.0 "The vector z1_ss must start from a value > 0, not $(z1_ss[begin])!"
     @assert all([s1_zs[i+1] > s1_zs[i] for i in 1:(length(s1_zs)-1)]) "s1_zs must be a float vector of increasing values!"

     @assert all(s_zs .≥ 0.0) "All s_zs must be ≥ 0.0!"
     @assert s_zs[begin] ≈ 0.0 "Why don't you start sampling from z=0 instead from z=$(s_zs[begin])?"
     @assert all([s_zs[i+1] > s_zs[i] for i in 1:(length(s_zs)-1)]) "s_zs must be a float vector of increasing values!"


     BD = BackgroundData(file_data, z_max; names=names_bg, h=h_0)
     s_of_z = Spline1D(BD.z, BD.comdist; bc="error")
     s1_ss = s_of_z.(s1_zs)
     s_ss = union([0.0], s_of_z.(s_zs[begin+1:end]))

     print_map_PhiTimesWindowF_multipole(s1_ss, s_ss,
          windowF, out; s_min=s_of_z(z_min), s_max=s_of_z(z_max), kwargs...)
end

function print_map_PhiTimesWindowF_multipole(
     windowF::Union{String,WindowF}, out::String,
     file_data::String; z_min, z_max,
     names_bg=NAMES_BACKGROUND, h_0=0.7, N_s1_ss::Int=100, N_s_ss::Int=100,
     m_s1::Float64=2.1, m_s::Float64=2.1, st::Float64=1.0, kwargs...)

     @assert 0.0 ≤ z_min < z_max "0.0 ≤ z_min < z_max must hold!"
     @assert N_s1_ss > 9 "N_s1_ss > 9 must hold!"
     @assert 0.0 < m_s1 < 10.0 "0.0 < m_s1 < 10.0 must hold!"
     @assert N_s_ss > 9 "N_s_ss > 9 must hold!"
     @assert 0.0 < m_s < 10.0 "0.0 < m_s < 10.0 must hold!"
     BD = BackgroundData(file_data, z_max; names=names_bg, h=h_0)
     s_of_z = Spline1D(BD.z, BD.comdist; bc="error")
     s_min, s_max = s_of_z(z_min), s_of_z(z_max)

     s1_ss = [s1 for s1 in range(st, m_s1 * s_max, length=N_s1_ss)]
     s_ss = union([0.0], [s for s in range(0.0, m_s * s_max, length=N_s_ss)][begin+1:end])

     print_map_PhiTimesWindowF_multipole(s1_ss, s_ss,
          windowF, out; s_min=s_min, s_max=s_max, kwargs...)
end






##########################################################################################92

##########################################################################################92

##########################################################################################92



function Q_multipole(
     s, windowf::WindowF;
     s_min, s_max,
     L::Int=0, alg::Symbol=:quad,
     llim=nothing, rlim=nothing,
     N_trap::Int=200,
     atol_quad::Float64=0.0, rtol_quad::Float64=1e-2,
     enhancer::Float64=1e6,
     in_alg::Symbol=:lobatto,
     in_N_lob::Int=100, in_N_trap::Int=200,
     in_atol_quad::Float64=0.0, in_rtol_quad::Float64=1e-2,
     in_enhancer::Float64=1e6, in_st::Float64=1.0,
     kwargs...)


     @assert N_trap > 2 "N_trap must be >2,  N_trap = $N_trap is not!"
     @assert atol_quad ≥ 0.0 "atol_quad must be ≥ 0.0,  atol_quad = $atol_quad is not!"
     @assert rtol_quad ≥ 0.0 "rtol_trap must be ≥ 0.0,  rtol_quad = $rtol_quad is not!"
     @assert L ≥ 0 "L must be ≥ 0, L = $L is not!"
     @assert isnothing(llim) || llim ≥ 0.0 "llim must be nothing or ≥ 0.0!"
     @assert isnothing(rlim) || rlim > 0.0 "rlim must be nothing or > 0.0!"
     @assert isnothing(llim) || isnothing(rlim) || rlim > llim "rlim must be > llim!"


     LLIM = isnothing(llim) ? 0.95 * s_min : llim
     RLIM = isnothing(rlim) ? 1.05 * s_max : isinf(rlim) ? 3.0 * s_max : rlim
     f(s1) = ϕ(s1, s_min, s_max) > 0 ? begin
          enhancer * s1^2 * PhiTimesWindowF_multipole(s1, s, windowf;
               s_min=s_min, s_max=s_max, L=L, alg=in_alg, N_lob=in_N_lob, N_trap=in_N_trap,
               atol_quad=in_atol_quad, rtol_quad=in_rtol_quad,
               enhancer=in_enhancer, st=in_st) * ϕ(s1, s_min, s_max)
     end : 0.0


     res = if alg == :trap
          ps = range(LLIM, RLIM, length=N_trap)
          trapz(ps, f.(ps))

     elseif alg == :quad
          quadgk(s -> f(s), LLIM, RLIM; rtol=rtol_quad, atol=atol_quad)[1]
     else
          throw(AssertionError("The value 'alg = :$alg' is not a valid algorithm; you must " *
                               "choose between ':trap' and ':quad' . "))
     end

     return res / enhancer
end



##########################################################################################92



function print_map_Q_multipole(
     s_ss::Vector{Float64},
     windowF::Union{String,WindowF}, out::String;
     s_min, s_max,
     pr::Bool=true, L_max::Int=4, kwargs...)

     check_parent_directory(out)
     check_namefile(out)

     @assert L_max ≥ 0 "L_max must be ≥ 0!"
     @assert 0.0 < s_min < s_max " 0.0 < s_min < s_max must hold!"
     @assert all(s_ss .≥ 0.0) "All s_ss must be ≥ 0.0!"
     #@assert s_ss[begin] ≈ 0.0 "Why don't you start sampling from s=0 instead from s=$(s_ss[begin])?"
     @assert all([s_ss[i+1] > s_ss[i] for i in 1:(length(s_ss)-1)]) "s_ss must be a float vector of increasing values!"

     WINDOWF = typeof(windowF) == String ? WindowF(windowF) : windowF

     t1 = time()
     Qs = [zeros(length(s_ss)) for L in 0:L_max]

     if pr == true
          for L in 0:L_max
               Qs[L+1] = @showprogress "calculating Q L=$L: " [
                    begin
                         res = Q_multipole(
                              s, WINDOWF; s_min=s_min, s_max=s_max,
                              L=L, kwargs...)

                         #println("s1 = $s1, s=$s, res = $res")
                         res
                    end for s in s_ss]
          end
     else
          for L in 0:L_max
               Qs[L+1] = [
                    Q_multipole(
                         s, WINDOWF; s_min=s_min, s_max=s_max,
                         L=L, kwargs...)
                    for s in s_ss]
          end
     end


     t2 = time()

     (pr) && println("\ntime needed for print_map_Q_multipole " *
                     "[in s] = $(@sprintf("%.5f", t2-t1)) \n")


     open(out, "w") do io

          println(io, BRAND)
          println(io, "# This is an integration map of the Q_{l_1} multipoles, defined as:")
          println(io, "#      Q_{l_1}(s_1, s \\mu) = \\int_0^{\\infty} \\mathrm{d}s_1 s_1^2 \\phi(s_1) F_{l_1}(s_1, \\mu)")
          println(io, "#      F_{l_1}(s_1, s \\mu) = \\int_{-1}^{+1} \\mathrm{d}\\mu \\mathcal{L}_{l_1}(\\mu) F(s_1, s, \\mu)")
          println(io, "#      F(s_1, s \\mu) =  \\phi(\\sqrt(s_1^2 + s^2 + 2 s_1 s \\mu)) F(s/s_1, \\mu)")
          println(io, "#                    =  \\sum_{l_1=0}^{\\infty} (2 l_1 + 1) \\mathcal{L}_{l_1}(\\mu) F_{l_1}(s_1, s) / 2 ")
          println(io, "# where \\mathcal{L}_{l_1}(\\mu) is thre Legendre polynomial if order l1 and")
          println(io, "# F(x, \\mu) is stored in a WindowF struct (for its analytical definition, check the code).\n#")

          println(io, "#\n# Time needed for this computation [in s]: $(t2-t1)")
          println(io, "# Range of interest:")
          println(io, "# \t s_min = $s_min h_0^{-1} Mpc")
          println(io, "# \t s_max = $s_max h_0^{-1} Mpc")
          println(io, "# The keyword arguments were:")

          if !isempty(kwargs)
               for key in keys(kwargs)
                    println(io, "# \t\t$(key) = $(kwargs[key])")
               end
          end

          println(io, "#\n# s [h_0^{-1} Mpc] \t " *
                      join(["Q_{l_1=$L} \t " for L in 0:L_max]))
          for (i, s) in enumerate(s_ss)
               println(io, "$s \t " *
                           join(["$(Qs[L+1][i]) \t " for L in 0:L_max]))
          end

     end

end


function print_map_Q_multipole(
     s_zs::Vector{Float64},
     windowF::Union{String,WindowF}, out::String,
     file_data::String; z_min, z_max,
     names_bg=NAMES_BACKGROUND, h_0=0.7, kwargs...)

     @assert 0.0 ≤ z_min < z_max "0.0 ≤ z_min < z_max must hold!"
     @assert all(s_zs .≥ 0.0) "All s_zs must be ≥ 0.0!"
     #@assert s_zs[begin] ≈ 0.0 "Why don't you start sampling from z=0 instead from z=$(s_zs[begin])?"
     @assert all([s_zs[i+1] > s_zs[i] for i in 1:(length(s_zs)-1)]) "s_zs must be a float vector of increasing values!"


     BD = BackgroundData(file_data, z_max; names=names_bg, h=h_0)
     s_of_z = Spline1D(BD.z, BD.comdist; bc="error")
     s_ss = s_of_z.(s1_zs)
     #s_ss = union([0.0], s_of_z.(s_zs[begin+1:end]))

     print_map_Q_multipole(s_ss,
          windowF, out; s_min=s_of_z(z_min), s_max=s_of_z(z_max), kwargs...)
end

function print_map_Q_multipole(
     windowF::Union{String,WindowF}, out::String,
     file_data::String; z_min, z_max,
     names_bg=NAMES_BACKGROUND, h_0=0.7, N::Int=100,
     st::Float64=1.0,
     m::Float64=2.1, kwargs...)

     @assert 0.0 ≤ z_min < z_max "0.0 ≤ z_min < z_max must hold!"
     @assert N > 9 "N_s_ss > 9 must hold!"
     @assert 0.0 < m < 10.0 "0.0 < m < 10.0 must hold!"
     BD = BackgroundData(file_data, z_max; names=names_bg, h=h_0)
     s_of_z = Spline1D(BD.z, BD.comdist; bc="error")
     s_min, s_max = s_of_z(z_min), s_of_z(z_max)

     s_ss = [s1 for s1 in range(st, m * s_max, length=N)]
     #s_ss = union([0.0], [s for s in range(0.0, m_s * s_max, length=N)][begin+1:end])

     print_map_Q_multipole(s_ss,
          windowF, out; s_min=s_min, s_max=s_max, kwargs...)
end