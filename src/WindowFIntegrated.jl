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


function integrated_F_quadgk(s, μ, s_min, s_max, windowF::WindowF; 
          llim = 0.0, rlim = Inf, rtol=1e-2, atol=0.0, kwargs...)
     RLIM = isinf(rlim) ? 3.0 * s_max : rlim 
     f(p) = ϕ(p, s_min, s_max) > 0 ? begin 
          spline_F(s/p, μ, windowF) * ϕ(p, s_min, s_max) * ϕ(s2(p, s, μ), s_min, s_max) * p^2
          end : 0.0
     quadgk(p->f(p), llim, RLIM; rtol=rtol, atol=atol, kwargs...)[1]
end


function integrated_F_trapz(s, μ, s_min, s_max, windowF::WindowF; 
          llim = 0.0, rlim = Inf, N::Integer = 1000)
     RLIM = isinf(rlim) ? 3.0 * s_max : rlim 
     f(p) = ϕ(p, s_min, s_max) > 0 ? begin 
          spline_F(s/p, μ, windowF) * ϕ(p, s_min, s_max) * ϕ(s2(p, s, μ), s_min, s_max) * p^2
          end : 0.0
     ps = range(llim, RLIM, length = N)
     trapz(ps, f.(ps))
end



"""
    WindowFIntegrated(
        ss::Vector{Float64}
        μs::Vector{Float64}
        IFs::Matrix{Float64}
        )

Struct containing ss, μs and IFs values of 


where ``F(x, μ)`` is the window function.

`ss` and `μs` are 1D vectors containing each value only once, while 
`IFs` values are contained in a matrix of size `(length(ss), length(μs))`, so:
- along a fixed column the changing value is `s`
- along a fixed row the changing value is `μ`

## Constructors

`WindowFIntegrated(s_min, s_max, s_eff, windowF::WindowF; WFI_opts::Dict=Dict{Symbol,Any}())` : 
take the 


See also: [`integrated_F_trapz`](@ref), [`integrated_F_quadgk`](@ref),
[`spline_integrF`](@ref), [`WindowF`](@ref)
"""
struct WindowFIntegrated
     ss::Vector{Float64}
     μs::Vector{Float64}
     IFs::Matrix{Float64}

     function WindowFIntegrated(s_min, s_max, windowF::WindowF;
               ss_start = 0.0, ss_stop = 0.0, ss_step = 21.768735478453323,
               trap::Bool = true, llim = 0.0, rlim = Inf, 
               rtol = 1e-2, atol = 0.0, N::Integer = 1000, pr::Bool = true)

          SS_STOP = iszero(ss_stop) ? 3.0 * s_max : ss_stop 
          ss = [s for s in ss_start:ss_step:SS_STOP]

          IFs = if trap==false
               pr ? begin 
                    @showprogress "calculating intF: " [
                         integrated_F_quadgk(s, μ, s_min, s_max, windowF; 
                              llim = llim, rlim = rlim, rtol=rtol, atol=atol)
                              for s in ss, μ in windowF.μs]
                    end : begin
                         [integrated_F_quadgk(s, μ, s_min, s_max, windowF; 
                              llim = llim, rlim = rlim, rtol=rtol, atol=atol)
                              for s in ss, μ in windowF.μs]
                    end
          else
               pr ? begin
                    @showprogress "calculating intF: " [
                         integrated_F_trapz(s, μ, s_min, s_max, windowF; 
                              llim = llim, rlim = rlim, N = N)
                              for s in ss, μ in windowF.μs]
               end : begin
                    [integrated_F_trapz(s, μ, s_min, s_max, windowF; 
                         llim = llim, rlim = rlim, N = N)
                         for s in ss, μ in windowF.μs]
               end
          end

          new(ss, windowF.μs, IFs)
     end

     function WindowFIntegrated(file::String)
          data = readdlm(file, comments=true)
          ss, μs, IFs = data[:, 1], data[:, 2], data[:, 3]
          @assert size(ss) == size(μs) == size(IFs) "ss, μs and IFs must have the same length!"
     
          new_ss = unique(ss)
          new_μs = unique(μs)
          new_IFs =
               if ss[2] == ss[1] && μs[2] ≠ μs[1]
                    transpose(reshape(IFs, (length(new_μs), length(new_ss))))
               elseif ss[2] ≠ ss[1] && μs[2] == μs[1]
                    reshape(IFs, (length(new_ss), length(new_μs)))
               else
                    throw(ErrorException("What kind of convenction for the file $file" *
                                         " are you using? I do not recognise it."))
               end
          new(new_ss, new_μs, new_IFs)
     end
end


"""
     spline_integrF(s, μ, str::WindowFIntegrated)::Float64

Return the 2-dim spline value of ``F`` in the given `(s,μ)`, where
``F`` is defined in the input `WindowF`.
The spline is obtained through the `interpolate` function of the 
[`GridInterpolations`](https://github.com/sisl/GridInterpolations.jl) Julia
package.

See also: [`WindowFIntegrated`](@ref)
"""
function spline_integrF(s, μ, str::WindowFIntegrated)
     grid = GridInterpolations.RectangleGrid(str.ss, str.μs)
     GridInterpolations.interpolate(grid, reshape(str.IFs, (:, 1)), [s, μ])
end


##########################################################################################92



function print_map_IntegratedF(out::String, windowFint::WindowFIntegrated)

     ss_grid = [s for s in windowFint.ss for μ in windowFint.μs]
     μs_grid = [μ for s in windowFint.ss for μ in windowFint.μs]
     IFs_grid = reshape(transpose(windowFint.IFs), (:,))

     open(out, "w") do io
          println(io, BRAND)
          println(io, "# This is an integration map of the function \\mathcal{F}(s, \\mu), defined as:")
          println(io, "# \\mathcal{F}(s, \\mu) = \\int_0^\\infty dp p^2 \\phi(p) \\phi(\\sqrt{p^2 + s^2 + 2 p s \\mu}) F(s/p, \\mu)")
          println(io, "# where F(x, \\mu) is stored in a WindowF struct (for its analytical definition, check the code.\n#")

          println(io, "#\n# s [h_0^{-1} Mpc] \t mu \t IF")
          for (s, μ, F) in zip(ss_grid, μs_grid, IFs_grid)
               println(io, "$s\t $μ \t $F")
          end
     end
end


function print_map_IntegratedF(in::String, out::String, s_min, s_max; kwargs...)
     windowF = WindowF(in)
     windowFint = WindowFIntegrated(s_min, s_max, windowF; kwargs...)
     print_map_IntegratedF(out, windowFint)
end
