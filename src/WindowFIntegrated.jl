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


function integrated_F_quadgk(s, μ, s_min, s_max, spline_F; 
          llim = 0.0, rlim = Inf, rtol=1e-2, atol=0.0, kwargs...)
     RLIM = isinf(rlim) ? 3.0 * s_max : rlim 
     f(p) = spline_F(s/p, μ) * ϕ(p, s_min, s_max) *
               ϕ(s2(p, s, μ), s_min, s_max) * p^2
     quadgk(p->f(p), llim, RLIM; rtol=rtol, atol=atol, kwargs...)[1]
end


function integrated_F_trapz(s, μ, s_min, s_max, spline_F; 
          llim = 0.0, rlim = Inf, N::Integer = 1000)
     RLIM = isinf(rlim) ? 3.0 * s_max : rlim 
     f(p) = spline_F(s/p, μ) * ϕ(p, s_min, s_max) *
               ϕ(s2(p, s, μ), s_min, s_max) * p^2
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
               trap::Bool = false, llim = 0.0, rlim = Inf, 
               rtol = 1e-2, atol = 0.0, N::Integer = 200)

          SS_STOP = iszero(ss_stop) ? 3.0 * s_max : ss_stop 
          ss = [s for s in ss_start:ss_step:SS_STOP]

          IFS = if trap==false
               [integrated_F_quadgk(s, μ, s_min, s_max, windowF.spline_F; 
                    llim = llim, rlim = rlim, rtol=rtol, atol=atol)
                    for s in ss, μ in windowF.μs]
          else
               [integrated_F_trapz(s, μ, s_min, s_max, windowF.spline_F; 
                    llim = llim, rlim = rlim, N = N)
                    for s in ss, μ in windowF.μs]
          end


          new(ss, windowF.μs, IFs)
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



function print_map_F(out::String, windowFint::WindowFIntegrated)
     time_1 = time()

     ss_grid = [x for x in windowFint.ss for μ in windowFint.μs]
     μs_grid = [μ for x in windowFint.ss for μ in windowFint.μs]
     IFS_grid = reshape(windowFint.IFs, (:,))

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
