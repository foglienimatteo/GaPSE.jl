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
     const DEFAULT_WFI_OPTS = Dict(
          :rlim=> nothing::Union{Float64, Nothing},
          :llim=> nothing::Union{Float64, Nothing},
          :N => 200::Int64,
          :trap => false::Bool,
          )

The default values to be used for `WindowFIntegrated`. 

See also: [`WindowFIntegrated`](@ref), [`Cosmology`](@ref)
"""
const DEFAULT_WFI_OPTS = Dict(
     :rlim=> nothing::Union{Float64, Nothing},
     :llim=> nothing::Union{Float64, Nothing},
     :N => 200::Int64,
     :trap => false::Bool,
)

function integrated_F_quadgk(s, μ, s_min, s_max, spline_F; 
          llim = nothing, rlim = nothing, rtol=1e-2, atol=0.0, kwargs...)
     LLIM = isnothing(llim) ? 0.0 : llim 
     RLIM = isnothing(rlim) ? 5.0 * s_max : rlim 
     f(p) = spline_F(s/p, μ) * ϕ(p, s_min, s_max) *
               ϕ(s2(p, s, μ), s_min, s_max) * p^2
     quadgk(p->f(p), LLIM, RLIM; rtol=rtol, atol=atol, kwargs...)[1]
end


function integrated_F_trapz(s, μ, s_min, s_max, spline_F; 
          llim = nothing, rlim = nothing, N::Integer = 1000)
     LLIM = isnothing(llim) ? 0.0 : llim 
     RLIM = isnothing(rlim) ? 5.0 * s_max : rlim 
     f(p) = spline_F(s/p, μ) * ϕ(p, s_min, s_max) *
               ϕ(s2(p, s, μ), s_min, s_max) * p^2
     ps = range(LLIM, RLIM, length = N)
     trapz(ps, f.(ps))
end


struct WindowFIntegrated
     ss::Vector{Float64}
     μs::Vector{Float64}
     IFs::Matrix{Float64}

     WindowFIntegrated(s_min, s_max, s_eff, windowF::WindowF;
               WFI_opts::Dict=Dict{Symbol,Any}())

          check_compatible_dicts(DEFAULT_WFI_OPTS, WFI_opts, "WFI_opts")
          WFI = merge(DEFAULT_WFI_OPTS, WFI_opts)

          IFS = if WFI.trap==false
               [integrated_F_quadgk(x * s_eff, μ, s_min, s_max, windowF.spline_F; 
                    llim = WFI.llim, rlim = WFI.rlim, rtol=WFI.rtol, atol=WFI.atol)
                    for x in windowF.xs, μ in windowF.μs]
          else
               [integrated_F_trapz(x * , μ, s_min, s_max, windowF.spline_F; 
                    llim = WFI.llim, rlim = WFI.rlim, N = WFI.N)
                    for x in windowF.xs, μ in windowF.μs]
          end


          new(windowF.xs .* s_eff, windowF.μs, IFs)
     end
end
