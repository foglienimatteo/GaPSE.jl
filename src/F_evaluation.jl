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



@doc raw"""
    integrand_F(θ, θ_1, x, μ, θ_max; tolerance=1e-8) :: Float64

Return the integrand of the function ``F(x,\mu; \theta_\mathrm{max})``, i.e the 
function ``f(x,\mu, \theta, \theta_1; \theta_\mathrm{max})``:

```math
\begin{split}
f(x,\mu, \theta, \theta_1; \theta_\mathrm{max}) = 
    \; &\Theta\left( \frac
        {x \cos \theta + \cos \theta_1}{\sqrt{x^1+2+2x\mu}} - 
        \cos(\theta_\mathrm{max}) 
        \right) 
    \; \times \; \Theta(\mu-\cos(\theta+\theta_1)) \; \, \times \\
    & \quad \Theta(\cos(\theta - \theta_1)-\mu) \; \times \;
    \frac{4\pi \sin\theta\sin\theta_1}
        {\sqrt{(\sin\theta\sin\theta_1)^2-(\cos\theta\cos\theta_1-\mu)^2}}
\end{split}

\begin{equation}
F(x,\mu; \theta_\mathrm{max}) = \int_0^{\theta_\mathrm{max}} 
        \mathrm{d}\theta_1 \int_0^\pi \mathrm{d} \theta 
        \; f(x,\mu, \theta, \theta_1; \theta_\mathrm{max})
\end{equation}
```

`tolerance` is a parameter needed in case of small negative denominator: the Heaviside
theta function mathematically prevent that 
``\mathrm{den}=(\sin\theta\sin\theta_1)^2-(\cos\theta\cos\theta_1-\mu)^2``
becomes negative, but computationally might happen that ``\mathrm{den}`` results as a
very small negative number (for instance `-1.2368946523-18`); in this case `tolerance`
solve the problem, returning 0 if ``0<-\mathrm{den}< \mathrm{tolerance}``

See also: [`F`](@ref), [`F_map`](@ref)
"""
function integrand_F(θ, θ_1, x, μ, θ_max; tolerance = 1e-8)
     if (x * cos(θ) + cos(θ_1)) / √(x^2 + 1 + 2 * x * μ) - cos(θ_max) > 0 &&
        (μ - cos(θ + θ_1)) > 0 &&
        (cos(θ - θ_1) - μ) > 0
          den = (sin(θ) * sin(θ_1))^2 - (cos(θ) * cos(θ_1) - μ)^2
          if den > 0
               return 4.0 * π * sin(θ_1) * sin(θ) / √den
          elseif -den < tolerance
               return 0.0
          else
               throw(ErrorException("negative denominator, greater than $(tolerance): $(den)"))
          end
     else
          return 0.0
     end
end



##########################################################################################92



@doc raw"""
     F(x, μ, θ_max; tolerance=1e-8) :: Tuple{Float64, Float64}

The function ``F(x,\mu; \theta_\mathrm{max})``, defined as
follows:

```math
\begin{split}
F(x,\mu; \theta_\mathrm{max}) = & \;4\pi 
    \int_0^{\theta_\mathrm{max}} \mathrm{d}\theta_1 \int_0^\pi \mathrm{d} \theta \; 
    \, \Theta\left(\frac
        {x \cos \theta + \cos \theta_1}{\sqrt{x^1+2+2x\mu}} - 
        \cos(\theta_\mathrm{max}) 
        \right) 
    \, \Theta(\mu-\cos(\theta+\theta_1)) \\
    &\Theta(\cos(\theta - \theta_1)-\mu) \;
    \frac{\sin\theta\sin\theta_1}
        {\sqrt{(\sin\theta\sin\theta_1)^2-(\cos\theta\cos\theta_1-\mu)^2}}
\end{split}
```

`tolerance` is a parameter needed in case of small negative denominator: the Heaviside
theta function mathematically prevent that 
``\mathrm{den}=(\sin\theta\sin\theta_1)^2-(\cos\theta\cos\theta_1-\mu)^2``
becomes negative, but computationally might happen that ``\mathrm{den}`` results as a
very small negative number (for instance `-1.2368946523-18`); in this case `tolerance`
solve the problem, returning 0 if ``0<-\mathrm{den}< \mathrm{tolerance}``.

The double integral is performed with [`hcubature`](@ref) function from the Julia
Package [`HCubature`](@ref); `rtol`, `atol` and all the `kwargs` insert into `F` 
are directly transferred to `hcubature`. 

PAY ATTENTION: do not set too small `atol` and `rtol`, or the computation
can easily become overwhelming! 

See also: [`F_map`](@ref), [`integrand_F`](@ref)
"""
function F(x, μ; θ_max = π / 2.0, tolerance = 1e-8, rtol = 1e-3, atol = 1e-5, kwargs...)
     @assert 0 < tolerance < 1 "tolerance must be inside (0,1), not $(tolerance)  "
     my_int(var) = integrand_F(var[1], var[2], x, μ, θ_max; tolerance = tolerance)
     a = [0.0, 0.0]
     b = [π, θ_max]
     return hcubature(my_int, a, b; rtol = rtol, atol = atol, kwargs...)
end



##########################################################################################92


@doc raw"""
     F_map(x_step = 0.01, μ_step = 0.01; out = "data/F_map.txt", 
          x1 = 0, x2 = 3, μ1 = -1, μ2 = 1, kwargs...) 

Evaluate the window function ``F(x,\mu; \theta_\mathrm{max})`` in a grid of ``\mu``
and ``x`` values.


See also: [`F_map`](@ref), [`integrand_F`](@ref)
"""
function F_map(x_step = 0.01, μ_step = 0.01;
     out = "data/F_map.txt", x1 = 0, x2 = 3, μ1 = -1, μ2 = 1, kwargs...)

     @assert x1 >= 0.0 "The lower limit of x must be >0, not $(x1)!"
     @assert x2 > x1 "The upper limit of x must be > than the lower one, not $(x2)<=$(x1)!"
     @assert μ1 >= -1.0 "The lower limit of μ must be >=-1, not $(μ1)!"
     @assert μ2 <= 1.0 "The upper limit of μ must be <=1, not $(μ2)!"
     @assert μ2 > μ1 "The upper limit of μ must be > than the lower one, not $(μ2)<=$(μ1)!"
     @assert 0 < x_step < 1 "The integration step of x must be 0<x_step<1, not $(x_step)!"
     @assert 0 < μ_step < 1 "The integration step of μ must be 0<μ_step<1, not $(μ_step)!"

     μs = μ1:μ_step:μ2
     xs = x1:x_step:x2
     xs_grid = [x for x = xs for μ = μs]
     μs_grid = [μ for x = xs for μ = μs]

     time_1 = time()
     new_F(x, μ) = F(x, μ; kwargs...)
     Fs_grid = @showprogress map(new_F, xs_grid, μs_grid)
     time_2 = time()

     #run(`rm $(out)`)
     open(out, "w") do io
          println(io, "# Parameters used in this integration map:")
          println(io, "# x_min = $(x1) \t x_max = $(x2) \t x_step = $(x_step)")
          println(io, "# mu_min = $(μ1) \t mu_max = $(μ2) \t mu_step = $(μ_step)")
          println(io, "# computational time (in s) : $(@sprintf("%.3f", time_2-time_1))")
          print(io, "# kwards passed: ")

          if isempty(kwargs)
               println(io, "none")
          else
               print(io, "\n")
               for key in keys(kwargs)
                    println(io, "# \t\t$(key) = $(kwargs[key])")
               end
          end

          println(io, "\nx \t mu \t F \t F_error")
          for (x, μ, F) in zip(xs_grid, μs_grid, Fs_grid)
               println(io, "$x\t $μ $(F[1])\t $(F[2])")
          end
     end
end

