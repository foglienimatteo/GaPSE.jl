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

# # spherical_bessels
#
# Reproduce the figures of the "Spherical Bessel Functions" page of the manual, and
# check numerically the two quantitative claims it makes:
#
#  - the small-argument expansion ``j_\ell(x) = x^\ell / (2\ell+1)!! + O(x^{\ell+2})``,
#    which is the one used to derive the small-``s`` limits of the ``I_\ell^n``;
#  - the linear regression for the first zero of ``j_\ell``.
#
# All the output files (plots and data) are saved in the `spherical_bessels/` directory.

using Pkg
Pkg.activate(@__DIR__)
using GaPSE

using Plots, LaTeXStrings, Printf
using SpecialFunctions: sphericalbesselj

pyplot() # if you do not have PyPlot/matplotlib installed, `gr()` works as well


##########################################################################################92
# Setup

const PATH_TO_GAPSE = normpath(joinpath(@__DIR__, ".."))

# Directory where the plots and the data will be saved
const DIR = joinpath(@__DIR__, "spherical_bessels")
@assert isdir(DIR) "ERROR: DIR=$DIR DOESN'T EXIST!!!"

# Set this to `true` in order to save a copy of the plots where the
# documentation expects to find them.
const SAVE_TO_DOCS = true
const DOCS_ASSETS = joinpath(PATH_TO_GAPSE, "docs", "src", "assets", "misc")

# the orders we show in the first figure, with the style used for each of them
const ORDERS = [
    (0, :red, :solid),
    (1, :blue, :solid),
    (2, :green, :dash),
    (3, :brown, :dash),
    (4, :black, :dot),
]

# the largest order for which we compute the first zero
const L_MAX = 100

# common keyword arguments, so that all the figures look the same
plot_kwargs() = Dict(
    :size => (1400, 800), :dpi => 300,
    :guidefontsize => 20, :tickfontsize => 16, :legendfontsize => 18,
    :grid => true,
)


##########################################################################################92
# The double factorial and the first zeros

"""
    dfact(n) ::Int

Return the double factorial ``n!!`` (with ``n!! = 1`` for ``n \\leq 0``).
"""
dfact(n) = n <= 0 ? 1 : prod(n:-2:1)

"""
    small_x(x, l) ::Float64

Return the leading term of the series expansion of ``j_\\ell`` near the origin:

```math
j_\\ell(x) = \\frac{x^\\ell}{(2\\ell+1)!!} \\left(1 + \\mathcal{O}(x^2)\\right)
    = x^\\ell \\left(
        \\frac{\\sqrt{\\pi}}{2^{\\ell+1} \\, \\Gamma(\\ell+3/2)} + \\mathcal{O}(x^2)
    \\right) \\; .
```

The two forms coincide because
``\\Gamma(\\ell+3/2) = \\sqrt{\\pi} \\, (2\\ell+1)!! \\, / \\, 2^{\\ell+1}``.
"""
small_x(x, l) = x^l / dfact(2 * l + 1)

"""
    first_zero(l; x_max_pad=30.0, N=200_000, rtol=1e-12) ::Float64

Return the first zero of ``j_\\ell(x)`` for ``x > 0``, found by scanning
``[\\ell + 0.1, \\ell + x_\\mathrm{max\\_pad}]`` for a sign change and then bisecting it.

There is no need for anything fancier: ``j_\\ell`` has no zero below ``\\ell`` (it is
still in its ``x^\\ell`` growth there), and its first zero sits around
``\\ell + 1.86 \\, \\ell^{1/3}``, so the scanned window always contains it.
"""
function first_zero(l; x_max_pad=30.0, N=200_000, rtol=1e-12)
    lo, hi = max(1e-6, float(l) + 0.1), float(l) + x_max_pad
    xs = range(lo, hi, length=N)
    vs = [sphericalbesselj(l, x) for x in xs]

    i = findfirst(k -> sign(vs[k]) != sign(vs[k+1]), 1:N-1)
    @assert !isnothing(i) "ERROR: no sign change of j_$l in [$lo, $hi] !!!"

    a, b = xs[i], xs[i+1]
    fa = sphericalbesselj(l, a)
    while (b - a) > rtol * b
        m = (a + b) / 2
        fm = sphericalbesselj(l, m)
        if sign(fm) == sign(fa)
            a, fa = m, fm
        else
            b = m
        end
    end
    return (a + b) / 2
end

"""
    linear_fit(xs, ys) ::Tuple{Float64,Float64}

Return the `(q, m)` of the least-squares straight line ``y = q + m \\, x``.
"""
function linear_fit(xs, ys)
    n = length(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    m = sum((xs .- mx) .* (ys .- my)) / sum((xs .- mx) .^ 2)
    return my - m * mx, m
end


##########################################################################################92
# Plots

"""
    plot_bessels()

Plot the first `length(ORDERS)` spherical Bessel functions of the first kind, and save
the figure as `spherical_bessels/spherical_bessels.png`.

This is the figure shown in the "Spherical Bessel Functions" page of the manual.
"""
function plot_bessels()
    xs = range(0, 20, length=2000)

    p = plot(; xlabel=L"\mathrm{adimensional \; argument} \; x",
        ylabel=L"\mathrm{spherical \; Bessel} \; j_{\ell}(x)",
        legend=:topright, plot_kwargs()...)
    for (l, col, sty) in ORDERS
        plot!(p, xs, [sphericalbesselj(l, x) for x in xs];
            label=L"\ell = %$l", color=col, ls=sty, lw=4)
    end

    savefig(p, joinpath(DIR, "spherical_bessels.png"))
    SAVE_TO_DOCS && savefig(p, joinpath(DOCS_ASSETS, "spherical_bessels.png"))
    return p
end

"""
    plot_small_x()

Plot ``j_\\ell(x)`` against its leading small-``x`` term in log-log scale, and save the
figure as `spherical_bessels/spherical_bessels_smallx.png`.

Every curve detaches from its asymptote around ``x \\simeq 1``, which is the reason why
the ``I_\\ell^n`` reach their own asymptotic limits only for
``s \\ll 1/k_\\mathrm{max}``: the argument of the Bessel function there is ``q s``, and
the expansion needs ``q s \\ll 1`` for every ``q`` that carries weight in the integral.
"""
function plot_small_x()
    xs = 10 .^ range(-3, 1.3, length=1000)

    p = plot(; xaxis=:log, yaxis=:log, xlabel=L"x",
        ylabel=L"|j_{\ell}(x)|", legend=:bottomright,
        ylims=(1e-18, 1e1), plot_kwargs()...)
    for (l, col, _) in ORDERS
        plot!(p, xs, abs.([sphericalbesselj(l, x) for x in xs]);
            label=L"j_{%$l}(x)", color=col, ls=:solid, lw=4)
        plot!(p, xs, [small_x(x, l) for x in xs];
            label=L"x^{%$l} / %$(dfact(2 * l + 1))", color=col, ls=:dash, lw=2)
    end

    savefig(p, joinpath(DIR, "spherical_bessels_smallx.png"))
    SAVE_TO_DOCS && savefig(p, joinpath(DOCS_ASSETS, "spherical_bessels_smallx.png"))
    return p
end

"""
    plot_first_zeros(zeros_l)

Plot the first zero of ``j_\\ell`` as a function of ``\\ell``, together with the linear
regression quoted in the manual and with the exact large-``\\ell`` expansion
``x \\simeq \\ell + 1.8557 \\, \\ell^{1/3}``, and save the figure as
`spherical_bessels/spherical_bessels_firstzeros.png`.

`zeros_l` is the vector of the first zeros for ``\\ell = 0, ..., L_\\mathrm{MAX}``, as
returned by `first_zero`.
"""
function plot_first_zeros(zeros_l)
    ls = 0:L_MAX
    q, m = linear_fit(collect(ls), zeros_l)

    p = plot(; xlabel=L"\ell", ylabel=L"\mathrm{first \; zero \; of} \; j_{\ell}(x)",
        legend=:topleft, plot_kwargs()...)
    scatter!(p, ls, zeros_l; label=L"\mathrm{numerical}", color=:black, ms=4, msw=0)
    plot!(p, ls, [q + m * l for l in ls];
        label=L"%$(round(q, digits=2)) + %$(round(m, digits=2)) \, \ell",
        color=:red, ls=:dash, lw=4)
    plot!(p, ls, [l + 1.8557 * l^(1 / 3) for l in ls];
        label=L"\ell + 1.8557 \, \ell^{1/3}", color=:blue, ls=:dot, lw=4)

    savefig(p, joinpath(DIR, "spherical_bessels_firstzeros.png"))
    SAVE_TO_DOCS && savefig(p, joinpath(DOCS_ASSETS, "spherical_bessels_firstzeros.png"))
    return p
end


##########################################################################################92
# Data

"""
    save_data(zeros_l)

Save in `spherical_bessels/first_zeros.txt` the first zero of ``j_\\ell`` for
``\\ell = 0, ..., L_\\mathrm{MAX}``, together with the two approximations plotted by
`plot_first_zeros`.
"""
function save_data(zeros_l)
    ls = 0:L_MAX
    q, m = linear_fit(collect(ls), zeros_l)

    out = joinpath(DIR, "first_zeros.txt")
    isfile(out) && rm(out)
    open(out, "w") do io
        println(io, GaPSE.BRAND)
        println(io, "#\n# The first zero of the spherical Bessel function j_l(x), for x > 0 .")
        println(io, "#")
        println(io, "# Least-squares straight line over 0 <= l <= $L_MAX :")
        println(io, "#   x = $q + $m * l")
        println(io, "#")
        println(io, "# l \t first_zero \t linear_fit \t l + 1.8557 l^(1/3)")
        for (l, z) in zip(ls, zeros_l)
            println(io, "$l \t $z \t $(q + m * l) \t $(l + 1.8557 * l^(1 / 3))")
        end
    end
    return out
end


##########################################################################################92
# Run everything

if abspath(PROGRAM_FILE) == @__FILE__
    SAVE_TO_DOCS && mkpath(DOCS_ASSETS)

    println("\nPlotting the first spherical Bessel functions ...")
    plot_bessels()

    println("\nPlotting their small-x behaviour ...")
    plot_small_x()

    println("\nComputing the first zero of j_l for l = 0, ..., $L_MAX ...")
    zeros_l = [first_zero(l) for l in 0:L_MAX]
    q, m = linear_fit(collect(0:L_MAX), zeros_l)
    @printf("\t linear regression: x = %.4f + %.4f l \n", q, m)
    @printf("\t l=0: %.4f (pi = %.4f) \t l=1: %.4f \t l=2: %.4f \t l=%d: %.4f \n",
        zeros_l[1], π, zeros_l[2], zeros_l[3], L_MAX, zeros_l[end])

    println("\nPlotting them ...")
    plot_first_zeros(zeros_l)

    println("\nSaving the data ...")
    println("\t saved in $(save_data(zeros_l))")

    println("\nAll the files are in $DIR .\n")
end
