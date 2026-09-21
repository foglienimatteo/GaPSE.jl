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

# # sigma_i
#
# The ``\sigma_i`` are the moments of the input matter Power Spectrum
#
# ```math
#   \sigma_i = \int_{k_\mathrm{min}}^{k_\mathrm{max}} \frac{\mathrm{d}q}{2\pi^2}
#              \, q^{\,2-i} \, P(q) \; ,
# ```
#
# and they are what every ``\Delta\chi \rightarrow 0`` limit of the TPCFs reduces to.
# Unlike the ``I_\ell^n``, they depend explicitly on where the integral is cut, and
# not all of them converge at the same rate: this script shows which ones are safe,
# which ones are not, and where it makes sense to put the cut.
#
# All the output files (plots and data) are saved in the `sigma_i/` directory.

using Pkg
Pkg.activate(@__DIR__)
using GaPSE

using Plots, LaTeXStrings, QuadGK, DelimitedFiles, Printf

pyplot() # if you do not have PyPlot/matplotlib installed, `gr()` works as well


##########################################################################################92
# Setup

const PATH_TO_GAPSE = normpath(joinpath(@__DIR__, ".."))

# Input matter Power Spectrum P(q) at z=0
const FILE_PS = joinpath(PATH_TO_GAPSE, "data", "WideA_ZA_pk.dat")

const DIR = joinpath(@__DIR__, "sigma_i")
@assert isdir(DIR) "ERROR: DIR=$DIR DOESN'T EXIST!!!"

const SAVE_TO_DOCS = true
const DOCS_ASSETS = joinpath(PATH_TO_GAPSE, "docs", "src", "assets", "sigma_i")

# The `InputPS` extrapolates with a power law outside the range of the input file;
# these are the extremes of the tabulated data, i.e. the region where P(q) is
# actually *measured* and not extrapolated.
const PS_TABLE = readdlm(FILE_PS, comments=true)
const K_DATA_MIN, K_DATA_MAX = extrema(PS_TABLE[:, 1])

# The three pairs of extremes that appear in GaPSE, and that we want to compare.
# - IPSTOOLS: the default `k_min`/`k_max` keywords of `IPSTools`, i.e. what the
#   stored `σ_i` are computed with, and therefore what every Δχ → 0 limit uses;
# - XICALC: what `IPSTools` hard-codes for the `xicalc` call that builds the
#   `I_l^n`, and therefore what the `J * I_l^n` branch effectively uses;
# - DATA: the extremes of the input file.
const KS_IPSTOOLS = (1e-6, 10.0)
const KS_XICALC = (1e-5, 1e3)
const KS_DATA = (K_DATA_MIN, K_DATA_MAX)

const IPS = GaPSE.InputPS(FILE_PS)

# The five moments stored by `IPSTools`. `i` is the index in σ_i, so the
# integrand carries q^(2-i).
const IS = [0, 1, 2, 3, 4]


##########################################################################################92
# The integrands and the moments

"""
    sigma_integrand(q, i) :: Float64

The integrand of ``\\sigma_i``, i.e. ``q^{2-i} P(q) / (2\\pi^2)``.
"""
sigma_integrand(q, i) = IPS(q) * q^(2 - i) / (2 * π^2)

"""
    sigma(i; kmin, kmax) :: Float64

Compute ``\\sigma_i`` between `kmin` and `kmax`. The integration is split at
`1e-1` and `1e1` because the integrand spans many decades and `quadgk` converges
much faster on the three pieces than on the whole range at once.
"""
function sigma(i; kmin=KS_IPSTOOLS[1], kmax=KS_IPSTOOLS[2])
    pts = filter(x -> kmin < x < kmax, [1e-1, 1e1])
    quadgk(q -> sigma_integrand(q, i), kmin, pts..., kmax)[1]
end

"""
    cumulative_sigma(i, qs; kmin) :: Vector{Float64}

``\\sigma_i(<q)`` for every `q` in `qs`, i.e. the moment accumulated up to `q`.
Dividing it by its last value shows what fraction of the integral has already
been collected at each `q` - which is the honest way to decide where the
integral can be cut.
"""
function cumulative_sigma(i, qs; kmin=KS_XICALC[1])
    out = zeros(Float64, length(qs))
    acc, prev = 0.0, kmin
    for (n, q) in enumerate(qs)
        q <= prev && (out[n] = acc; continue)
        acc += quadgk(x -> sigma_integrand(x, i), prev, q)[1]
        out[n] = acc
        prev = q
    end
    out
end


##########################################################################################92
# Plots

plot_kwargs() = Dict(
    :size => (1400, 800), :dpi => 300,
    :legendfontsize => 11, :guidefontsize => 14, :tickfontsize => 11,
    :titlefontsize => 16, :margin => 8Plots.mm,
)

logticks(lo, hi) = 10.0 .^ (floor(Int, log10(lo)):ceil(Int, log10(hi)))

"""
    vlines_of_interest!(p)

Draw the three pairs of extremes discussed in the header, so that one can see at
a glance how much of each integrand each of them keeps.
"""
function vlines_of_interest!(p)
    vline!(p, [KS_DATA...], ls=:solid, lw=2, c=:black, alpha=0.55,
        label=@sprintf("input file: [%.1e, %.1e]", KS_DATA...))
    vline!(p, [KS_IPSTOOLS...], ls=:dash, lw=2, c=:red,
        label=@sprintf("IPSTools σ_i: [%.0e, %.0e]", KS_IPSTOOLS...))
    vline!(p, [KS_XICALC...], ls=:dot, lw=2, c=:blue,
        label=@sprintf("xicalc / I_l^n: [%.0e, %.0e]", KS_XICALC...))
    p
end

"""
    plot_integrands(; qs) :: Plots.Plot

The five integrands ``q^{2-i} P(q) / (2\\pi^2)`` on a log-log scale, with the
extremes of interest marked. The steeper the integrand at large ``q``, the more
the moment depends on where it is cut.
"""
function plot_integrands(; qs=10 .^ range(-7, 4, length=1500))
    p = plot(; xscale=:log10, yscale=:log10, xlabel=L"q \; [h \, \mathrm{Mpc}^{-1}]",
        ylabel=L"q^{2-i} \, P(q) \, / \, 2\pi^2", legend=:bottomleft,
        title="Integrands of the " * L"\sigma_i", plot_kwargs()...)
    for i in IS
        ys = [sigma_integrand(q, i) for q in qs]
        plot!(p, qs, ys, lw=2, label=L"i = %$i")
    end
    vlines_of_interest!(p)
    xticks!(p, logticks(1e-7, 1e4))
    p
end

"""
    plot_cumulatives(; qs) :: Plots.Plot

``\\sigma_i(<q) / \\sigma_i(<q_\\mathrm{max})`` for each `i`: the fraction of the
moment already collected at `q`. A curve that reaches 1 well inside the plotted
range is a moment that can be cut safely; one that is still climbing at the right
edge is a moment whose value *is* the cut.
"""
function plot_cumulatives(; qs=10 .^ range(-5, 3, length=400))
    p = plot(; xscale=:log10, xlabel=L"q \; [h \, \mathrm{Mpc}^{-1}]",
        ylabel=L"\sigma_i(<q) \, / \, \sigma_i", legend=:topleft, ylims=(-0.05, 1.15),
        title="How much of each " * L"\sigma_i" * " is collected below " * L"q",
        plot_kwargs()...)
    for i in IS
        cs = cumulative_sigma(i, qs)
        plot!(p, qs, cs ./ cs[end], lw=2, label=L"i = %$i")
    end
    hline!(p, [1.0], ls=:dash, lw=1, c=:gray, label="")
    vlines_of_interest!(p)
    xticks!(p, logticks(1e-5, 1e3))
    p
end


##########################################################################################92
# Convergence tables

"""
    convergence_table(; kmins, kmaxs, io)

Print ``\\sigma_i`` for a grid of integration extremes, so that one can see when
(and whether) each moment converges. The `kmin` scan is done at the reference
`kmax` and vice versa.
"""
function convergence_table(; kmins=[1e-7, 1e-6, 1e-5, 1e-4, 1e-3],
    kmaxs=[1e0, 1e1, 2e1, 1e2, 1e3, 1e4], io=stdout)

    println(io, "\n### sigma_i as a function of k_max   (k_min = $(KS_XICALC[1]) fixed)\n")
    print(io, @sprintf("%-12s", "k_max"))
    for i in IS
        print(io, @sprintf("%14s", "sigma_$i"))
    end
    println(io)
    for kmax in kmaxs
        print(io, @sprintf("%-12.1e", kmax))
        for i in IS
            print(io, @sprintf("%14.6e", sigma(i; kmin=KS_XICALC[1], kmax=kmax)))
        end
        println(io)
    end

    println(io, "\n### sigma_i as a function of k_min   (k_max = $(KS_IPSTOOLS[2]) fixed)\n")
    print(io, @sprintf("%-12s", "k_min"))
    for i in IS
        print(io, @sprintf("%14s", "sigma_$i"))
    end
    println(io)
    for kmin in kmins
        print(io, @sprintf("%-12.1e", kmin))
        for i in IS
            print(io, @sprintf("%14.6e", sigma(i; kmin=kmin, kmax=KS_IPSTOOLS[2])))
        end
        println(io)
    end
end

"""
    compare_the_two_ranges(; io)

The table that matters for the ``\\Delta\\chi \\rightarrow 0`` limits: the same
``\\sigma_i`` computed over the range `IPSTools` uses for them and over the range
`xicalc` uses for the ``I_\\ell^n``. Any ratio far from 1 is a place where the
limit branch and the `J * I_l^n` branch do not describe the same Power Spectrum.
"""
function compare_the_two_ranges(; io=stdout)
    println(io, "\n### the two ranges that GaPSE actually uses\n")
    a_lab = @sprintf("[%.0e, %.0e]", KS_IPSTOOLS...)
    b_lab = @sprintf("[%.0e, %.0e]", KS_XICALC...)
    println(io, @sprintf("%-10s%18s%18s%12s", "sigma_i", "IPSTools", "xicalc", "ratio"))
    println(io, @sprintf("%-10s%18s%18s%12s", "", a_lab, b_lab, "xicalc/IPSTools"))
    for i in IS
        a = sigma(i; kmin=KS_IPSTOOLS[1], kmax=KS_IPSTOOLS[2])
        b = sigma(i; kmin=KS_XICALC[1], kmax=KS_XICALC[2])
        println(io, @sprintf("%-10s%18.6e%18.6e%12.4f", "sigma_$i", a, b, b / a))
    end
end

"""
    my_sigmas(kmin, kmax; io)

Compute all the ``\\sigma_i`` between any two extremes you like, and compare them
with the `IPSTools` defaults. This is the cell to edit when you want to try a
range of your own.
"""
function my_sigmas(kmin, kmax; io=stdout)
    println(io, @sprintf("\n### sigma_i over [%.3e, %.3e]\n", kmin, kmax))
    println(io, @sprintf("%-10s%16s%16s%12s", "sigma_i", "your range", "IPSTools", "ratio"))
    for i in IS
        a = sigma(i; kmin=kmin, kmax=kmax)
        b = sigma(i; kmin=KS_IPSTOOLS[1], kmax=KS_IPSTOOLS[2])
        println(io, @sprintf("%-10s%16.6e%16.6e%12.4f", "sigma_$i", a, b, a / b))
    end
end


##########################################################################################92
# Saving

function save_plot(p, name)
    savefig(p, joinpath(DIR, name))
    SAVE_TO_DOCS && isdir(DOCS_ASSETS) && savefig(p, joinpath(DOCS_ASSETS, name))
    nothing
end

"""
    save_data(; qs)

Write the integrands and the cumulative fractions to `sigma_i/`, plus the two
tables, so that the numbers quoted in the documentation can be checked.
"""
function save_data(; qs=10 .^ range(-7, 4, length=1500))
    open(joinpath(DIR, "sigma_i_integrands.txt"), "w") do io
        println(io, "# Integrands of the sigma_i: q^(2-i) * P(q) / (2 pi^2)")
        println(io, "# Input Power Spectrum: $(basename(FILE_PS))")
        println(io, "# q \t " * join(["integrand_i=$i" for i in IS], " \t "))
        for q in qs
            println(io, join([q, (sigma_integrand(q, i) for i in IS)...], " \t "))
        end
    end
    open(joinpath(DIR, "sigma_i_tables.txt"), "w") do io
        println(io, "# sigma_i = int_kmin^kmax dq/(2 pi^2) q^(2-i) P(q)")
        println(io, "# Input Power Spectrum: $(basename(FILE_PS))")
        println(io, "# tabulated range of the input file: [$(K_DATA_MIN), $(K_DATA_MAX)]")
        compare_the_two_ranges(io=io)
        convergence_table(io=io)
    end
    nothing
end


##########################################################################################92
# Run everything

function main()
    compare_the_two_ranges()
    convergence_table()
    save_plot(plot_integrands(), "sigma_i_integrands.png")
    save_plot(plot_cumulatives(), "sigma_i_cumulatives.png")
    save_data()
    println("\nDone. Plots and data are in $(DIR).")
end

main()
