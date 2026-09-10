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

# # Iln_terms
#
# Plot all the ``I_\ell^n`` integrals used by the GaPSE TPCFs, together with
# their small-``s`` asymptotic behaviours.
#
# All the output files (plots and data) are saved in the `Iln_terms/` directory.

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

# Directory where the plots and the data will be saved
const DIR = joinpath(@__DIR__, "Iln_terms")
@assert isdir(DIR) "ERROR: DIR=$DIR DOESN'T EXIST!!!"

# Set this to `true` in order to save a copy of the plots where the
# documentation expects to find them.
const SAVE_TO_DOCS = true
const DOCS_ASSETS = joinpath(PATH_TO_GAPSE, "docs", "src", "assets", "Iln_terms")

# Integration extremes of the sigma_i, the same defaults used by `IPSTools`
const K_MIN, K_MAX = 1e-6, 10.0

# Comoving separations where the I_l^n will be evaluated
const SS = 10 .^ range(-4, 4, length=600)


##########################################################################################92
# The input Power Spectrum and the I_l^n

ips = GaPSE.InputPS(FILE_PS)
tools = GaPSE.IPSTools(ips; k_min=K_MIN, k_max=K_MAX, N=1024,
    fit_min=0.05, fit_max=0.5, con=true)


"""
    sigma(i) ::Float64

Return the moment of the input Power Spectrum

```math
\\sigma_i = \\int_{k_\\mathrm{min}}^{k_\\mathrm{max}}
    \\frac{\\mathrm{d}q}{2 \\pi^2} \\, q^{2-i} \\, P(q) \\; .
```

`IPSTools` stores only ``\\sigma_0, ..., \\sigma_4``, while the asymptotic limits of
``I_2^0``, ``I_4^0`` and ``I_3^1`` need the negative-index ones, so we recompute them here.
"""
sigma(i) = quadgk(q -> ips(q) * q^(2 - i) / (2 * π^2), K_MIN, K_MAX)[1]

"""
    dfact(n) ::Int

Return the double factorial ``n!!`` (with ``n!! = 1`` for ``n \\leq 0``).
"""
dfact(n) = n <= 0 ? 1 : prod(n:-2:1)

# name, l, n, the IntegralIPS stored in `tools`
const ILN = [
    ("I00", 0, 0, tools.I00),
    ("I20", 2, 0, tools.I20),
    ("I40", 4, 0, tools.I40),
    ("I02", 0, 2, tools.I02),
    ("I22", 2, 2, tools.I22),
    ("I31", 3, 1, tools.I31),
    ("I13", 1, 3, tools.I13),
    ("I11", 1, 1, tools.I11),
]

"""
    asymptote(s, l, n) ::Float64

Return the leading small-``s`` behaviour of ``I_\\ell^n``:

```math
I_\\ell^n(s) \\; \\xrightarrow[s \\rightarrow 0]{} \\;
    \\frac{\\sigma_{n-\\ell}}{(2\\ell+1)!!} \\, s^{\\,\\ell-n} \\; .
```
"""
asymptote(s, l, n) = sigma(n - l) * s^(l - n) / dfact(2 * l + 1)

"""
    asymptote_tilde(s) ::Float64

Return the leading small-``s`` behaviour of ``\\tilde{I}_0^4``, i.e.
``-\\sigma_2 / (6 \\, s^2)``.
"""
asymptote_tilde(s) = -sigma(2) / (6 * s^2)


##########################################################################################92
# Plots

# common keyword arguments, so that all the figures look the same
const LOGTICKS = (
    vcat([a * 10.0^b for b in -4:3 for a in 1:9], 10.0^4),
    vcat([a == 1 ? L"10^{%$b}" : nothing for b in -4:3 for a in 1:9], L"10^{4}")
)

plot_kwargs() = Dict(
    :xaxis => :log, :yaxis => :log,
    :xlabel => L"s \quad [h_0^{-1}\mathrm{Mpc}]",
    :xticks => LOGTICKS,
    :legend => :bottomleft,
    :size => (600, 400),
)

"""
    plot_single(name, l, n, f; tilde=false)

Plot ``|I_\\ell^n(s)|`` in log-log scale together with its small-``s`` asymptote, and
save the figure as `Iln_terms/<name>.png`.

We plot the absolute value because all these integrals oscillate and change sign at
large ``s``, where a logarithmic vertical axis would not be defined.
"""
function plot_single(name, l, n, f; tilde=false)
    ys = [f(s) for s in SS]
    lab = tilde ? L"|\tilde{I}_0^4(s)|" : L"|I_{%$l}^{%$n}(s)|"
    asy = tilde ? [asymptote_tilde(s) for s in SS] : [asymptote(s, l, n) for s in SS]
    asylab = tilde ? L"|-\sigma_2 / (6 s^2)|" :
             L"|\sigma_{%$(n-l)} \, s^{%$(l-n)} / (2 \cdot %$l + 1)!!|"

    p = plot(SS, abs.(ys); label=lab, lw=2,
        ylabel=tilde ? L"|\tilde{I}_0^4(s)|" : L"|I_{\ell}^{n}(s)|",
        title=tilde ? L"\tilde{I}_0^4" : L"I_{%$l}^{%$n}", plot_kwargs()...)
    plot!(p, SS, abs.(asy); label=asylab, ls=:dash, lw=2, color=:black)

    savefig(p, joinpath(DIR, name * ".png"))
    SAVE_TO_DOCS && savefig(p, joinpath(DOCS_ASSETS, name * ".png"))
    return p
end

"""
    plot_all()

Plot all the ``|I_\\ell^n(s)|`` together in a single log-log figure, and save it as
`Iln_terms/all_Iln.png`.
"""
function plot_all()
    p = plot(; ylabel=L"|I_{\ell}^{n}(s)|", title=L"\mathrm{All \; the} \; I_{\ell}^{n}",
        plot_kwargs()...)
    for (name, l, n, f) in ILN
        plot!(p, SS, abs.([f(s) for s in SS]); label=L"I_{%$l}^{%$n}", lw=2)
    end
    plot!(p, SS, abs.([tools.I04_tilde(s) for s in SS]);
        label=L"\tilde{I}_0^4", lw=2, ls=:dot)

    savefig(p, joinpath(DIR, "all_Iln.png"))
    SAVE_TO_DOCS && savefig(p, joinpath(DOCS_ASSETS, "all_Iln.png"))
    return p
end

"""
    save_data()

Save in `Iln_terms/Iln_values.txt` a table with the comoving separations `s` and the
values of all the ``I_\\ell^n`` (and of ``\\tilde{I}_0^4``) there evaluated.
"""
function save_data()
    out = joinpath(DIR, "Iln_values.txt")
    isfile(out) && rm(out)
    open(out, "w") do io
        println(io, GaPSE.BRAND)
        println(io, "#\n# The I_l^n integrals evaluated in the following comoving separations.")
        println(io, "# Input Power Spectrum file: $FILE_PS")
        println(io, "# k_min = $K_MIN , k_max = $K_MAX")
        println(io, "#\n# sigma_0 = $(tools.σ_0) \t sigma_2 = $(tools.σ_2)")
        println(io, "# sigma_-2 = $(sigma(-2)) \t sigma_-4 = $(sigma(-4))")
        println(io, "#")
        println(io, "# s [h_0^{-1} Mpc] \t " * join([n for (n, _, _, _) in ILN], " \t ") * " \t I04_tilde")
        for s in SS
            vals = [f(s) for (_, _, _, f) in ILN]
            println(io, "$s \t " * join(vals, " \t ") * " \t $(tools.I04_tilde(s))")
        end
    end
    return out
end


##########################################################################################92
# Run everything

if abspath(PROGRAM_FILE) == @__FILE__
    SAVE_TO_DOCS && mkpath(DOCS_ASSETS)

    println("\nPlotting the single I_l^n ...")
    for (name, l, n, f) in ILN
        plot_single(name, l, n, f)
        println("\t $name done")
    end
    plot_single("I04_tilde", 0, 4, tools.I04_tilde; tilde=true)
    println("\t I04_tilde done")

    println("\nPlotting all of them together ...")
    plot_all()

    println("\nSaving the data ...")
    out = save_data()
    println("\t saved in $out")

    println("\nAll the files are in $DIR .\n")
end
