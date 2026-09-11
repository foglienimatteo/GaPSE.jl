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
using SpecialFunctions: sphericalbesselj

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

# CAREFUL: these are NOT the integration extremes of the I_l^n!
# They are only the ones used by `IPSTools` for its own sigma_0, ..., sigma_4.
const K_MIN, K_MAX = 1e-6, 10.0

# THESE are the integration extremes of the I_l^n: `IPSTools` hard-codes
# `kmin, kmax, s0 = 1e-5, 1e3, 1e-3` (see `src/IPSTools.jl`) and hands them over to
# `xicalc`, no matter what `k_min` and `k_max` you pass to it.
# The sigma_i that appear in the asymptotic limits of the I_l^n must be computed with
# THESE extremes, otherwise the comparison is meaningless: sigma_2 barely notices the
# difference, but sigma_{-2} and sigma_{-4} change by 5 and 8 orders of magnitude.
const XICALC_KMIN, XICALC_KMAX = 1e-5, 1e3

# Comoving separations where the I_l^n will be evaluated
const SS = 10 .^ range(-6, 5, length=800)

# Comoving separations where the I_l^n will also be computed by direct quadrature.
# We stop at 1e-2 because there q*s <= 10 over the whole integration range, so the
# integrand is still perfectly resolved by `QUAD_GRID` below.
const SS_DIRECT = 10 .^ range(-6, -2, length=50)

# log-spaced grid used by `I_direct`
const QUAD_GRID = 10 .^ range(log10(XICALC_KMIN), log10(XICALC_KMAX), length=50_000)


##########################################################################################92
# The input Power Spectrum and the I_l^n

ips = GaPSE.InputPS(FILE_PS)
tools = GaPSE.IPSTools(ips; k_min=K_MIN, k_max=K_MAX, N=1024,
    fit_min=0.05, fit_max=0.5, con=true)

# P(q) evaluated once and for all on `QUAD_GRID`, so that `I_direct` is cheap
const PQ_GRID = [ips(q) for q in QUAD_GRID]


const SIGMA_CACHE = Dict{Int,Float64}()

"""
    sigma(i; kmin=XICALC_KMIN, kmax=XICALC_KMAX) ::Float64

Return the moment of the input Power Spectrum

```math
\\sigma_i = \\int_{k_\\mathrm{min}}^{k_\\mathrm{max}}
    \\frac{\\mathrm{d}q}{2 \\pi^2} \\, q^{2-i} \\, P(q) \\; .
```

`IPSTools` stores only ``\\sigma_0, ..., \\sigma_4``, while the asymptotic limits of
``I_2^0``, ``I_4^0`` and ``I_3^1`` need the negative-index ones, so we recompute them here.

The default extremes are the ones `IPSTools` uses for the ``I_\\ell^n`` themselves, and
NOT the `k_min`/`k_max` it uses for its own ``\\sigma_i``: with `WideA_ZA_pk.dat` the
difference is negligible for ``\\sigma_2`` (0.1%) but it is a factor ``5 \\times 10^4``
for ``\\sigma_{-2}`` and ``5 \\times 10^8`` for ``\\sigma_{-4}``, because those integrals
are completely dominated by their upper extreme.
"""
sigma(i; kmin=XICALC_KMIN, kmax=XICALC_KMAX) =
    get!(SIGMA_CACHE, i) do
        quadgk(q -> ips(q) * q^(2 - i) / (2 * π^2), kmin, 1e-1, 1e1, kmax)[1]
    end

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
    I_direct(l, n, s) ::Float64

Compute

```math
I_\\ell^n(s) = \\int_{k_\\mathrm{min}}^{k_\\mathrm{max}} \\frac{\\mathrm{d}q}{2\\pi^2}
    \\, q^2 \\, P(q) \\, \\frac{j_\\ell(qs)}{(qs)^n}
```

by brute-force trapezoidal quadrature on the log-spaced `QUAD_GRID`, with the same
extremes `xicalc` is given inside `IPSTools`.

This is slow and it is only usable for ``s \\lesssim 10^{-2}``, where the integrand still
has no oscillation, but it is the only way to see the true ``s \\rightarrow 0`` behaviour:
an `IntegralIPS` cannot show it, because below its `left` field it does not evaluate the
integral at all (see `plot_single`).
"""
function I_direct(l, n, s)
    integrand = [PQ_GRID[i] * QUAD_GRID[i]^3 * sphericalbesselj(l, QUAD_GRID[i] * s) /
                 (QUAD_GRID[i] * s)^n / (2 * π^2) for i in eachindex(QUAD_GRID)]
    lqs = log.(QUAD_GRID)
    return sum((integrand[i] + integrand[i+1]) * (lqs[i+1] - lqs[i]) / 2
               for i in 1:length(QUAD_GRID)-1)
end

"""
    I04_tilde_direct(s) ::Float64

Compute, by the same brute-force quadrature of `I_direct`,

```math
\\tilde{I}_0^4(s) = \\int_{k_\\mathrm{min}}^{k_\\mathrm{max}} \\frac{\\mathrm{d}q}{2\\pi^2}
    \\, q^2 \\, P(q) \\, \\frac{j_0(qs) - 1}{(qs)^4} \\; ,
```

i.e. the very same integral of `GaPSE.func_I04_tilde`.

Note that we CANNOT obtain it as `I_direct(0, 4, s) - sigma(4) / s^4`: the two terms are
equal up to ``\\mathcal{O}(s^2)``, so for small ``s`` the subtraction cancels every
significant digit. For the same reason the ratio ``(j_0(x)-1)/x^4`` is evaluated through
its series for ``x < 10^{-2}``, since `j_0(x) - 1` is pure round-off there.
"""
function I04_tilde_direct(s)
    # (j_0(x) - 1) / x^4 = -1/(6 x^2) + 1/120 - x^2/5040 + O(x^4)
    ratio(x) = x < 1e-2 ? -1 / (6 * x^2) + 1 / 120 - x^2 / 5040 :
               (sphericalbesselj(0, x) - 1) / x^4

    integrand = [PQ_GRID[i] * QUAD_GRID[i]^3 * ratio(QUAD_GRID[i] * s) / (2 * π^2)
                 for i in eachindex(QUAD_GRID)]
    lqs = log.(QUAD_GRID)
    return sum((integrand[i] + integrand[i+1]) * (lqs[i+1] - lqs[i]) / 2
               for i in 1:length(QUAD_GRID)-1)
end

"""
    asymptote(s, l, n) ::Float64

Return the leading small-``s`` behaviour of ``I_\\ell^n``:

```math
I_\\ell^n(s) \\; \\xrightarrow[s \\rightarrow 0]{} \\;
    \\frac{\\sigma_{n-\\ell}}{(2\\ell+1)!!} \\, s^{\\,\\ell-n} \\; .
```

It is reached only for ``s \\ll 1/k_\\mathrm{max} = 10^{-3}\\, h_0^{-1}\\mathrm{Mpc}``,
which is where every term of the ``j_\\ell`` series but the first becomes negligible.
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

"""
    logticks(lo, hi)

Return the decade ticks (major labelled, minor unlabelled) that fall inside `[lo, hi]`.

Generating them from the plotted range, instead of once and for all, is what keeps the
labels of the left edge from piling up on each other when a figure does not span all the
decades.
"""
function logticks(lo, hi)
    b_min, b_max = floor(Int, log10(lo)), ceil(Int, log10(hi))
    ts = [a * 10.0^b for b in b_min:b_max for a in 1:9]
    ls = [a == 1 ? L"10^{%$b}" : nothing for b in b_min:b_max for a in 1:9]
    keep = lo .<= ts .<= hi
    return (ts[keep], ls[keep])
end

# common keyword arguments, so that all the figures look the same
plot_kwargs(lo, hi) = Dict(
    :xaxis => :log, :yaxis => :log,
    :xlabel => L"s \quad [h_0^{-1}\mathrm{Mpc}]",
    :xticks => logticks(lo, hi),
    :xlims => (lo, hi),
    :legend => :bottomleft,
    :legendfontsize => 7,
    :size => (700, 470),
)

"""
    shade_extrapolations!(p, iln)

Grey out the two regions where the `iln::IntegralIPS` does NOT evaluate the integral.

An `IntegralIPS` is a spline only between its `left` (`= fit_min = 0.05` for all the
``I_\\ell^n``, `0.1` for ``\\tilde{I}_0^4``) and its `right` fields; outside them it
returns a power law ``a + b \\, s^{\\,s_i}`` whose coefficients are fitted on
``[\\mathrm{fit\\_min}, \\mathrm{fit\\_max}] = [0.05, 0.5]`` (left) and on the last 16
points of the `xicalc` grid (right).

That extrapolation has nothing to do with the true ``s \\rightarrow 0`` behaviour, and
this is the single reason why a naive plot of an `IntegralIPS` down to
``s = 10^{-4}`` looks nothing like the analytic asymptote.
"""
function shade_extrapolations!(p, iln)
    SS[begin] < iln.left &&
        vspan!(p, [SS[begin], iln.left]; color=:gray, alpha=0.13, label="")
    iln.right < SS[end] &&
        vspan!(p, [iln.right, SS[end]]; color=:gray, alpha=0.13, label="")
    vline!(p, [iln.left, iln.right]; color=:gray, ls=:dot, lw=1, label="")
    return p
end

"""
    plot_single(name, l, n, f; tilde=false)

Plot ``|I_\\ell^n(s)|`` in log-log scale and save the figure as `Iln_terms/<name>.png`.

Each figure carries four things:

 1. the `IntegralIPS` stored in `IPSTools` (solid), i.e. what GaPSE actually uses;
 2. the direct quadrature `I_direct` for ``s \\leq 10^{-2}`` (dotted), i.e. the true
    value of the integral there;
 3. the analytic small-``s`` asymptote (dashed, black);
 4. two grey bands, marking where the `IntegralIPS` is a power-law extrapolation and
    not the integral.

We plot the absolute value because all these integrals oscillate and change sign at
large ``s``, where a logarithmic vertical axis would not be defined.
"""
function plot_single(name, l, n, f; tilde=false)
    ys = [f(s) for s in SS]
    dir = tilde ? [I04_tilde_direct(s) for s in SS_DIRECT] :
          [I_direct(l, n, s) for s in SS_DIRECT]

    # The asymptote is drawn only up to s = 1: it is a pure power law, so over the whole
    # 11 decades of `SS` it would span 25 of them and squash everything else.
    ss_asy = SS[SS.<=1.0]
    asy = tilde ? [asymptote_tilde(s) for s in ss_asy] :
          [asymptote(s, l, n) for s in ss_asy]

    lab = tilde ? L"|\tilde{I}_0^4(s)| \;\; \mathrm{(IPSTools)}" :
          L"|I_{%$l}^{%$n}(s)| \;\; \mathrm{(IPSTools)}"
    asylab = tilde ? L"|-\sigma_2 / (6 s^2)|" :
             L"|\sigma_{%$(n-l)} \, s^{%$(l-n)} / %$(dfact(2 * l + 1))|"

    # the vertical range is set by the data alone, letting the asymptote clip
    vals = filter(v -> isfinite(v) && v > 0, abs.(vcat(ys, dir)))

    p = plot(; ylabel=tilde ? L"|\tilde{I}_0^4(s)|" : L"|I_{\ell}^{n}(s)|",
        title=tilde ? L"\tilde{I}_0^4" : L"I_{%$l}^{%$n}",
        ylims=(minimum(vals) / 30, maximum(vals) * 30),
        plot_kwargs(SS[begin], SS[end])...)
    plot!(p, ss_asy, abs.(asy); label=asylab, ls=:dash, lw=2, color=:black)
    plot!(p, SS, abs.(ys); label=lab, lw=2)
    plot!(p, SS_DIRECT, abs.(dir); label=L"\mathrm{direct \; quadrature}", lw=4, ls=:dot)
    shade_extrapolations!(p, f)

    savefig(p, joinpath(DIR, name * ".png"))
    SAVE_TO_DOCS && savefig(p, joinpath(DOCS_ASSETS, name * ".png"))
    return p
end

"""
    plot_all()

Plot all the ``|I_\\ell^n(s)|`` together in a single log-log figure, and save it as
`Iln_terms/all_Iln.png`.

Only the region where the `IntegralIPS` are splines is shown, i.e.
``[\\mathrm{left}, \\mathrm{right}]``: outside it they are power-law extrapolations,
and plotting them there would be misleading.
"""
function plot_all()
    left = maximum(f.left for (_, _, _, f) in ILN)
    right = minimum(f.right for (_, _, _, f) in ILN)
    ss = SS[left .<= SS .<= right]

    p = plot(; ylabel=L"|I_{\ell}^{n}(s)|", title=L"\mathrm{All \; the} \; I_{\ell}^{n}",
        plot_kwargs(left, right)...)
    for (name, l, n, f) in ILN
        plot!(p, ss, abs.([f(s) for s in ss]); label=L"I_{%$l}^{%$n}", lw=2)
    end
    plot!(p, ss, abs.([tools.I04_tilde(s) for s in ss]);
        label=L"\tilde{I}_0^4", lw=2, ls=:dot)

    savefig(p, joinpath(DIR, "all_Iln.png"))
    SAVE_TO_DOCS && savefig(p, joinpath(DOCS_ASSETS, "all_Iln.png"))
    return p
end

"""
    plot_ratios()

Plot, for each ``I_\\ell^n``, the ratio between the directly-computed integral and its
analytic small-``s`` asymptote, and save it as `Iln_terms/ratios.png`.

Every curve must tend to 1 for ``s \\rightarrow 0``: this is the actual numerical check
of the limits derived in the manual. The convergence sets in only for
``s \\lesssim 1/k_\\mathrm{max} = 10^{-3} \\, h_0^{-1}\\mathrm{Mpc}``.
"""
function plot_ratios()
    p = plot(; xaxis=:log, yaxis=:identity, ylims=(0, 1.3),
        xlabel=L"s \quad [h_0^{-1}\mathrm{Mpc}]",
        ylabel=L"I_{\ell}^{n}(s) \; / \; \mathrm{asymptote}(s)",
        title=L"\mathrm{Convergence \; to \; the} \; s \rightarrow 0 \; \mathrm{limits}",
        xticks=logticks(SS_DIRECT[begin], SS_DIRECT[end]),
        xlims=(SS_DIRECT[begin], SS_DIRECT[end]),
        legend=:bottomleft, legendfontsize=7, size=(700, 470))
    for (name, l, n, f) in ILN
        plot!(p, SS_DIRECT, [I_direct(l, n, s) / asymptote(s, l, n) for s in SS_DIRECT];
            label=L"I_{%$l}^{%$n}", lw=2)
    end
    hline!(p, [1.0]; color=:black, ls=:dash, lw=2, label="")

    savefig(p, joinpath(DIR, "ratios.png"))
    SAVE_TO_DOCS && savefig(p, joinpath(DOCS_ASSETS, "ratios.png"))
    return p
end

"""
    save_data()

Save in `Iln_terms/Iln_values.txt` a table with the comoving separations `s` and the
values of all the ``I_\\ell^n`` (and of ``\\tilde{I}_0^4``) there evaluated.

The header records both sets of ``\\sigma_i``: the ones over
``[k_\\mathrm{min}, k_\\mathrm{max}]`` that `IPSTools` stores, and the ones over
``[10^{-5}, 10^3]`` that enter the asymptotic limits.
"""
function save_data()
    out = joinpath(DIR, "Iln_values.txt")
    isfile(out) && rm(out)
    open(out, "w") do io
        println(io, GaPSE.BRAND)
        println(io, "#\n# The I_l^n integrals evaluated in the following comoving separations.")
        println(io, "# Input Power Spectrum file: $(basename(FILE_PS))")
        println(io, "#")
        println(io, "# The I_l^n are integrated over [$XICALC_KMIN, $XICALC_KMAX] (the extremes")
        println(io, "# `IPSTools` hands over to `xicalc`), and the sigma_i of the asymptotic")
        println(io, "# limits must use the same ones:")
        for i in [-4, -2, 0, 2, 4]
            println(io, "#   sigma_$i = $(sigma(i))")
        end
        println(io, "#")
        println(io, "# For comparison, the sigma_i stored by `IPSTools` (over [$K_MIN, $K_MAX]):")
        println(io, "#   sigma_0 = $(tools.σ_0) \t sigma_2 = $(tools.σ_2) \t sigma_4 = $(tools.σ_4)")
        println(io, "#")
        println(io, "# CAREFUL: outside [$(tools.I00.left), $(tools.I00.right)] " *
                    "(and [$(tools.I04_tilde.left), $(tools.I04_tilde.right)] for I04_tilde)")
        println(io, "# these values are power-law extrapolations, NOT the integrals.")
        println(io, "#")
        println(io, "# s [h_0^{-1} Mpc] \t " * join([n for (n, _, _, _) in ILN], " \t ") * " \t I04_tilde")
        for s in SS
            vals = [f(s) for (_, _, _, f) in ILN]
            println(io, "$s \t " * join(vals, " \t ") * " \t $(tools.I04_tilde(s))")
        end
    end
    return out
end

"""
    save_direct_data()

Save in `Iln_terms/Iln_direct_values.txt` the directly-computed ``I_\\ell^n(s)`` and the
ratio to their analytic asymptote, for the `SS_DIRECT` separations. Every ratio must tend
to 1 for ``s \\rightarrow 0``.
"""
function save_direct_data()
    out = joinpath(DIR, "Iln_direct_values.txt")
    isfile(out) && rm(out)
    open(out, "w") do io
        println(io, GaPSE.BRAND)
        println(io, "#\n# The I_l^n computed by direct quadrature (`I_direct`), and their ratio")
        println(io, "# to the analytic small-s asymptote sigma_{n-l} s^{l-n} / (2l+1)!! .")
        println(io, "# All the ratios must tend to 1 for s -> 0 .")
        println(io, "#")
        println(io, "# s [h_0^{-1} Mpc] \t " *
                    join([n for (n, _, _, _) in ILN], " \t ") * " \t " *
                    join([n * "_ratio" for (n, _, _, _) in ILN], " \t "))
        for s in SS_DIRECT
            vals = [I_direct(l, n, s) for (_, l, n, _) in ILN]
            rats = [v / asymptote(s, l, n) for (v, (_, l, n, _)) in zip(vals, ILN)]
            println(io, "$s \t " * join(vals, " \t ") * " \t " * join(rats, " \t "))
        end
    end
    return out
end


##########################################################################################92
# Run everything

if abspath(PROGRAM_FILE) == @__FILE__
    SAVE_TO_DOCS && mkpath(DOCS_ASSETS)

    println("\nThe I_l^n are splines only for $(tools.I00.left) <= s <= $(tools.I00.right) ,")
    println("and I04_tilde only for $(tools.I04_tilde.left) <= s <= $(tools.I04_tilde.right) :")
    println("outside, they are power-law extrapolations.\n")

    println("Plotting the single I_l^n ...")
    for (name, l, n, f) in ILN
        plot_single(name, l, n, f)
        println("\t $name done")
    end
    plot_single("I04_tilde", 0, 4, tools.I04_tilde; tilde=true)
    println("\t I04_tilde done")

    println("\nPlotting all of them together ...")
    plot_all()

    println("\nPlotting the convergence to the s -> 0 limits ...")
    plot_ratios()

    println("\nSaving the data ...")
    println("\t saved in $(save_data())")
    println("\t saved in $(save_direct_data())")

    println("\nAll the files are in $DIR .\n")
end
