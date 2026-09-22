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


# # input_ps
#
# The input matter Power Spectrum ``P(q)``, its two asymptotic power laws, and what
# they imply for the convergence of the moments
#
# ```math
#   \sigma_i = \int_{k_\mathrm{min}}^{k_\mathrm{max}} \frac{\mathrm{d}q}{2\pi^2}
#              \, q^{\,2-i} \, P(q) \; .
# ```
#
# `InputPS` reads a tabulated ``P(q)`` and, outside the tabulated range, continues
# it with a power law fitted on the first/last few points. Everything the
# ``\Delta\chi \rightarrow 0`` limits rely on - whether a given ``\sigma_i`` is a
# number or is set by where the integral is cut - is decided by those two
# exponents, which is why they deserve a page of their own.
#
# All the output files (plots and data) are saved in the `input_ps/` directory.

using Pkg
Pkg.activate(@__DIR__)

# This environment is (re)built here, so that a fresh clone - or a different
# Julia version - needs no manual setup step.
#  - GaPSE is NOT a registered package: an environment reaches it only through
#    the `path` entry that `Pkg.develop` writes into `Manifest.toml`. Checking
#    `Base.identify_package` alone is not enough: it resolves the UUID from
#    `Project.toml` and succeeds even when the source cannot be found.
#  - `Project.toml` holds names and UUIDs only, is Julia-version INdependent and
#    is tracked by git; it is rebuilt from `THEORY_DEPS` if it is missing.
#  - `Manifest.toml` holds the resolved versions, IS Julia-version specific and
#    is gitignored; `Pkg.resolve()` + `Pkg.instantiate()` rebuild it for the
#    Julia that is running. `resolve` first, or `instantiate` refuses whenever
#    `Project.toml` declares more than the manifest knows about.
# Nothing is downloaded when the environment is already complete.
let THEORY_DEPS = ["Plots", "PyPlot", "LaTeXStrings", "QuadGK",
                   "DelimitedFiles", "Printf", "SpecialFunctions"]
    declared = collect(keys(Pkg.project().dependencies))
    id = Base.identify_package("GaPSE")
    if !("GaPSE" in declared) || isnothing(id) || isnothing(Base.locate_package(id))
        @info "theory/: making the GaPSE of this repository available"
        Pkg.develop(path=dirname(@__DIR__))
    end
    todo = filter(d -> !(d in declared), THEORY_DEPS)
    isempty(todo) || (@info "theory/: adding the missing dependencies" todo; Pkg.add(todo))
    Pkg.resolve()
    Pkg.instantiate()
end

using GaPSE

using Plots, LaTeXStrings, QuadGK, DelimitedFiles, Printf

# `pyplot()` needs a working matplotlib; `gr()` is used instead where it is absent
try
    pyplot()
catch err
    @warn "theory/: PyPlot cannot be initialised; falling back to the GR backend. " *
          "The figures will look slightly different." exception = err
    gr()
end


##########################################################################################92
# Setup

const PATH_TO_GAPSE = normpath(joinpath(@__DIR__, ".."))
const FILE_PS = joinpath(PATH_TO_GAPSE, "data", "WideA_ZA_pk.dat")

const DIR = joinpath(@__DIR__, "input_ps")
@assert isdir(DIR) "ERROR: DIR=$DIR DOESN'T EXIST!!!"

const SAVE_TO_DOCS = true
const DOCS_ASSETS = joinpath(PATH_TO_GAPSE, "docs", "src", "assets", "input_ps")

const IPS = GaPSE.InputPS(FILE_PS)

const PS_TABLE = readdlm(FILE_PS, comments=true)
const K_TAB, P_TAB = PS_TABLE[:, 1], PS_TABLE[:, 2]
const K_MIN_TAB, K_MAX_TAB = extrema(K_TAB)

# The two power laws `InputPS` uses outside the tabulated range, P(q) = a + b q^s.
# `l_` is the q -> 0 side, `r_` the q -> +infinity one.
left_powerlaw(q) = IPS.l_a + IPS.l_b * q^IPS.l_si
right_powerlaw(q) = IPS.r_a + IPS.r_b * q^IPS.r_si

# The range over which GaPSE's `IPSTools` actually integrates: `xicalc` is called
# with these, whatever `k_min`/`k_max` are passed for the stored sigma_i.
const KS_XICALC = (1e-5, 1e3)


##########################################################################################92
# Plots

plot_kwargs() = Dict(
    :size => (1400, 800), :dpi => 300,
    :legendfontsize => 11, :guidefontsize => 14, :tickfontsize => 11,
    :titlefontsize => 16, :margin => 8Plots.mm,
)

logticks(lo, hi) = 10.0 .^ (floor(Int, log10(lo)):ceil(Int, log10(hi)))

"""
    plot_input_ps(; qs) :: Plots.Plot

``P(q)`` over many decades, with the two asymptotic power laws drawn through it
and the tabulated range marked. Outside that range the curve *is* the power law:
the dashed lines and the solid one lie on top of each other there, which is the
point of the figure.
"""
function plot_input_ps(; qs=10 .^ range(-8, 4, length=2000))
    p = plot(; xscale=:log10, yscale=:log10,
        xlabel=L"q \; [h \, \mathrm{Mpc}^{-1}]",
        ylabel=L"P(q) \; [h^{-3} \, \mathrm{Mpc}^3]",
        legend=:bottomleft, title="The input matter Power Spectrum",
        plot_kwargs()...)

    plot!(p, qs, [IPS(q) for q in qs], lw=2.5, c=:black, label="InputPS")
    scatter!(p, K_TAB, P_TAB, ms=1.6, mc=:orange, msw=0, label="tabulated data")

    ls = 10 .^ range(-8, log10(K_MIN_TAB) + 1.5, length=200)
    rs = 10 .^ range(log10(K_MAX_TAB) - 1.5, 4, length=200)
    plot!(p, ls, left_powerlaw.(ls), lw=2, ls=:dash, c=:red,
        label=@sprintf("left: P ~ q^{%+.3f}", IPS.l_si))
    plot!(p, rs, right_powerlaw.(rs), lw=2, ls=:dash, c=:blue,
        label=@sprintf("right: P ~ q^{%+.3f}", IPS.r_si))

    vline!(p, [K_MIN_TAB, K_MAX_TAB], lw=2, c=:gray, alpha=0.6,
        label=@sprintf("tabulated range [%.1e, %.1e]", K_MIN_TAB, K_MAX_TAB))
    xticks!(p, logticks(1e-8, 1e4))
    p
end

"""
    plot_local_slope(; qs) :: Plots.Plot

The local logarithmic slope ``\\mathrm{d}\\ln P / \\mathrm{d}\\ln q``. It is what
decides the convergence of each ``\\sigma_i``: the integrand of ``\\sigma_i`` goes
as ``q^{\\,2-i+\\mathrm{d}\\ln P/\\mathrm{d}\\ln q}``, so the moment converges at
the UV end when that exponent is below ``-1``, and at the IR end when it is above.
"""
function plot_local_slope(; qs=10 .^ range(-8, 4, length=2000))
    sl = [(log(IPS(q * 1.01)) - log(IPS(q / 1.01))) / (2 * log(1.01)) for q in qs]
    p = plot(; xscale=:log10, xlabel=L"q \; [h \, \mathrm{Mpc}^{-1}]",
        ylabel=L"\mathrm{d}\ln P / \mathrm{d}\ln q", legend=:bottomleft,
        title="Local slope of the input Power Spectrum", plot_kwargs()...)
    plot!(p, qs, sl, lw=2.5, c=:black, label="InputPS")
    hline!(p, [IPS.l_si], ls=:dash, lw=2, c=:red,
        label=@sprintf("left fit: %+.3f", IPS.l_si))
    hline!(p, [IPS.r_si], ls=:dash, lw=2, c=:blue,
        label=@sprintf("right fit: %+.3f", IPS.r_si))
    hline!(p, [-3.0], ls=:dot, lw=2, c=:green,
        label=L"-3\;:\;\mathrm{the\;CDM\;tail}\;P \sim k^{n_s-4}\ln^2 k")
    # the exponent at which the sigma_0 integrand q^2 P(q) stops converging
    hline!(p, [-3.0], ls=:dot, lw=0, label="")
    vline!(p, [K_MIN_TAB, K_MAX_TAB], lw=2, c=:gray, alpha=0.6, label="tabulated range")
    xticks!(p, logticks(1e-8, 1e4))
    p
end


##########################################################################################92
# What the two exponents imply for the sigma_i

"""
    convergence_of_the_sigmas(; io)

For each ``\\sigma_i``, the exponent of its integrand at the two ends and the
verdict. The integrand of ``\\sigma_i`` is ``q^{2-i} P(q)``, so with
``P \\sim q^{s}`` it behaves as ``q^{2-i+s}``: the integral converges at ``q \\to 0``
when ``2-i+s > -1`` and at ``q \\to \\infty`` when ``2-i+s < -1``.
"""
function convergence_of_the_sigmas(; io=stdout)
    println(io, "\n### does sigma_i converge?   (integrand q^(2-i) P(q))\n")
    @printf(io, "%-10s%14s%14s%14s%14s\n", "sigma_i",
        "IR exponent", "IR verdict", "UV exponent", "UV verdict")
    for i in 0:4
        eir, euv = 2 - i + IPS.l_si, 2 - i + IPS.r_si
        @printf(io, "%-10s%14.3f%14s%14.3f%14s\n", "sigma_$i",
            eir, eir > -1 ? "converges" : "DIVERGES",
            euv, euv < -1 ? "converges" : "DIVERGES")
    end
    println(io, "\nWith the true CDM tail P ~ k^(n_s-4) ln^2(k) ~ k^-3 the UV exponent of")
    println(io, "sigma_0 would be exactly -1, i.e. sigma_0 would be LOGARITHMICALLY divergent.")
    println(io, "The fitted right slope here is $(round(IPS.r_si, digits=4)), shallower than -3, so in")
    println(io, "GaPSE sigma_0 diverges as a power: sigma_0(<K) ~ K^$(round(3 + IPS.r_si, digits=3)).")
end

"""
    fitted_exponents(; io)

Print the two fitted power laws and the slope measured on the tabulated data.
"""
function fitted_exponents(; io=stdout)
    sl(a, b) = log(P_TAB[b] / P_TAB[a]) / log(K_TAB[b] / K_TAB[a])
    n = length(K_TAB)
    println(io, "\n### the two power laws of InputPS\n")
    @printf(io, "  left : P(q) = %.5e + %.5e * q^%+.5f\n", IPS.l_a, IPS.l_b, IPS.l_si)
    @printf(io, "  right: P(q) = %.5e + %.5e * q^%+.5f\n", IPS.r_a, IPS.r_b, IPS.r_si)
    println(io, "\n### slope measured on the tabulated data\n")
    @printf(io, "  first decade [%.3e, %.3e] : %+.4f\n", K_TAB[1], K_TAB[10], sl(1, 10))
    @printf(io, "  last decade  [%.3e, %.3e] : %+.4f\n", K_TAB[n-10], K_TAB[n], sl(n - 10, n))
end


##########################################################################################92
# Saving

function save_plot(p, name)
    savefig(p, joinpath(DIR, name))
    SAVE_TO_DOCS && isdir(DOCS_ASSETS) && savefig(p, joinpath(DOCS_ASSETS, name))
    nothing
end

function save_data(; qs=10 .^ range(-8, 4, length=2000))
    open(joinpath(DIR, "input_ps.txt"), "w") do io
        println(io, "# The input matter Power Spectrum as InputPS returns it,")
        println(io, "# together with its two asymptotic power laws.")
        println(io, "# Input file: $(basename(FILE_PS)), tabulated on [$(K_MIN_TAB), $(K_MAX_TAB)]")
        println(io, "# q \t P(q) \t left_powerlaw \t right_powerlaw")
        for q in qs
            println(io, join([q, IPS(q), left_powerlaw(q), right_powerlaw(q)], " \t "))
        end
    end
    open(joinpath(DIR, "input_ps_slopes.txt"), "w") do io
        println(io, "# Asymptotic behaviour of $(basename(FILE_PS)) and its consequences")
        fitted_exponents(io=io)
        convergence_of_the_sigmas(io=io)
    end
    nothing
end


##########################################################################################92
# Run everything

function main()
    fitted_exponents()
    convergence_of_the_sigmas()
    save_plot(plot_input_ps(), "input_ps.png")
    save_plot(plot_local_slope(), "input_ps_slope.png")
    save_data()
    println("\nDone. Plots and data are in $(DIR).")
end

main()
