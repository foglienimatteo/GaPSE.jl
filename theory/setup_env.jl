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

# # setup_env
#
# Build the environment of `theory/`, so that every script and notebook of this
# directory runs on a fresh clone, and under any Julia version, with no manual
# setup step. Use it as
#
# ```julia
#   include(joinpath(@__DIR__, "setup_env.jl"))
#   setup_theory_env()
# ```
#
# Why it is needed at all:
#
# - `GaPSE` is **not a registered package**. An environment can only reach it
#   through a `path` entry, and only `Pkg.develop` writes one.
# - `Manifest.toml` holds the *resolved versions*, so it is the file that is
#   specific to a Julia version. It is gitignored, hence absent on a fresh
#   clone, and it has to be rebuilt anyway when switching between Julia
#   versions - which is exactly what `Pkg.instantiate()` below does.
# - `Project.toml` holds only names and UUIDs and is Julia-version independent,
#   so it is kept under version control: when it is there, this file is a cheap
#   no-op and nothing is downloaded. It is rebuilt from `THEORY_DEPS` when it is
#   missing, so deleting it is not fatal either.

using Pkg

const THEORY_DIR = @__DIR__
const GAPSE_ROOT = dirname(THEORY_DIR)

# What this directory needs, beside GaPSE itself. Keep it in sync with the
# `[deps]` of `theory/Project.toml`; it is the fallback used to rebuild that
# file when it is missing.
const THEORY_DEPS = ["Plots", "PyPlot", "LaTeXStrings", "QuadGK",
    "DelimitedFiles", "Printf", "SpecialFunctions"]

"""
    setup_theory_env(; deps = THEORY_DEPS) :: Nothing

Activate `theory/` and make sure that GaPSE and `deps` are usable from it.

Nothing is downloaded when the environment is already complete, so this is
cheap to call at the top of every script and notebook.
"""
function setup_theory_env(; deps=THEORY_DEPS)
    Pkg.activate(THEORY_DIR)
    declared = collect(keys(Pkg.project().dependencies))

    # GaPSE has to be both *declared* and *locatable*: it is declared by
    # `Project.toml`, but what makes it locatable is the `path` entry that
    # `Pkg.develop` writes into `Manifest.toml`. Checking only the first is not
    # enough - `Base.identify_package` resolves the UUID from `Project.toml`
    # and succeeds even when the source cannot be found.
    id = Base.identify_package("GaPSE")
    if !("GaPSE" in declared) || isnothing(id) || isnothing(Base.locate_package(id))
        @info "theory/: making the GaPSE of this repository available" GAPSE_ROOT
        Pkg.develop(path=GAPSE_ROOT)
    end

    todo = filter(d -> !(d in declared), deps)
    if !isempty(todo)
        @info "theory/: adding the dependencies that are not declared yet" todo
        Pkg.add(todo)
    end

    # `resolve` first: when `Project.toml` declares more than the manifest knows
    # about - e.g. after restoring it from git over a manifest built by an
    # earlier run - `instantiate` alone refuses with "`X` is a direct
    # dependency, but does not appear in the manifest".
    Pkg.resolve()
    # then install, and resolve the versions for *this* Julia
    Pkg.instantiate()
    return nothing
end

"""
    use_pyplot_or_gr() :: Nothing

Select `pyplot()` as the `Plots` backend, falling back to `gr()` when PyPlot
cannot start - it needs a working `matplotlib`, which is not there on every
machine. Call it *after* `using Plots`.
"""
function use_pyplot_or_gr()
    try
        pyplot()
    catch err
        @warn "theory/: PyPlot cannot be initialised (it needs a working matplotlib); " *
              "falling back to the GR backend. The figures will look slightly different." exception = (err, catch_backtrace())
        gr()
    end
    return nothing
end
