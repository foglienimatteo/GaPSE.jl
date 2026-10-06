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

# The image ships GaPSE together with a JupyterLab that can run both the notebooks
# of `theory/` and the examples of `ipynbs/`.
#
# The base image comes from the Jupyter Docker Stacks. Since 2023 they publish on
# Quay, so the tag has to be taken from quay.io/jupyter/julia-notebook and not from
# the `jupyter/*` repositories on Docker Hub, which are no longer updated. It already
# ships Julia, IJulia and a registered Julia kernel, so only GaPSE and the plotting
# extras are added here. The Julia version is pinned to match the `julia = "1.12"`
# of `Project.toml`.

FROM quay.io/jupyter/julia-notebook:julia-1.12.7

LABEL org.opencontainers.image.title="GaPSE" \
      org.opencontainers.image.description="Galaxy Power Spectrum Estimator - a Julia package for the two-point correlation functions and power spectra of relativistic Galaxy Number Counts" \
      org.opencontainers.image.source="https://github.com/foglienimatteo/GaPSE.jl" \
      org.opencontainers.image.licenses="GPL-3.0-or-later" \
      org.opencontainers.image.version="0.10.0"

USER root

COPY --chown=${NB_UID}:${NB_GID} . /home/${NB_USER}/GaPSE
WORKDIR /home/${NB_USER}/GaPSE

USER ${NB_UID}

# PyPlot draws through matplotlib, and PyCall has to bind to the interpreter of the
# stack's conda environment instead of building a private one.
ENV PYTHON=/opt/conda/bin/python3
ENV JULIA_NUM_THREADS=auto

RUN pip install --no-cache-dir matplotlib

# The package itself. `Pkg.instantiate()` resolves `test/` as well, since the
# `[workspace]` table of `Project.toml` declares it as a member.
RUN julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'

# The extras the notebooks need. They are deliberately absent from the `[deps]` of
# `Project.toml`: `Plots`, `LaTeXStrings` and `PyPlot` are needed to redraw the
# figures, never by the library itself.
RUN julia --project=. -e 'using Pkg; Pkg.add(["Plots", "LaTeXStrings", "PyPlot"]); \
                          Pkg.build("PyCall"); Pkg.build("PyPlot"); Pkg.precompile()'

# The base image already defines the entrypoint and the command that start JupyterLab
# on port 8888, so neither is overridden here. To run the test suite instead:
#
#   docker run --rm matteofoglieni/gapse:0.10.0a \
#       julia --project=/home/jovyan/GaPSE -e 'using Pkg; Pkg.test("GaPSE")'
