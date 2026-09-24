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

# The image ships GaPSE together with a JupyterLab that can run both the
# notebooks of `theory/` and the examples of `ipynbs/`.
#
# It is built on the official Julia image rather than on `jupyter/julia-notebook`,
# which is frozen since October 2023 and offers no Julia newer than 1.9.3: GaPSE
# declares `julia = "1.12"` in its `Project.toml`, so the Julia version has to be
# pinned here explicitly.

FROM julia:1.12-bookworm

LABEL org.opencontainers.image.title="GaPSE" \
      org.opencontainers.image.description="Galaxy Power Spectrum Estimator - a Julia package for the two-point correlation functions and power spectra of relativistic Galaxy Number Counts" \
      org.opencontainers.image.source="https://github.com/foglienimatteo/GaPSE.jl" \
      org.opencontainers.image.licenses="GPL-3.0-or-later" \
      org.opencontainers.image.version="0.10.0"

# `matplotlib` is what PyPlot draws through; `jupyterlab` provides the interface.
# `--break-system-packages` is needed because Debian bookworm marks its Python
# installation as externally managed (PEP 668), and this is a single-purpose image.
RUN apt-get update \
 && apt-get install -y --no-install-recommends python3 python3-pip git \
 && rm -rf /var/lib/apt/lists/* \
 && pip3 install --no-cache-dir --break-system-packages jupyterlab matplotlib

# Run as a non-root user, as the previous images did.
ARG NB_USER=gapse
ARG NB_UID=1000
RUN useradd --create-home --uid ${NB_UID} ${NB_USER}

ENV HOME=/home/${NB_USER}
ENV JULIA_DEPOT_PATH=${HOME}/.julia
ENV JULIA_NUM_THREADS=auto
# PyCall/PyPlot must bind to the system python we just installed, not build their own.
ENV PYTHON=/usr/bin/python3

COPY --chown=${NB_UID}:${NB_UID} . ${HOME}/GaPSE
WORKDIR ${HOME}/GaPSE
USER ${NB_USER}

# The package itself.
RUN julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'

# The extras the notebooks need, installed into the shared v1.12 environment rather
# than into `Project.toml`: `Plots`, `LaTeXStrings` and `PyPlot` are needed to redraw
# the figures, never by the library, and `IJulia` is what gives JupyterLab its Julia
# kernel. Being on the load path behind the active project, they are importable from
# any notebook without becoming a dependency of GaPSE.
RUN julia -e 'using Pkg; Pkg.add(["Plots", "LaTeXStrings", "PyPlot", "IJulia"]); \
              Pkg.build("PyCall"); Pkg.build("PyPlot"); Pkg.precompile()' \
 && julia -e 'using IJulia; IJulia.installkernel("Julia", "--project=@.")'

EXPOSE 8888
CMD ["jupyter", "lab", "--ip=0.0.0.0", "--port=8888", "--no-browser"]
