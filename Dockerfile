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
# of "theory/" and the examples of "ipynbs/".
#
# The base image comes from the Jupyter Docker Stacks. 
# Since 2023 they publish on Quay, so the tag has to be taken from quay.io/jupyter/julia-notebook 
# and not from the "jupyter/*" repositories on Docker Hub, which are no longer updated. 
# This image already ships Julia, IJulia and a registered Julia kernel, so 
# only GaPSE and the plotting extras are added here.

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

# "auto" is the size of the CPU affinity mask, not the host core count: 
# it follows "docker run --cpuset-cpus" but NOT "docker run --cpus", which is a CFS quota and leaves the mask untouched
# => with "--cpus" (if you really cannot use "--cpuset-cpus") pass "-e JULIA_NUM_THREADS=<N>" as well, 
#    or the threads oversubscribe the cores the container was actually given
ENV JULIA_NUM_THREADS=auto


# GaPSE precompilation. 
# "Pkg.instantiate()" resolves "test/" too, since the "[workspace]" table of
# "Project.toml" declares it as a member, so the test suite can be run inside the container.
RUN julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'

# NOTE: The notebooks do NOT run in the environment above!
# The kernel of the base image carries "--project=@."
# => a notebook in "theory/" activates "theory/Project.toml" and one in "ipynbs/" activates "ipynbs/Project.toml"
# => those are the environments to install, and they are already complete: "Plots", "LaTeXStrings", "Dierckx" etc 
#    are declared there, and deliberately absent from the "[deps]" of GaPSE, which never plots anything.
# "Pkg.resolve()" is the safety net for a "Manifest.toml" that
# entered the build context anyway: see the ".dockerignore".
RUN julia --project=ipynbs -e 'using Pkg; Pkg.resolve(); Pkg.instantiate(); Pkg.precompile()'
RUN julia --project=theory -e 'using Pkg; Pkg.resolve(); Pkg.instantiate(); Pkg.precompile()'

# The shared environment of this Julia (no "--project", i.e. "~/.julia/environments/v1.12"),
# for the notebooks the user writes OUTSIDE the two directories above. 
# With "--project=@.", a notebook in a directory with no "Project.toml" above it lands here, 
# and "@v#.#" IS the fallback for a package the active project does not declare
# => "using GaPSE, Plots" works out of the box in, say, "/home/jovyan/work/", with no environment to set up first
RUN julia -e 'using Pkg; Pkg.develop(path=pwd()); Pkg.add(["Plots", "LaTeXStrings"]); Pkg.precompile()'

# No Python here on purpose: every notebook selects "gr()", the GR backend needs no Python at
# all, and "PyPlot"/"pythonplot()" were dropped under Julia 1.12 (see the CHANGELOG). Should a
# notebook ever need matplotlib again uncomment the following lines:
#     # PyPlot and PythonPlot draws through matplotlib, and they have to bind to the interpreter of the
#     # stack's conda environment instead of building a private one.
#     ENV PYTHON=/opt/conda/bin/python3
#     RUN pip install --no-cache-dir matplotlib
#     RUN julia -e 'Pkg.add("PythonPlot"); Pkg.build("PythonPlot"); Pkg.instantiate(); Pkg.precompile()'

# The base image already defines the entrypoint and the command that start JupyterLab
# on port 8888, so neither has to be overridden here.
# To run this container:
#
#   $ docker run --detach --env JUPYTER_TOKEN=test --name gapse --publish 10000:8888 foglienimatteo/gapse:latest
#
# To run the test suite instead (with a limited number of cores):
#
#   $ docker run --rm --cpuset-cpus="0-1" foglienimatteo/gapse:latest julia --project=/home/jovyan/GaPSE -e 'using Pkg; Pkg.test("GaPSE")'
#
# "JULIA_NUM_THREADS=auto" is in the environment of the image and "Pkg.test" forwards it to the
# test worker, so it will run the suite on 2 threads.

