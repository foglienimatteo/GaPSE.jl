FROM jupyter/julia-notebook:julia-1.9.1
#FROM julia:1.8.5-bullseye

USER root

COPY . /home/jovyan/GaPSE
WORKDIR /home/jovyan/GaPSE

USER root
RUN julia --project=. --eval 'using Pkg; Pkg.instantiate(); Pkg.precompile()'
RUN julia --project=. --eval 'using Pkg; for p in [  \
        "Plots", "LaTeXStrings", "PyPlot" \
    ]; \ 
    Pkg.add(p); end; Pkg.resolve()'
USER jovyan

#ENTRYPOINT ["/bin/bash"]


