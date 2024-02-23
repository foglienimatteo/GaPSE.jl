FROM jupyter/julia-notebook:julia-1.9.1
#FROM julia:1.8.5-bullseye

USER root

COPY . /home/jovyan/GaPSE
WORKDIR /home/jovyan/GaPSE
RUN chown -R jovyan:users /home/jovyan/GaPSE

USER jovyan
RUN pip3 install matplotlib
RUN julia --project=. --eval 'using Pkg; Pkg.instantiate(); Pkg.precompile()'
RUN julia --project=. --eval 'using Pkg; for p in [  \
        "Plots", "LaTeXStrings", "PyPlot" \
    ]; \ 
    Pkg.add(p); end; Pkg.resolve()'

#ENTRYPOINT ["/bin/bash"]


