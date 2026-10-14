# Written by theory/threads_scaling.ipynb - do not edit by hand.
#
# One process of the scaling study: it builds the Cosmology, warms the JIT up, and only
# then times `map_ξ_GNC_multipole`, i.e. the function that `map_over_ss` parallelises.
#
# usage:  julia -t <N> --project=<theory> worker.jl <outdir> <outfile> <n_ss> <effect...>
#
#   <outdir>   : where to append the results
#   <outfile>  : name of the file to append to
#   <n_ss>     : comma separated sizes of the `s` grid, e.g. "3,16,64,256"
#   <effect>   : one or more GNC effect names
#
# It appends one row per (effect, n_ss) pair
#
#     threads   effect   n_ss   time[s]
#
# and, when a single effect is asked for, it also writes the whole `ξ` of the largest
# grid to `xis_<N>threads.txt`, so that the results of different thread counts can be
# compared digit by digit.

using GaPSE, Printf

const OUTDIR = abspath(ARGS[1])
const OUTFILE = ARGS[2]
const N_SS = parse.(Int, split(ARGS[3], ','))
const EFFECTS = ARGS[4:end]
const NT = Threads.nthreads()
const N_REP = 3
const PATH_TO_GAPSE = normpath(joinpath(@__DIR__, "..", ".."))

params = GaPSE.CosmoParams(0.05, 0.20, pi / 2.0;
    Ω_b = 0.0489, Ω_cdm = 0.251020, 
    h_0 = 0.70, s_lim = 1e-2, z_spline_lim = 1000.0,
    b1 = 1.0, s_b1 = 0.0, 𝑓_evo1 = 0.0, 
    b2 = nothing, s_b2 = nothing, 𝑓_evo2 = nothing,
    IPS_opts = Dict(
        :fit_left_min => 1e-6, :fit_left_max => 3e-6,
        :fit_right_min => 1e1, :fit_right_max => 2e1),
    IPSTools_opts = Dict(
        :N => 1024, :fit_min => 0.05, :fit_max => 0.5, 
        :con => true, :k_min => 1e-8, :k_max => 10.0),
);

COSMO = GaPSE.Cosmology(
    params,
    joinpath(PATH_TO_GAPSE, "data", "WideA_ZA_background.dat"),
    joinpath(PATH_TO_GAPSE, "data", "WideA_ZA_pk.dat"),
    joinpath(PATH_TO_GAPSE, "data", "F_REFERENCE_pi2.txt"),
    joinpath(PATH_TO_GAPSE, "data", "IntegrF_REFERENCE_pi2_z005020.txt")
);

const BASE_KW = Dict(
    :use_windows => false, :enhancer => 1e8, 
    :obs => :noobsvel,
    :N_trap => 30, :N_lob => 30, :N_χs => 40, :N_χs_2 => 20,
    :suit_sampling => false, :alg => :lobatto
);

# NOTE `pr = false`: the ProgressMeter is thread-safe by taking a lock as soon as it is
# advanced from more than one thread, which is exactly the kind of serialisation that a
# threading benchmark must not measure. The progress bars would go to /dev/null anyway,
# `run_scan` piping the worker's stdout there.
run_map(eff, kw, v_ss) = GaPSE.map_ξ_GNC_multipole(COSMO, eff, v_ss; L=0, pr=false, kw...)[2]

grid(n) = 10 .^ range(0, log10(2 * COSMO.s_max), length=n)

open(joinpath(OUTDIR, OUTFILE), "a") do io
    for eff in EFFECTS
        kw = GaPSE.specif_kwargs_GNC(eff, BASE_KW)
        run_map(eff, kw, grid(4))   # warm-up: everything below is measured on compiled code
        for n in N_SS
            v_ss = grid(n)
            xis = run_map(eff, kw, v_ss)
            t = minimum([(s = time(); run_map(eff, kw, v_ss); time() - s) for _ in 1:N_REP])
            @printf(io, "%3d %24s %6d %14.6f\n", NT, eff, n, t)
            flush(io)
            if length(EFFECTS) == 1 && n == maximum(N_SS)
                open(joinpath(OUTDIR, "xis_$(NT)threads.txt"), "w") do f
                    for x in xis
                        @printf(f, "%.17g\n", x)
                    end
                end
            end
        end
    end
end
println("worker done: nthreads = ", NT)
