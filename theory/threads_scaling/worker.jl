# Written by theory/threads_scaling.ipynb - do not edit by hand.
#
# One process of the scaling scan: it builds the Cosmology, warms the JIT up, and only
# then times `map_ξ_GNC_multipole`, the function `map_over_ss` parallelises.
#
# usage:  julia -t <N> --project=<theory> worker.jl <outdir> <n_ss...>

using GaPSE, Printf

const OUTDIR = abspath(ARGS[1])
const N_SS = parse.(Int, ARGS[2:end])
const NT = Threads.nthreads()
const N_REP = 5
const EFFECT = "auto_lensing"
const PATH_TO_GAPSE = normpath(joinpath(@__DIR__, "..", ".."))

params = GaPSE.CosmoParams(0.05, 0.20, pi / 2.0;
    Ω_b=0.0489, Ω_cdm=0.251020, h_0=0.7, s_lim=1e-2,
    b1=1.0, b2=nothing, s_b1=0.0, s_b2=nothing, 𝑓_evo1=0.0, 𝑓_evo2=nothing,
    IPS_opts=Dict(:fit_left_min => 1e-6, :fit_left_max => 3e-6,
        :fit_right_min => 1e1, :fit_right_max => 2e1),
    IPSTools_opts=Dict(:N => 1024, :fit_min => 0.05, :fit_max => 0.5,
        :con => true, :k_min => 1e-8, :k_max => 10.0))
COSMO = GaPSE.Cosmology(params,
    joinpath(PATH_TO_GAPSE, "data", "WideA_ZA_background.dat"),
    joinpath(PATH_TO_GAPSE, "data", "WideA_ZA_pk.dat"),
    joinpath(PATH_TO_GAPSE, "data", "F_REFERENCE_pi2.txt"),
    joinpath(PATH_TO_GAPSE, "data", "IntegrF_REFERENCE_pi2_z005020.txt"))

const KW = GaPSE.specif_kwargs_GNC(EFFECT, Dict(
    :use_windows => false, :enhancer => 1e8, :obs => :no,
    :N_trap => 30, :N_lob => 30, :N_χs => 40, :N_χs_2 => 20,
    :suit_sampling => true, :alg => :trap))

run_map(v_ss) = GaPSE.map_ξ_GNC_multipole(COSMO, EFFECT, v_ss; L=0, pr=false, KW...)[2]
grid(n) = 10 .^ range(0, log10(2 * COSMO.s_max), length=n)

run_map(grid(4))   # warm-up: everything below is measured on compiled code

open(joinpath(OUTDIR, "times.txt"), "a") do io
    for n in N_SS
        v_ss = grid(n)
        xis = run_map(v_ss)
        t = minimum([(s = time(); run_map(v_ss); time() - s) for _ in 1:N_REP])
        @printf(io, "%3d %5d %14.6f\n", NT, n, t)
        flush(io)
        # the largest grid is also written out in full, to check that the thread count
        # changes nothing but the wall time
        if n == maximum(N_SS)
            open(joinpath(OUTDIR, "xis_$(NT)threads.txt"), "w") do f
                for x in xis
                    @printf(f, "%.17g\n", x)
                end
            end
        end
    end
end
println("worker done: nthreads = ", NT)
