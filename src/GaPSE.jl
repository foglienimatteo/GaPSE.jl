module GaPSE

using TwoFAST, Dierckx
using HCubature, QuadGK, LegendrePolynomials
using GridInterpolations, PyCall, SciPy
using ProgressMeter, Printf, DelimitedFiles

h_0 = 0.70
z_MIN = 0.05
z_MAX = 0.2
θ_MAX = π / 2

FILE_PS = "data/WideA_ZA_pk.dat"
NAMES_PS = ["k (h/Mpc)", "P (Mpc/h)^3"]
FILE_BACKGROUND = "data/WideA_ZA_background.dat"
NAMES_BACKGROUND = ["z", "proper time [Gyr]", "conf. time [Mpc]", "H [1/Mpc]",
     "comov. dist.", "ang.diam.dist.", "lum. dist.", "comov.snd.hrz.",
     "(.)rho_g", "(.)rho_b", "(.)rho_cdm", "(.)rho_lambda", "(.)rho_ur",
     "(.)rho_crit", "gr.fac. D", "gr.fac. f"]

include("F-evaluation.jl")
include("Background_functions.jl")
include("Auto-doppler.jl")

end # module
