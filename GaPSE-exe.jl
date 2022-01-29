


using Pkg
Pkg.activate(normpath(@__DIR__))
using GaPSE

FILE_NAME = split(PROGRAM_FILE, "/")[end]

#main(x::Union{String, Float64, Int64}...) = main([string(var) for var in [x...]])
function main()

     xs = [x for x in 0:0.1:3]
     μs = vcat([μ for μ in -1:0.01:-0.91], [μ for μ in -0.9:0.1:0.9], [μ for μ in 0.91:0.01:1.0])
     GaPSE.F_map(xs, μs; out = "data/F_REFERENCE.txt", rtol=5e-3, atol=1e-2)

     #GaPSE.print_map_int_on_mu_lensing("outputs/xi_lensing.txt"; 
     #     χ_atol = 5e-3, χ_rtol=1e-3, atol = 1e-3, rtol = 1e-3, tol=1)
     #GaPSE.print_PS_multipole("outputs/P_lensing.txt", "outputs/xi_lensing.txt")

     #GaPSE.print_map_int_on_mu_doppler("xi_doppler.txt",
     #     atol = 1e-3, rtol = 1e-2, tol = 0.1, enhancer=1)
     #GaPSE.print_PS("P_doppler.txt", "xi_doppler.txt")
     #GaPSE.print_PS("P_doppler.txt", "auto_doppler")

end

if (ARGS == String[])
     #println("\nwithout input arguments/commands, show this help message and exit\n")
     #main(["--help"])
     main()
else
     #main(ARGS)
     return 0
end
