


using Pkg
Pkg.activate(normpath(@__DIR__))
using GaPSE

FILE_NAME = split(PROGRAM_FILE, "/")[end]

#main(x::Union{String, Float64, Int64}...) = main([string(var) for var in [x...]])
function main()

     #GaPSE.print_map_int_on_mu_lensing("xi_lensing.txt"; atol=1e-5, rtol=1e-4)
     #GaPSE.print_PS("P_lensing.txt", "xi_lensing.txt")

     #GaPSE.print_map_int_on_mu_doppler("xi_doppler.txt")
     #GaPSE.print_PS("P_doppler.txt", "xi_doppler.txt")
     GaPSE.print_PS("P_doppler.txt", "auto_dopplerdtyj")

end

if (ARGS == String[])
     #println("\nwithout input arguments/commands, show this help message and exit\n")
     #main(["--help"])
     main()
else
     #main(ARGS)
     return 0
end
