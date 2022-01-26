


using Pkg
Pkg.activate(normpath(@__DIR__))
using GaPSE

FILE_NAME = split(PROGRAM_FILE, "/")[end]

#main(x::Union{String, Float64, Int64}...) = main([string(var) for var in [x...]])
function main()
     #=
     data = GaPSE.map_integral_on_mu_lensing()
     run(`rm xi_lensing.txt`)
     open("xi_lensing.txt", "w") do io
          [println(io, s, " \t ", xi) for (s, xi) in zip(data[1], data[2])]
     end

     new_data = GaPSE.PS_lensing()
     #run(`rm P_lensing.txt`)
     open("P_lensing.txt", "w") do io
          [println(io, k, " \t ", pk) for (k, pk) in zip(new_data[1], new_data[2])]
     end
     =#


     xi_doppler = "xi_doppler.txt"
     P_doppler = "P_doppler.txt"

     GaPSE.print_map_int_on_mu_doppler(xi_doppler)
     GaPSE.print_PS_doppler(P_doppler)

end

if (ARGS == String[])
     #println("\nwithout input arguments/commands, show this help message and exit\n")
     #main(["--help"])
     main()
else
     #main(ARGS)
     return 0
end
