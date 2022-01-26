


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

     rtol = sqrt(eps())
     data = GaPSE.map_integral_on_mu_doppler()
     isfile(xi_doppler) && run(`rm $xi_doppler`)
     open(xi_doppler, "w") do io
          [println(io, s, " \t ", xi) for (s, xi) in zip(data[1], data[2])]
     end

     time1 = time()
     new_data = GaPSE.PS_doppler()
     time2 = time()
     diff = time2-time1
     println("\ntime needed [in s] = $diff\n")
     isfile(P_doppler) && run(`rm $P_doppler`)
     open(P_doppler, "w") do io
          [println(io, k, " \t ", pk) for (k, pk) in zip(new_data[1], new_data[2])]
     end
end

if (ARGS == String[])
     #println("\nwithout input arguments/commands, show this help message and exit\n")
     #main(["--help"])
     main()
else
     #main(ARGS)
     return 0
end
