


using Pkg
Pkg.activate(normpath(@__DIR__))
using GaPSE

FILE_NAME = split(PROGRAM_FILE, "/")[end]

#main(x::Union{String, Float64, Int64}...) = main([string(var) for var in [x...]])
function main()
     data = GaPSE.PS_doppler()
     open("my_first_P_doppler.txt", "w") do io
          [println(io, k, " \t ", pk) for (k, pk) in zip(data[1], data[2])]
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