


using Pkg
Pkg.activate(normpath(@__DIR__))
using GaPSE

FILE_NAME = split(PROGRAM_FILE, "/")[end]

#main(x::Union{String, Float64, Int64}...) = main([string(var) for var in [x...]])
function main(args)
     println("hello!")
     return nothing
end

if (ARGS == String[])
     println("\nwithout input arguments/commands, show this help message and exit\n")
     #main(["--help"])
else
     main(ARGS)
end