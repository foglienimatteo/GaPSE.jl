push!(LOAD_PATH, "../src/")

using Documenter
using GaPSE

Documenter.makedocs(
     format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
     modules = [GaPSE],
     sitename = "GaPSE.jl",
     pages = [
          "Introduction" => "index.md",
          "The F window function" => "F_evaluation.md",
          "Tool functions" => "Tool_functions.md",
          "Power Spectrum Multipoles" => "Power_Spectrum.md",
     ],
)

deploydocs(repo = "github.com/cosmofico97/GaPSE.git")
