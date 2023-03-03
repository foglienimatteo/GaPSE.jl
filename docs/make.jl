push!(LOAD_PATH, "../src/")

using Documenter
using GaPSE

Documenter.makedocs(
     format=Documenter.HTML(prettyurls=get(ENV, "CI", nothing) == "true"),
     modules=[GaPSE],
     sitename="GaPSE.jl",
     pages=[
          "Introduction" => "index.md",
          "The basic structure" => [
               "Background Data" => "BackgroundData.md",
               "Cosmology Parameters" => "CosmoParams.md",
               "Cosmology Struct" => "Cosmology.md",
          ],
          "The window function F and its integration" => [
               "Window F" => "WindowF.md",
               "Integrated Window F" => "WindowFIntegrated.md",
          ],
          "Calculating TPCFs multipoles" => [
               "GNC" => "GNC_Correlations.md",
               "LD" => "LD_Correlations.md",
               "GNCxLD" => "GNCxLD_Correlations.md",
               "LDxGNC" => "GNCxLD_Correlations.md",
          ],
          "Calculating TPCFs with the PP Approximation" => "PlaneParallelApprox.md",
          "Calculating Power Spectra" => "PowerSpectra.md",
          "implication on PNG" => "PNG.md",
          "Utilities" => [
               "Dictionaries and names" => "Dicts.md",
               "Mathematical Utilities" => "MathUtils.md",
               "Cosmology Utilities" => "CosmoUtils.md",
               "Input Power Spectrum Tools" => "IPSTools.md",
               "Other Utilities" => "OtherUtils.md",
          ],
     ],
)

deploydocs(repo = "github.com/cosmofico97/GaPSE.jl.git")
