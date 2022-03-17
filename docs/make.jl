push!(LOAD_PATH, "../src/")

using Documenter
using GaPSE

Documenter.makedocs(
     format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
     modules = [GaPSE],
     sitename = "GaPSE.jl",
     pages = [
          "Introduction" => "index.md",
          "The basic structure" => [
               "Background Data" => "BackgroundData.md",
               "Cosmology Parameters" => "CosmoParams.md",
               "Cosmology Struct" => "Cosmology.md",
          ],
          "The F window function" => "WindowF.md",
          "GR effect TPCFs" => [
               "Auto Correlations" => "AutoCorrelations.md",
               "CrossCorrelations" => "CrossCorrelations.md",
          ],
          "Calculating TPCFs and PSs" =>[
                    "Xi Multipoles" => "XiMultipoles.md",
                    "Sum Xi Multipoles" => "SumXiMultipoles.md",
                    "Power Spectrum Multipoles" => "PowerSpectrum.md",
          ],
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
