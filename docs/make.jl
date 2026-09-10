push!(LOAD_PATH, "../src/")

using Documenter
using GaPSE

Documenter.makedocs(
     format=Documenter.HTML(prettyurls=get(ENV, "CI", nothing) == "true"),
     modules=[GaPSE],
     sitename="GaPSE.jl",
     pages=[
          "Introduction" => "index.md",
          "Theory" => [
               "Spline Theory" => "SplineTheory.md"
          ],
          "The code basic structures" => [
               "Background Data" => "BackgroundData.md",
               "Cosmology Parameters" => "CosmoParams.md",
               "Cosmology Struct" => "Cosmology.md",
          ],
          "The window function F and its integration" => [
               "Window F" => "WindowF.md",
               "Integrated Window F" => "WindowFIntegrated.md",
          ],
          "TPCFs multipoles" => [
               "GNC correlations" => "GNC_Correlations_1-2.md",
               "GNC integrands and multipoles" => "GNC_Correlations_3.md",

               "LD correlations" => "LD_Correlations_1.md",
               "LD integrands and multipoles" => "LD_Correlations_2.md",

               "GNCxLD correlations" => "GNCxLD_Correlations_1.md",
               "GNCxLD integrands and multipoles" => "GNCxLD_Correlations_2.md",

               "LDxGNC correlations" => "LDxGNC_Correlations_1.md",
               "LDxGNC integrands and multipoles" => "LDxGNC_Correlations_2.md",

               #"GNC" => [
               #     "Auto-correlations" => "GNC_Correlations_1.md",
               #     "Cross-correlations" => "GNC_Correlations_2.md",
               #     "Integrands and multipoles" => "GNC_Correlations_3.md",
               #],
               #"LD" => [
               #     "Correlations" => "LD_Correlations_1.md",
               #     "Integrands and multipoles" => "LD_Correlations_2.md",
               #],
               #"GNCxLD" => [
               #     "Cross-correlations" => "GNCxLD_Correlations_1.md",
               #     "Integrands and multipoles" => "GNCxLD_Correlations_2.md",
               #],
               #"LDxGNC" => [
               #     "Cross-correlations" => "LDxGNC_Correlations_1.md",
               #     "Multipoles" => "LDxGNC_Correlations_2.md",
               #],

               "The Plain Parralel Approximation" => "PlaneParallelApprox.md",
          ],
          "Calculating Power Spectra" => "PowerSpectra.md",
          "Power Spectra for a generic window" => "PowerSpectraGenWin.md",
          "implication on PNG" => "PNG.md",
          "Utilities" => [
               "MySpline" => "Spline.md",
               "Kernels" => "Kernels.md",
               "Dictionaries and names" => "Dicts.md",
               "Mathematical Utilities" => "MathUtils.md",
               "Cosmology Utilities" => "CosmoUtils.md",
               "Input Power Spectrum Tools" => "IPSTools.md",
               "Other Utilities" => "OtherUtils.md",
          ],
     ],
)

deploydocs(
     repo = "github.com/foglienimatteo/GaPSE.jl.git",
     devbranch = "main",
     push_preview = true,
)
