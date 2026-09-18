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
               "Spline Theory" => "SplineTheory.md",
               "Spherical Bessel Functions" => "SphericalBesselFunctions.md",
               "The I_l^n integrals" => "IlnIntegrals.md",
               "The Δχ → 0 limits" => [
                    "Introduction and results" => "DeltaChiLimits.md",
                    "Family 1: Lensing x Lensing" => "DeltaChiLimits_1_LensingLensing.md",
                    "Family 2: Lensing x Doppler" => "DeltaChiLimits_2_LensingDoppler.md",
                    "Family 3: Newtonian x Lensing" => "DeltaChiLimits_3_NewtonianLensing.md",
                    "Family 4: Lensing x Local GP" => "DeltaChiLimits_4_LensingLocalGP.md",
                    "Family 5: Newtonian x Integrated GP" => "DeltaChiLimits_5_NewtonianIntegratedGP.md",
                    "Family 6: the Δχ⁴ Ĩ₀⁴ terms" => "DeltaChiLimits_6_Ichi4Tilde.md",
                    "Family 7: the vanishing-factor terms" => "DeltaChiLimits_7_VanishingFactor.md",
                    "Family 8: the J₂₂I₂² + J₃₁I₁³ terms" => "DeltaChiLimits_8_J22J31.md",
               ],
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
