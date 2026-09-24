push!(LOAD_PATH, "../src/")

using Documenter
using Documenter.JSON
using GaPSE

Documenter.makedocs(
     format=Documenter.HTML(prettyurls=get(ENV, "CI", nothing) == "true"),
     modules=[GaPSE],
     sitename="GaPSE.jl",
     pages=[
          "Introduction" => "index.md",
          "Theory" => [
               "Spline Theory" => "theory_SplineTheory.md",
               "Spherical Bessel Functions" => "theory_SphericalBesselFunctions.md",
               "The input Power Spectrum" => "theory_InputPowerSpectrum.md",
               "The I_l^n integrals" => "theory_IlnIntegrals.md",
               "The Δχ → 0 limits" => [
                    "Introduction and results" => "theory_DeltaChiLimits.md",
                    "Family 1: Lensing x Lensing" => "theory_DeltaChiLimits_1_LensingLensing.md",
                    "Family 2: Lensing x Doppler" => "theory_DeltaChiLimits_2_LensingDoppler.md",
                    "Family 3: Newtonian x Lensing" => "theory_DeltaChiLimits_3_NewtonianLensing.md",
                    "Family 4: Lensing x Local GP" => "theory_DeltaChiLimits_4_LensingLocalGP.md",
                    "Family 5: Newtonian x Integrated GP" => "theory_DeltaChiLimits_5_NewtonianIntegratedGP.md",
                    "Family 6: the Δχ⁴ Ĩ₀⁴ terms" => "theory_DeltaChiLimits_6_Ichi4Tilde.md",
                    "Family 7: the vanishing-factor terms" => "theory_DeltaChiLimits_7_VanishingFactor.md",
                    "Family 8: the J₂₂I₂² + J₃₁I₁³ terms" => "theory_DeltaChiLimits_8_J22J31.md",
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

# Pages that live in `docs/src/` but are absent from the `pages` list above are still
# built and deployed by Documenter; they simply get no entry in the navigation menu.
# That is how the draft pages below are published: reachable by URL, invisible otherwise.
# Documenter offers no per-page switch to keep them out of the site search, so their
# records are dropped from `search_index.js` here, after `makedocs` has written it.
const UNLISTED_PAGES = ["theory_Ilnintegrals-mellin"]

let file = joinpath(@__DIR__, "build", "search_index.js")
     if isfile(file) && !isempty(UNLISTED_PAGES)
          text = read(file, String)
          # the file is `var documenterSearchIndex = {"docs": [ ... ]\n}`
          m = match(r"^(var documenterSearchIndex = \{\"docs\":\n)(.*)(\n\}\n?)$"s, text)
          if isnothing(m)
               @warn "make.jl: unexpected search_index.js layout, leaving it untouched." file
          else
               records = JSON.parse(m[2])
               kept = filter(r -> !any(p -> startswith(r["location"], p * "."), UNLISTED_PAGES), records)
               write(file, m[1], JSON.json(kept), m[3])
               @info "make.jl: dropped $(length(records) - length(kept)) search records" UNLISTED_PAGES
          end
     end
end

deploydocs(
     repo = "github.com/foglienimatteo/GaPSE.jl.git",
     devbranch = "main",
     push_preview = true,
)
