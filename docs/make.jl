push!(LOAD_PATH,"../src/")

using Documenter
using GaPSE

Documenter.makedocs(
	format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
	modules = [GaPSE],
	sitename = "GaPSE.jl",
	pages = [
			"Introduction" => "index.md",
			],
)

deploydocs(repo = "github.com/cosmofico97/GaPSE.git")
