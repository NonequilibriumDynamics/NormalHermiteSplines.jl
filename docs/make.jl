using Documenter, NormalHermiteSplines

makedocs(
    sitename = "NormalHermiteSplines.jl",
	format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
	authors = "Igor Kohanovsky",
    pages = [
				"Home" => "index.md",
				"Example Usage" => "Usage.md",
				"Numerical Tests" => "Numerical-Tests.md",
				"Public API" => "Public-API.md",
				"Interpolating Normal Hermite Splines" => "Interpolating-Normal-Hermite-Splines.md",
				"Relation to Polyharmonic Splines"  => "Relation-to-Polyharmonic-Splines.md",
				"Reproducing Kernel of Bessel Potential space" => "Reproducing-Kernel-of-Bessel-Potential-space.md",
			]
)

deploydocs(
    repo = "github.com/IgorKohan/NormalHermiteSplines.jl.git",
    devurl = "v0.2.0",
    versions = ["v0.2.0" => "v^", "v#.#"],
)