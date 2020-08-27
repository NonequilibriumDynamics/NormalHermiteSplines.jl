using Documenter, NormalHermiteSplines

makedocs(
    sitename = "NormalHermiteSplines.jl",
	format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
	authors = "Igor Kohanovsky",
    pages = [
				"Home" => "index.md",
				"Public API" => "Public-API.md",
				"Example Usage" => "Usage.md",
				"Numerical Tests" => "Numerical-Tests.md",
				"Tests with real data" => "Tests-with-real-data.md",
				"Interpolating Normal Hermite Splines" => "Interpolating-Normal-Hermite-Splines.md",
				"Reproducing Kernel of Bessel Potential space" => "Reproducing-Kernel-of-Bessel-Potential-space.md",
				"Relation to Polyharmonic Splines"  => "Relation-to-Polyharmonic-Splines.md",
				"Updating Cholesky Factorization" => "Updating-Cholesky-Factorization.md",
			]
)

deploydocs(
    repo = "github.com/IgorKohan/NormalHermiteSplines.jl.git",
    devurl = "v0.3.0",
    versions = ["v0.3.0" => "v^", "v#.#"],
)
