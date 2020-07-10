using Documenter, NormalHermiteSplines

makedocs(
    sitename = "NormalHermiteSplines.jl",
	format = Documenter.HTML(),
	authors = "Igor Kohanovsky",	
    pages = [
				"Home" => "index.md",
				"Example Usage" => "Usage.md",
				"Public API" => "Public-API.md",
				"Interpolating Normal Hermite Splines" => "Interpolating-Normal-Hermite-Splines.md",
				"Reproducing Kernel of Bessel Potential space" => "Reproducing-Kernel-of-Bessel-Potential-space.md",
			]
)

deploydocs(
	       repo = "github.com/IgorKohan/NormalHermiteSplines.jl.git",
           devurl = "v0.1.0",
           versions = ["v0.1.0" => "v^", "v#.#"],
	      )
