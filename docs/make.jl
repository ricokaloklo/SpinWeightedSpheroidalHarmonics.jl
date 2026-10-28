using Pkg; Pkg.add("Documenter")
using Documenter, SpinWeightedSpheroidalHarmonics

makedocs(
    sitename="SpinWeightedSpheroidalHarmonics.jl",
    pages = [
        "Home" => "index.md",
        "Arbitrary precision" => "precision.md",
        "Examples" => "examples.md",
        "APIs" => "APIs.md",
    ],
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true"
    )
)

deploydocs(
    repo = "github.com/ricokaloklo/SpinWeightedSpheroidalHarmonics.jl.git",
)
