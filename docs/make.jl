using Pkg; Pkg.add("Documenter")
using Documenter, SpinWeightedSpheroidalHarmonics

makedocs(
    sitename="SpinWeightedSpheroidalHarmonics.jl",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true"
    )
)

deploydocs(
    repo = "github.com/ricokaloklo/SpinWeightedSpheroidalHarmonics.jl.git",
)
