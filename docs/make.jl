using Documenter, SpinWeightedSpheroidalHarmonics

makedocs(
    sitename="SpinWeightedSpheroidalHarmonics.jl",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true"
    )
)