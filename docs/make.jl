using DiscreteIntegration
using Documenter

DocMeta.setdocmeta!(DiscreteIntegration, :DocTestSetup, :(using DiscreteIntegration); recursive=true)

makedocs(;
    modules=[DiscreteIntegration],
    authors="Kai-Jonathan Kuhn <kai-jonathan.kuhn@uni-bayreuth.de> and contributors",
    sitename="DiscreteIntegration.jl",
    format=Documenter.HTML(;
        canonical="https://kaikuhn.github.io/DiscreteIntegration.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/kaikuhn/DiscreteIntegration.jl",
    devbranch="main",
)
