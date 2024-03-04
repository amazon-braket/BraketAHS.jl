using BraketAHS
using Documenter

DocMeta.setdocmeta!(BraketAHS, :DocTestSetup, :(using BraketAHS); recursive=true)

makedocs(;
    modules=[BraketAHS],
    authors="Yaroslav Kharkov, Katharine Hyatt",
    sitename="BraketAHS.jl",
    format=Documenter.HTML(;
        canonical="https://ykharkov.github.io/BraketAHS.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/ykharkov/BraketAHS.jl",
    devbranch="main",
)
