using PadeDiagonals
using Documenter

DocMeta.setdocmeta!(PadeDiagonals, :DocTestSetup, :(using PadeDiagonals); recursive=true)

makedocs(;
    modules=[PadeDiagonals],
    authors="Joris Pinkse <pinkse@gmail.com> and contributors",
    sitename="PadeDiagonals.jl",
    format=Documenter.HTML(;
        canonical="https://NittanyLion.github.io/PadeDiagonals.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/NittanyLion/PadeDiagonals.jl",
    devbranch="main",
)
