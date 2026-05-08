using Nbody2D
using Documenter

DocMeta.setdocmeta!(Nbody2D, :DocTestSetup, :(using Nbody2D); recursive=true)

makedocs(;
    modules=[Nbody2D],
    authors="Job Feldbrugge <14946916+jfeldbrugge@users.noreply.github.com> and contributors",
    sitename="Nbody2D.jl",
    format=Documenter.HTML(;
        canonical="https://jfeldbrugge.github.io/Nbody2D.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Theory" => "theory.md",
        "Tutorial" => "tutorial.md",
        "Citations" => "citations.md",
        "Reference" => "reference.md",
    ],
)

deploydocs(;
    repo="github.com/jfeldbrugge/Nbody2D.jl",
    devbranch="main",
)
