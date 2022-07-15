using Pertussis
using Documenter

DocMeta.setdocmeta!(Pertussis, :DocTestSetup, :(using Pertussis); recursive=true)

makedocs(;
    modules=[Pertussis],
    authors="Pengfei Song",
    repo="https://github.com/Song921012/Pertussis.jl/blob/{commit}{path}#{line}",
    sitename="Pertussis.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://Song921012.github.io/Pertussis.jl",
        edit_link="master",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/Song921012/Pertussis.jl",
    devbranch="master",
)
