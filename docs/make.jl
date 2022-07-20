using Pertussis
using Documenter

DocMeta.setdocmeta!(Pertussis, :DocTestSetup, :(using Pertussis); recursive=true)

makedocs(;
    modules=[Pertussis],
    authors="Pengfei Song",
    repo="https://github.com/LIAM-COVID-19-Forecasting/Pertussis.jl/blob/{commit}{path}#{line}",
    sitename="Pertussis.jl",
    format = Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://LIAM-COVID-19-Forecasting.github.io/Pertussis.jl",
        edit_link="master",
        assets=String[]
    ),
    pages=[
        "Home" => "index.md",
    ]
)

deploydocs(;
    repo="github.com/LIAM-COVID-19-Forecasting/Pertussis.jl",
    devbranch="master"
)
