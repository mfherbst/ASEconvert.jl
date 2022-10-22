using ASEconvert
using Documenter

DocMeta.setdocmeta!(ASEconvert, :DocTestSetup, :(using ASEconvert); recursive=true)

makedocs(;
    modules=[ASEconvert],
    authors="Michael F. Herbst <info@michael-herbst.com> and contributors",
    repo="https://github.com/mfherbst/ASEconvert.jl/blob/{commit}{path}#{line}",
    sitename="ASEconvert.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://mfherbst.github.io/ASEconvert.jl",
        edit_link="master",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/mfherbst/ASEconvert.jl",
    devbranch="master",
)
