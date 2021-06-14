using Documenter
using BeamPropagationMethod

makedocs(    
    modules = [BeamPropagationMethod],
    authors="Oliver Kliebisch <oliver@kliebisch.net> and contributors",
    repo="https://github.com/tachawkes/BeamPropagationMethod.jl/blob/{commit}{path}#L{line}",
    sitename="BeamPropagationMethod.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        assets=String[],
    ),
    pages=[
        "Introduction" => "index.md"        
        "Examples" => [
            "Simple multimode fiber" => "examples/example1.md"
        ]
    ],
)

deploydocs(
    repo = "github.com/TacHawkes/BeamPropagationMethod.jl.git"    
)
