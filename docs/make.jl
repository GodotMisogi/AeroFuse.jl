using AeroMDAO
using Documenter
using DocumenterTools: Themes
using Literate

## Generate theme
for w in ("light", "dark")
    header = read(joinpath(@__DIR__, "theme/style.scss"), String)
    theme = read(joinpath(@__DIR__, "theme/$(w)defs.scss"), String)
    write(joinpath(@__DIR__, "theme/$(w).scss"), header*"\n"*theme)
end

Themes.compile(joinpath(@__DIR__, "theme/light.scss"), joinpath(@__DIR__, "src/assets/themes/documenter-light.css"))
Themes.compile(joinpath(@__DIR__, "theme/dark.scss"), joinpath(@__DIR__, "src/assets/themes/documenter-dark.css"))

## Generate Markdown files using Literate.jl
Literate.markdown("src/tutorials.jl")
Literate.markdown("src/howto.jl")
Literate.markdown("src/theory.jl")

## Generate documentation
makedocs(
    # modules = [AeroMDAO, AeroMDAO.VortexLattice],
    sitename = "AeroMDAO",
    authors  = "Arjit Seth, Stephane Redonnet, and Rhea P. Liem",
    # repo = "https://github.com/GodotMisogi/AeroMDAO.jl",
    pages = [
                "Home"          => "index.md"
                "Tutorials"     => "tutorials.md"
                "How-to Guide"  => "howto.md"
                "Theory"        => "theory.md"
                "Reference"     => [
                                    "Geometry API"      => "geometry.md"
                                    "Aerodynamics API"  => "aerodynamics.md"
                                    "Structures API"    => "structures.md"
                                    "In-Progress API"   => "development.md"
                                    ]
            ],
    format = Documenter.HTML(
                            # /prettyurls = CI,
                            assets = [
                                # "assets/logo.ico",
                                asset("https://fonts.googleapis.com/css?family=Montesserat|Fira+Code&display=swap", class=:css),
                                ],
                            # highlightjs = "theme/highlight.js",
                            ),
)

## Deployment
deploydocs(
    repo = "github.com/GodotMisogi/AeroMDAO.jl.git",
    devbranch = "develop",
    versions = ["stable" => "v^", "v#.#", "dev" => "dev"],
)