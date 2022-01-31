using AeroMDAO
using Documenter

## Generate theme
using DocumenterTools: Themes
for w in ("light",)
    header = read(joinpath(@__DIR__, "theme/style.scss"), String)
    theme = read(joinpath(@__DIR__, "theme/$(w)defs.scss"), String)
    write(joinpath(@__DIR__, "theme/$(w).scss"), header*"\n"*theme)
end

Themes.compile(joinpath(@__DIR__, "theme/light.scss"), joinpath(@__DIR__, "src/assets/themes/documenter-light.css"))
Themes.compile(joinpath(@__DIR__, "theme/dark.scss"), joinpath(@__DIR__, "src/assets/themes/documenter-dark.css"))

## Generate Markdown files using Literate.jl
using Literate

src = joinpath(@__DIR__, "src")
lit = joinpath(@__DIR__, "lit")

for (root, _, files) ∈ walkdir(lit), file ∈ files
    splitext(file)[2] == ".jl" || continue
    ipath = joinpath(root, file)
    opath = splitdir(replace(ipath, lit => src))[1]
    Literate.markdown(ipath, opath)
end

## Generate documentation
makedocs(
    # modules = [AeroMDAO, AeroMDAO.VortexLattice],
    sitename = "AeroMDAO",
    authors  = "Arjit Seth, Stephane Redonnet, and Rhea P. Liem",
    # repo = "https://github.com/GodotMisogi/AeroMDAO.jl",
    pages = [
                "Home"          => "index.md"
                "Tutorials"     => [
                                    "Airfoil Aerodynamic Analysis" => "tutorials-airfoil.md",
                                    "Aircraft Aerodynamic Analysis" => "tutorials-aircraft.md"
                                    ]
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
                                "assets/logo.ico",
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