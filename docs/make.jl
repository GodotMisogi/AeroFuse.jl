using AeroFuse
using Documenter
using Literate

# NOTE: The old custom Documenter theme (theme/*.scss) targeted Documenter 0.27's Bulma SCSS
# variables and no longer compiles under Documenter 1.x, so the build uses Documenter's
# built-in light/dark themes. Re-port the SCSS to the 1.x variables to restore it.

## Generate Markdown files using Literate.jl
src = joinpath(@__DIR__, "src")
lit = joinpath(@__DIR__, "lit")

for (root, _, files) ∈ walkdir(lit), file ∈ files
    splitext(file)[2] == ".jl" || continue
    ipath = joinpath(root, file)
    opath = splitdir(replace(ipath, lit => src))[1]
    Literate.markdown(ipath, opath)
    # Literate.notebook(ipath, opath)
end

## Generate documentation
makedocs(
    # modules = [AeroFuse, AeroFuse.VortexLattice],
    sitename = "AeroFuse",
    authors  = "Arjit Seth and Rhea P. Liem",
    # repo = "https://github.com/GodotMisogi/AeroFuse.jl",
    pages = [
        "Home"          => "index.md"
        "Tutorials"     =>  [
            "Airfoil Aerodynamic Analysis"  => "tutorials-airfoil.md",
            "Aircraft Aerodynamic Analysis" => "tutorials-aircraft.md",
            "Aerodynamic Stability Analysis" => "tutorials-stability.md",
        ]
        "How-to Guide"  => "howto.md"
        "Theory"        => "theory.md"
        "Reference"     =>  [
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
    checkdocs = :exports,
    # Documenter 1.x makes @example failures fatal by default. Some tutorials download airfoil
    # data over HTTP, so a flaky network would abort the whole build; keep those non-fatal
    # (0.27 behaviour) while leaving doctests, cross-references and docstring checks strict.
    warnonly = [:example_block],
)

## Deployment
deploydocs(
    repo = "github.com/GodotMisogi/AeroFuse.jl.git",
    devbranch = "develop",
    versions = ["stable" => "v^", "v#.#", "dev" => "dev"],
)