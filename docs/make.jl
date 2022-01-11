using Documenter, AeroMDAO
using DocumenterTools: Themes

##
for w in ("light", "dark")
    header = read(joinpath(@__DIR__, "theme/style.scss"), String)
    theme = read(joinpath(@__DIR__, "theme/$(w)defs.scss"), String)
    write(joinpath(@__DIR__, "theme/$(w).scss"), header*"\n"*theme)
end
using DocumenterTools: Themes
Themes.compile(joinpath(@__DIR__, "theme/light.scss"), joinpath(@__DIR__, "src/assets/themes/documenter-light.css"))
Themes.compile(joinpath(@__DIR__, "theme/dark.scss"), joinpath(@__DIR__, "src/assets/themes/documenter-dark.css"))

## Generate documentation
makedocs(
    # modules = [AeroMDAO, AeroMDAO.VortexLattice],
    sitename = "AeroMDAO",
    authors  = "Arjit Seth, Rhea P. Liem, and Stephane Redonnet",
    # repo = "https://github.com/GodotMisogi/AeroMDAO.jl",
    pages = [
                "Home"   => "index.md",
                "Guide"  => [ 
                                "Geometry"      => "geometry.md"
                                "Aerodynamics"  => "aerodynamics.md"
                                "Structures"    => "structures.md"
                                "Theory"        => "theory.md"
                            ],
                "API"    => "api.md"
            ],
    format = Documenter.HTML(
                            # /prettyurls = CI,
                            assets = [
                                # "assets/logo.ico",
                                asset("https://fonts.googleapis.com/css?family=Montesserat|Fira+Code&display=swap", class=:css),
                                ],
                            # highlightjs = "theme/highlight.js",
                            )
)

## Deployment
deploydocs(
    # root = "<current-directory>",
    target = "build",
    # dirname = "",
    repo = "github.com/GodotMisogi/AeroMDAO.jl.git",
    branch = "gh-pages",
    # deps = nothing | <Function>,
    # make = nothing | <Function>,
    devbranch = "develop",
    devurl = "dev",
    versions = ["stable" => "v^", "v#.#", devurl => devurl],
    # forcepush = false,
    # deploy_config = auto_detect_deploy_system(),
    # push_preview = false,
    # repo_previews = repo,
    # branch_previews = branch,
)