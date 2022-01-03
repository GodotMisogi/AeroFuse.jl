using Documenter, AeroMDAO

makedocs(
    # modules = [AeroMDAO, AeroMDAO.VortexLattice],
    sitename = "AeroMDAO",
    # repo = "https://github.com/GodotMisogi/AeroMDAO.jl",
    pages = [
                "Home"   => "index.md",
                "Guide" => [ 
                                # "Guide"         => "guide.md"
                                "Aerodynamics"  => "aerodynamics.md"
                                "Structures"    => "structures.md"
                                "Theory"        => "theory.md"
                            ],
                "API"    => "api.md"
            ],
)