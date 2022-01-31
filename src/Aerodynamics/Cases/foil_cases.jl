## 2D cases
#==========================================================================================#

"""
    solve_case(foil :: Foil, uniform :: Uniform2D; sources = false, wake_length = 1e3, num_panels :: Integer = 60)

Evaluate a doublet-source case given a `Foil` with a `Uniform2D`, with optional named arguments to specify whether the source terms are non-zero, the length of the wake, and the number of panels for the analysis.
"""
function solve_case(foil :: Foil, freestream :: Uniform2D; viscous = false, sources = false, wake_length = 1e5, num_panels :: Integer = 60, num_wake = 15)
    panels = make_panels(foil, num_panels)
    solve_system(panels, freestream, num_wake, wake_length)
end