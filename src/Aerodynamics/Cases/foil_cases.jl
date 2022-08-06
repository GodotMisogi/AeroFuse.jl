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

function solve_case(wing :: Wing, α, β, U, wake_length = 1.0e5, npt_span :: Vector{Integer} = [40], npt_chord :: Integer = 60)
	wing_mesh = WingMesh(wing, npt_span, npt_chord)
	surf_pts  = surface_coordinates(wing_mesh)
	surf_pans = make_panels(surf_pts)
	foil_pans = @view permutedims(surf_pans)[:]

	fs = Freestream(α, β, zeros(3))
	# V∞ = U * velocity(fs)

	solve_system(foil_pans, U, fs, wake_length)
end