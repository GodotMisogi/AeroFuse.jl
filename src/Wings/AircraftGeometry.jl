#-------------------------AIRCRAFT GUFF---------------------#

abstract type Aircraft end

## Foil geometry
#==========================================================================================#

include("foil_geometry.jl")

export Foil, read_foil, kulfan_CST, naca4

## Wing geometry
#==========================================================================================#

include("wing_geometry.jl")

export HalfWing, Wing, mean_aerodynamic_chord, span, aspect_ratio, projected_area, taper_ratio, info, print_info, lead_wing, wing_bounds, paneller, mesh_horseshoes, mesh_wing, mesh_cambers, make_panels, vlmesh_wing