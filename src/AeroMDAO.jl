module AeroMDAO

#----------------------IMPORTS--------------------------------#
using LinearAlgebra
using StaticArrays
using CoordinateTransformations, Rotations
using ForwardDiff, DiffResults
using PrettyTables

## Math tools
#==========================================================================================#

include("Tools/MathTools.jl")
import .MathTools: tupvector, fwdsum, fwddiff, fwddiv, weighted_vector, vectarray, slope, splitat, adj3, columns, extend_yz, reflect_mapper, cosine_dist, sine_dist, cosine_interp, structtolist, inverse_rotation, rotation, affine_2D, Point2D, Point3D, x, y, z, reshape_array, midpair_map, partition

export Point2D, Point3D, tupvector, tuparray, midpair_map, cosine_dist, cosine_interp, sine_dist, rotation, inverse_rotation, affine_2D, fwdsum, fwddiff, fwddiv, partition, weighted_vector


## Non-dimensionalization
#==========================================================================================#

include("Tools/NonDimensional.jl")
using .NonDimensional

export dynamic_pressure, force_coefficient, moment_coefficient, rate_coefficient, pressure_coefficient, aerodynamic_coefficients, print_coefficients, reynolds_number, print_derivatives

## Panels
#===========================================================================#

include("Geometry/PanelGeometry.jl")
using .PanelGeometry

export Panel, Panel2D, WakePanel2D, Panel3D, Point2D, collocation_point, point1, point2, point3, point4, transform, panel_normal, midpoint, panel_location, panel_tangent, panel_points, panel_dist, wake_panel, wake_panels, panel_area, reverse_panel, panel_length, transform_panel, panel_angle, panel_vector, panel_velocity, panel_scalar, trailing_edge_panel, get_surface_values


## Wing geometry
#==========================================================================================#

include("Geometry/AircraftGeometry.jl")
using .AircraftGeometry

export Aircraft, Foil, kulfan_CST, naca4, camber_CST, paneller, read_foil, split_foil, foil_camthick, camthick_foil, camber_thickness, cosine_foil, camthick_to_CST, coords_to_CST, max_thickness_to_chord_ratio_location, # 2D setups
Fuselage, projected_area, length, cosine_distribution, # Fuselage
HalfWing, Wing, mean_aerodynamic_chord, span, aspect_ratio, projected_area, taper_ratio, info, print_info, leading_edge, leading_chopper, trailing_chopper, wing_chopper, wing_bounds, paneller, mesh_horseshoes, mesh_wing, mesh_cambers, make_panels, vlmesh_wing, max_tbyc_sweeps, mean_aerodynamic_center, panel_wing, wetted_area

## Laplace
#==========================================================================================#

include("Tools/Laplace.jl")
import .Laplace: Uniform2D, Freestream, velocity, potential, stream, aircraft_velocity, cartesian_to_freestream

export Uniform2D, velocity, Freestream, aircraft_velocity, stream, vortex_stream_1, vortex_stream_2, source_stream, cartesian_to_freestream

## Aerodynamic analyses
#==========================================================================================#

## Doublet-source panel method

include("Aerodynamics/DoubletSource/DoubletSource.jl")
import .DoubletSource: solve_problem, doublet_matrix, source_matrix, boundary_vector, wake_panels, source_strengths, tangential_velocities, solve_strengths, eval_coefficients, lift_coefficient

export doublet_matrix, source_matrix, boundary_vector, wake_panels, source_strengths, tangential_velocities, solve_strengths

## Linear-strength source and vorticity panel method

include("Aerodynamics/LinearVortexSource/LinearVortexSource.jl")
import .LinearVortexSource: total_velocity, source_velocity, vortex_velocity, vortex_influence_matrix, source_influence_matrix, neumann_boundary_condition, kutta, two_point_matrix, linear_source_matrix, linear_vortex_matrix, constant_source_matrix, constant_source_boundary_condition

export total_velocity, source_velocity, vortex_velocity, vortex_influence_matrix, source_influence_matrix, neumann_boundary_condition, kutta, two_point_matrix, linear_source_matrix, linear_vortex_matrix, constant_source_matrix, constant_source_boundary_condition

## Vorticity-streamfunction panel method

# include("Aerodynamics/VorticityStream/VorticityStream.jl")
# import .VorticityStream: influence_matrix, boundary

# export influence_matrix, boundary

## Viscous-inviscid analysis
# include("Aerodynamics/ViscFoil/ViscFoil.jl")
# import .ViscFoil: solve_bl_case, defect_block, edge_velocities, solve_inviscid_vortices, solve_inviscid_doublets

# export solve_bl_case, defect_block, edge_velocities, solve_inviscid_vortices, solve_inviscid_doublets

## Vortex lattice

include("Aerodynamics/VortexLattice/VortexLattice.jl")
using .VortexLattice

export Horseshoe, streamlines, solve_horseshoes, transform, horseshoe_line, horseshoe_point, nearfield_drag, evaluate_coefficients, body_to_stability_axes, stability_to_body_axes, body_to_wind_axes, wind_to_body_axes

## Profile drag estimation

include("Aerodynamics/profile_drag.jl")

export wetted_area_drag

## Cases

include("Aerodynamics/Cases/cases.jl")
include("Aerodynamics/Cases/stability_cases.jl")
include("Aerodynamics/Cases/foil_cases.jl")

export solve_case, solve_stability_case, streamlines, print_case

## Post-processing
#==========================================================================================#

include("Tools/plot_tools.jl")

export plot_panels, plot_streams, plot_wing, plot_surface

end