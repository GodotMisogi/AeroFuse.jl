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
import .MathTools: tupvector, fwdsum, fwddiff, fwddiv, weighted_vector, vectarray, slope, splitat, adj3, columns, extend_yz, reflect_mapper, cosine_interp, structtolist, inverse_rotation, rotation, affine_2D, Point2D, Point3D, x, y, z, reshape_array, midpair_map, partition, uniform_spacing, linear_spacing, cosine_spacing, sine_spacing

export tupvector, fwdsum, fwddiff, fwddiv, weighted_vector, vectarray, slope, splitat, adj3, columns, extend_yz, reflect_mapper, cosine_interp, structtolist, inverse_rotation, rotation, affine_2D, Point2D, Point3D, reshape_array, midpair_map, partition, uniform_spacing, linear_spacing, cosine_spacing, sine_spacing


## Non-dimensionalization
#==========================================================================================#

include("Tools/NonDimensional.jl")
import .NonDimensional: dynamic_pressure, force_coefficient, moment_coefficient, rate_coefficient, pressure_coefficient, aerodynamic_coefficients, print_coefficients, reynolds_number, print_derivatives, force, moment

export dynamic_pressure, force_coefficient, moment_coefficient, rate_coefficient, pressure_coefficient, print_coefficients, reynolds_number, print_derivatives, force, moment

## Panels
#===========================================================================#

include("Geometry/PanelGeometry/PanelGeometry.jl")
import .PanelGeometry: AbstractPanel, AbstractPanel2D, Panel2D, WakePanel2D, AbstractPanel3D, Panel3D, p1, p2, p3, p4, zero, collocation_point, panel_length, transform_panel, transform_panel_points, panel_angle, panel_tangent, panel_normal, panel_location, panel_area, panel_coords, transform, midpoint, panel_points, wake_panel, wake_panels, reverse_panel, panel_velocity, panel_scalar, trailing_edge_panel, get_surface_values, panel_vector, panel_dist, average_chord, average_width, wetted_area

export AbstractPanel, AbstractPanel2D, Panel2D, WakePanel2D, AbstractPanel3D, Panel3D, Point2D, collocation_point, p1, p2, p3, p4, transform, panel_normal, midpoint, panel_location, panel_tangent, panel_points, panel_dist, wake_panel, wake_panels, panel_area, reverse_panel, panel_length, transform_panel, panel_angle, panel_vector, panel_velocity, panel_scalar, trailing_edge_panel, get_surface_values, average_chord, average_width, wetted_area


## Wing geometry
#==========================================================================================#

include("Geometry/AircraftGeometry/AircraftGeometry.jl")
import .AircraftGeometry: Aircraft, Foil, kulfan_CST, naca4, camber_CST, paneller, read_foil, split_foil, foil_camthick, camthick_foil, camber_thickness, cosine_foil, camthick_to_CST, coords_to_CST, max_thickness_to_chord_ratio_location, Fuselage, projected_area, length, cosine_spacing, HalfWing, HalfWingSection, Wing, WingSection, mean_aerodynamic_chord, span, aspect_ratio, projected_area, taper_ratio, info, leading_edge, trailing_edge, chop_leading_edge, chop_trailing_edge, chop_wing, chop_sections, chop_coordinates, chop_spanwise_sections, chop_chords, chop_spans, wing_bounds, paneller, mesh_horseshoes, mesh_wing, mesh_cambers, make_panels, vlmesh_wing, max_tbyc_sweeps, mean_aerodynamic_center, panel_wing, wetted_area, number_of_spanwise_panels, spanwise_spacing, coordinates, chord_coordinates, camber_coordinates, surface_coordinates, foils, chords, twists, spans, dihedrals, sweeps, left, right

export Aircraft, Foil, kulfan_CST, naca4, camber_CST, paneller, read_foil, split_foil, foil_camthick, camthick_foil, camber_thickness, cosine_foil, camthick_to_CST, coords_to_CST, max_thickness_to_chord_ratio_location, # 2D setups
Fuselage, projected_area, length, cosine_spacing, # Fuselage
HalfWing, HalfWingSection, Wing, WingSection, mean_aerodynamic_chord, span, aspect_ratio, projected_area, taper_ratio, info, leading_edge, trailing_edge, chop_leading_edge, chop_trailing_edge, chop_wing, chop_sections, chop_coordinates, chop_spanwise_sections, chop_chords, chop_spans, wing_bounds, paneller, mesh_horseshoes, mesh_wing, mesh_cambers, make_panels, vlmesh_wing, max_tbyc_sweeps, mean_aerodynamic_center, panel_wing, wetted_area, number_of_spanwise_panels, spanwise_spacing, coordinates, chord_coordinates, camber_coordinates, surface_coordinates, foils, chords, twists, spans, dihedrals, sweeps, left, right

## Laplace
#==========================================================================================#

include("Tools/Laplace.jl")
import .Laplace: Uniform2D, Freestream, velocity, potential, stream, aircraft_velocity, cartesian_to_freestream, freestream_to_cartesian

export Uniform2D, velocity, Freestream, aircraft_velocity, stream, vortex_stream_1, vortex_stream_2, source_stream, cartesian_to_freestream, freestream_to_cartesian

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

## Vortex lattice

include("Aerodynamics/VortexLattice/VortexLattice.jl")
import .VortexLattice: Horseshoe, streamlines, solve_system, transform, bound_leg, bound_leg_center, bound_leg_vector, r1, r2, points, horseshoe_line, horseshoe_point, make_horseshoes, nearfield_drag, evaluate_coefficients, body_to_wind_axes, body_to_stability_axes, stability_to_body_axes, wind_to_body_axes, case_dynamics, VLMState, VLMSystem, VLMSurface, compute_influence_matrix!, compute_boundary_condition!, compute_horseshoes!, solve_system!, update_velocity!, evaluate_residual!, compute_surface_forces!, compute_surface_moments!, compute_farfield_forces!, generate_system!, update_circulations!, solve_case!, compute_wake_properties!, normals, horseshoes, collocation_points, AIC, RHS, circulations, surface_forces, surface_moments, surface_force_coefficients, surface_moment_coefficients, nearfield_coefficients, farfield_coefficients, aerodynamic_coefficients, build_system, rate_coefficient

export Horseshoe, streamlines, solve_system, transform, bound_leg, bound_leg_center, bound_leg_vector, r1, r2, points, horseshoe_line, horseshoe_point, make_horseshoes, nearfield_drag, evaluate_coefficients, body_to_wind_axes, body_to_stability_axes, stability_to_body_axes, wind_to_body_axes, case_dynamics,  # Pure
VLMState, VLMSystem, VLMSurface, compute_influence_matrix!, compute_boundary_condition!, compute_horseshoes!, solve_system!, update_velocity!, evaluate_residual!, compute_surface_forces!, compute_surface_moments!, compute_farfield_forces!, generate_system!, update_circulations!, solve_case!, compute_wake_properties!, normals, horseshoes, collocation_points, AIC, RHS, circulations, surface_forces, surface_moments, surface_force_coefficients, surface_moment_coefficients, nearfield_coefficients, farfield_coefficients, aerodynamic_coefficients, build_system, rate_coefficient # Impure

## Profile drag estimation

include("Aerodynamics/profile_drag.jl")

export wetted_area_drag, profile_drag_coefficient

## Cases

include("Aerodynamics/Cases/cases.jl")
include("Aerodynamics/Cases/stability_cases.jl")
include("Aerodynamics/Cases/foil_cases.jl")

export solve_case, solve_stability_case, streamlines, print_case, print_info

## Structural analyses
#==========================================================================================#

include("Structures/Beams.jl")
import .Beams: Material, Beam, Tube, radii, area, moment_of_inertia, polar_moment_of_inertia, J_coeffs, Iy_coeffs, Iz_coeffs, tube_stiffness_matrix, deflection_stiffness_matrix, torsional_stiffness_matrix

export Material, Beam, Tube, radii, area, moment_of_inertia, polar_moment_of_inertia, J_coeffs, Iy_coeffs, Iz_coeffs, tube_stiffness_matrix, deflection_stiffness_matrix, torsional_stiffness_matrix

## Post-processing
#==========================================================================================#

include("Tools/plot_tools.jl")

export plot_panels, plot_streams, plot_wing, plot_surface

end