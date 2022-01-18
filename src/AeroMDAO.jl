module AeroMDAO

## Libraries
#==========================================================================================#

using LinearAlgebra
using StaticArrays
using CoordinateTransformations
using Rotations
using ForwardDiff: jacobian!
using DiffResults: JacobianResult, jacobian, value
using PrettyTables
using StructArrays

using Statistics: mean

using SplitApplyCombine: combinedimsview, combinedims
export combinedimsview, combinedims

using ComponentArrays: ComponentVector, ComponentArray, valkeys
export ComponentVector, ComponentArray, valkeys

using Setfield
using LabelledArrays

## Methods to be extended in submodules
#==========================================================================================#

function velocity end

function solve_linear end

function properties end

## Math tools
#==========================================================================================#

include("Tools/MathTools.jl")
import .MathTools: fwdsum, fwddiff, fwddiv, weighted_vector, vectarray, slope, splitat, adj3, columns, extend_yz, reflect_mapper, cosine_interp, structtolist, inverse_rotation, rotation, affine_2D, Point2D, Point3D, x, y, z, reshape_array, midpair_map, partition, uniform_spacing, linear_spacing, cosine_spacing, sine_spacing

# export fwdsum, fwddiff, fwddiv, weighted_vector, vectarray, slope, splitat, adj3, columns, extend_yz, reflect_mapper, cosine_interp, structtolist, inverse_rotation, rotation, affine_2D, Point2D, Point3D, reshape_array, midpair_map, partition, uniform_spacing, linear_spacing, cosine_spacing, sine_spacing


## Non-dimensionalization
#==========================================================================================#

include("Tools/NonDimensional.jl")
import .NonDimensional: dynamic_pressure, force_coefficient, moment_coefficient, rate_coefficient, pressure_coefficient, aerodynamic_coefficients, reynolds_number, force, moment

export dynamic_pressure, force_coefficient, moment_coefficient, rate_coefficient, pressure_coefficient, reynolds_number, force, moment

## Panels
#==========================================================================================#

include("Geometry/PanelGeometry/PanelGeometry.jl")
import .PanelGeometry: AbstractPanel, AbstractPanel2D, Panel2D, WakePanel2D, AbstractPanel3D, Panel3D, p1, p2, p3, p4, zero, collocation_point, panel_length, transform_panel, transform_panel_points, panel_angle, panel_tangent, panel_normal, panel_location, panel_area, panel_coordinates, transform, midpoint, panel_points, wake_panel, wake_panels, reverse_panel, panel_velocity, panel_scalar, trailing_edge_panel, get_surface_values, panel_vector, panel_dist, average_chord, average_width, wetted_area, make_panels

export AbstractPanel, AbstractPanel2D, Panel2D, WakePanel2D, AbstractPanel3D, Panel3D, collocation_point, p1, p2, p3, p4, transform, panel_normal, midpoint, panel_location, panel_tangent, panel_points, panel_dist, wake_panel, wake_panels, panel_area, reverse_panel, panel_length, transform_panel, panel_angle, panel_vector, panel_velocity, panel_scalar, trailing_edge_panel, get_surface_values, average_chord, average_width, wetted_area, make_panels


## Wing geometry
#==========================================================================================#

include("Geometry/AircraftGeometry/AircraftGeometry.jl")
import .AircraftGeometry: AbstractAircraft, AbstractWing, AbstractFoil, Foil, arc_length, kulfan_CST, naca4, camber_CST, paneller, read_foil, split_foil, coordinates_to_camber_thickness, camber_thickness_to_coordinates, camber_thickness, camber_thickness_to_coordinates, cosine_foil, camber_thickness_to_CST, coordinates_to_CST, max_thickness_to_chord_ratio_location, Fuselage, projected_area, length, cosine_spacing, HalfWing, HalfWingSection, Wing, WingSection, affine_transformation, mean_aerodynamic_chord, span, aspect_ratio, projected_area, taper_ratio, leading_edge, trailing_edge, chop_leading_edge, chop_trailing_edge, chop_wing, chop_sections, chop_coordinates, chop_spanwise_sections, chop_chords, chop_spans, wing_bounds, paneller, mesh_horseshoes, mesh_wing, mesh_cambers, max_thickness_to_chord_ratio_sweeps, mean_aerodynamic_center, panel_wing, number_of_spanwise_panels, symmetric_spacing, coordinates, chord_coordinates, camber_coordinates, surface_coordinates, foils, chords, twists, spans, dihedrals, sweeps, position, orientation, WingMesh, chord_panels, camber_panels, normal_vectors, surface_panels, AbstractSpacing, Sine, Cosine, Uniform, properties

export AbstractAircraft, AbstractWing, AbstractFoil, Foil, arc_length, kulfan_CST, naca4, camber_CST, paneller, read_foil, split_foil, coordinates_to_camber_thickness, camber_thickness_to_coordinates, camber_thickness, camber_thickness_to_coordinates, cosine_foil, camber_thickness_to_CST, coordinates_to_CST, max_thickness_to_chord_ratio_location, # 2D setups
Fuselage, projected_area, length, cosine_spacing, # Fuselage
HalfWing, HalfWingSection, Wing, WingSection, affine_transformation, mean_aerodynamic_chord, span, aspect_ratio, projected_area, taper_ratio, leading_edge, trailing_edge, chop_leading_edge, chop_trailing_edge, chop_wing, chop_sections, chop_coordinates, chop_spanwise_sections, chop_chords, chop_spans, wing_bounds, paneller, mesh_horseshoes, mesh_wing, mesh_cambers, max_thickness_to_chord_ratio_sweeps, mean_aerodynamic_center, panel_wing, number_of_spanwise_panels, symmetric_spacing, coordinates, chord_coordinates, camber_coordinates, surface_coordinates, foils, chords, twists, spans, dihedrals, sweeps, position, orientation, WingMesh, chord_panels, camber_panels, normal_vectors, surface_panels, AbstractSpacing, Sine, Cosine, Uniform, properties

make_horseshoes(wing :: WingMesh) = StructArray(Horseshoe.(chord_panels(wing), normal_vectors(wing)))

make_vortex_rings(wing :: WingMesh) = StructArray(VortexRing.(camber_panels(wing)))

export make_horseshoes, make_vortex_rings

## Laplace
#==========================================================================================#

include("Tools/Laplace.jl")
import .Laplace: Uniform2D, potential, stream,  cartesian_to_freestream, freestream_to_cartesian

export Uniform2D, stream, vortex_stream_1, vortex_stream_2, source_stream, cartesian_to_freestream, freestream_to_cartesian

## Aerodynamic analyses
#==========================================================================================#

## Doublet-source panel method

include("Aerodynamics/DoubletSource/DoubletSource.jl")
import .DoubletSource: solve_problem, doublet_matrix, source_matrix, boundary_vector, wake_panels, source_strengths, tangential_velocities, solve_strengths, eval_coefficients, lift_coefficient

export doublet_matrix, source_matrix, boundary_vector, wake_panels, source_strengths, tangential_velocities, solve_strengths

## Linear-strength source and vorticity panel method

include("Aerodynamics/LinearVortexSource/LinearVortexSource.jl")
import .LinearVortexSource: total_velocity, source_velocity, vortex_velocity, vortex_influence_matrix, source_influence_matrix, neumann_boundary_condition, kutta_condition, two_point_matrix, linear_source_matrix, linear_vortex_matrix, constant_source_matrix, constant_source_boundary_condition

export total_velocity, source_velocity, vortex_velocity, vortex_influence_matrix, source_influence_matrix, neumann_boundary_condition, kutta_condition, two_point_matrix, linear_source_matrix, linear_vortex_matrix, constant_source_matrix, constant_source_boundary_condition

## Vortex lattice

include("Aerodynamics/VortexLattice/VortexLattice.jl")
import .VortexLattice: Horseshoe, VLMSystem, References, AbstractAxisSystem, Stability, Wind, Body, Geometry, streamlines, influence_coefficient, influence_matrix, boundary_condition, solve_system, transform, bound_leg, horseshoe_point, bound_leg_center, bound_leg_vector, r1, r2, points, Horseshoe, horseshoe_point, surface_velocity, surface_forces, surface_moments, nearfield_drag, geometry_to_wind_axes, geometry_to_stability_axes, stability_to_geometry_axes, wind_to_geometry_axes,  rate_coefficient, nearfield, farfield, farfield_forces, surface_velocities, surface_forces, surface_dynamics, surface_coefficients, nearfield_coefficients, farfield_coefficients, VortexRing, Freestream, velocity, body_frame_velocity, kinematic_viscosity, mach_number

export Horseshoe, VLMSystem, References, AbstractAxisSystem, Stability, Wind, Body, Geometry, streamlines, influence_coefficient, influence_matrix, boundary_condition, solve_system, transform, bound_leg, horseshoe_point, bound_leg_center, bound_leg_vector, r1, r2, points, Horseshoe, horseshoe_point, surface_velocity, surface_forces, surface_moments, nearfield_drag, geometry_to_wind_axes, geometry_to_stability_axes, stability_to_geometry_axes, wind_to_geometry_axes,  rate_coefficient, nearfield, farfield, farfield_forces, surface_velocities, surface_forces, surface_dynamics, surface_coefficients, nearfield_coefficients, farfield_coefficients, VortexRing, Freestream, velocity, body_frame_velocity, kinematic_viscosity, mach_number

## Profile drag estimation

include("Aerodynamics/profile_drag.jl")

export wetted_area_drag, profile_drag_coefficient, local_dissipation_drag

## Viscous airfoil analysis

include("Aerodynamics/ViscFoil/ViscFoil.jl")

import .ViscFoil: BoundaryLayer2D, solve_inviscid_doublets, solve_inviscid_vortices, defect_block, edge_velocities, solve_viscous_case

export BoundaryLayer2D, solve_inviscid_doublets, solve_inviscid_vortices, defect_block, edge_velocities, solve_viscous_case

## Cases

include("Aerodynamics/Cases/cases.jl")
include("Aerodynamics/Cases/stability_cases.jl")
include("Aerodynamics/Cases/foil_cases.jl")

export solve_case, solve_case_derivatives, streamlines, print_case, print_info, print_coefficients, print_derivatives, span_loads,
triangle_connectivities, extrapolate_point_mesh

## Structural analyses
#==========================================================================================#

include("Structures/Beams.jl")
import .Beams: Material, Tube, Beam, radii, area, moment_of_inertia, polar_moment_of_inertia, J_coeffs, Iyy_coeffs, Izz_coeffs, tube_stiffness_matrix, bending_stiffness_matrix, axial_stiffness_matrix, build_stiffness_matrix, solve_cantilever_beam, elastic_modulus, shear_modulus, yield_stress, density, principal_stress, torsional_stress, von_mises_stress

export Material, Tube, Beam, radii, area, moment_of_inertia, polar_moment_of_inertia, J_coeffs, Iyy_coeffs, Izz_coeffs, tube_stiffness_matrix, bending_stiffness_matrix, axial_stiffness_matrix, build_stiffness_matrix, solve_cantilever_beam, elastic_modulus, shear_modulus, yield_stress, density, principal_stress, torsional_stress, von_mises_stress

## Aerostructural analyses
#==========================================================================================#

include("Aerostructural/Aerostructural.jl")
import .Aerostructural: AerostructWing, make_beam_mesh, transform_stiffy, permute_stiffy, build_big_stiffy, adjacent_adder, section_moments, compute_loads, fem_load_vector, rotation_matrix, transfer_displacements, mesh_translation, mesh_rotation, new_horseshoes, solve_coupled_residual!, aerostruct_gauss_seidel

export AerostructWing, make_beam_mesh, transform_stiffy, permute_stiffy, build_big_stiffy, adjacent_adder, section_moments, compute_loads, fem_load_vector, rotation_matrix, transfer_displacements, mesh_translation, mesh_rotation, new_horseshoes, solve_coupled_residual!, aerostruct_gauss_seidel

## Post-processing
#==========================================================================================#

include("Tools/plot_tools.jl")

export plot_panels, plot_streamlines, plot_wing, plot_surface

end