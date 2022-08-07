module AeroMDAO

## Libraries
#==========================================================================================#

using LinearAlgebra
using StaticArrays
using CoordinateTransformations
using Rotations
using PrettyTables
using RecipesBase
# using ReusePatterns

using Statistics: mean
using StatsBase: weights

using SplitApplyCombine: combinedimsview, combinedims, splitdimsview, splitdims
export combinedimsview, combinedims, splitdimsview, splitdims

using ComponentArrays

export ComponentVector, ComponentArray, valkeys

using Setfield
using LabelledArrays

## Methods to be extended in submodules
#==========================================================================================#
function velocity end

function solve_linear end

function solve_nonlinear end

function solve_linear! end

function solve_nonlinear! end

function solve_system end

function properties end

function surface_velocities end

function surface_coefficients end

export solve_linear, solve_nonlinear, solve_linear!, solve_nonlinear!, reynolds_number

## Math tools
#==========================================================================================#

include("Tools/MathTools.jl")
import .MathTools: forward_sum, forward_difference, forward_division, weighted_vector, vectarray, slope, splitat, adj3, columns, extend_yz, reflect_mapper, cosine_interp, structtolist, inverse_rotation, rotation, affine_2D, Point2D, Point3D, x, y, z, reshape_array, midpair_map, partition, uniform_spacing, linear_spacing, cosine_spacing, sine_spacing

# export forward_sum, forward_difference, forward_division, weighted_vector, vectarray, slope, splitat, adj3, columns, extend_yz, reflect_mapper, cosine_interp, structtolist, inverse_rotation, rotation, affine_2D, Point2D, Point3D, reshape_array, midpair_map, partition, uniform_spacing, linear_spacing, cosine_interpolation, sine_spacing

## Non-dimensionalization
#==========================================================================================#

include("Tools/NonDimensional.jl")
import .NonDimensional: dynamic_pressure, reynolds_number, force_coefficient, moment_coefficient, rate_coefficient, pressure_coefficient, aerodynamic_coefficients, force, moment

export dynamic_pressure, reynolds_number, force_coefficient, moment_coefficient, rate_coefficient, pressure_coefficient, force, moment

## Laplace
#==========================================================================================#

include("Tools/Laplace.jl")
import .Laplace: Freestream, Uniform2D, potential, stream, cartesian_to_freestream, freestream_to_cartesian

export Uniform2D, stream, vortex_stream_1, vortex_stream_2, source_stream, cartesian_to_freestream, freestream_to_cartesian

import .Laplace: DoubletLine3D, velocity

export DoubletLine3D, velocity

## Panels
#==========================================================================================#

include("Geometry/PanelGeometry/PanelGeometry.jl")
import .PanelGeometry: AbstractPanel, AbstractPanel2D, Panel2D, WakePanel2D, AbstractPanel3D, Panel3D, panel_length, transform_panel, transform_panel_points, panel_angle, tangent_vector, normal_vector, panel_location, panel_area, panel_coordinates, transform, midpoint, panel_points, wake_panel, wake_panels, reverse_panel, panel_velocity, panel_scalar, trailing_edge_panel, get_surface_values, panel_vector, distance, average_chord, average_width, wetted_area, make_panels, local_coordinate_system, trailing_edge_info

export AbstractPanel, AbstractPanel2D, Panel2D, WakePanel2D, AbstractPanel3D, Panel3D, transform, normal_vector, midpoint, panel_location, tangent_vector, panel_points, distance, wake_panel, wake_panels, panel_area, reverse_panel, panel_length, transform_panel, panel_angle, panel_vector, panel_velocity, panel_scalar, trailing_edge_panel, get_surface_values, average_chord, average_width, wetted_area, make_panels, local_coordinate_system, trailing_edge_info

export AbstractPanel, AbstractPanel2D, Panel2D, WakePanel2D, AbstractPanel3D, Panel3D, WakePanel3D, p1, p2, p3, p4, xs, ys, zs, midpoint, panel_location, panel_tangent, panel_points, distance, wake_panel, wake_panels, panel_area, reverse_panel, panel_length, transform_panel, panel_angle, panel_vector, panel_velocity, panel_scalar, trailing_edge_panel, get_surface_values, average_chord, average_width, wetted_area, make_panels, local_coordinate_system, panel_coordinates, collocation_point

## Aircraft geometry
#==========================================================================================#

include("Geometry/AircraftGeometry/AircraftGeometry.jl")

# Abstract types
import .AircraftGeometry: AbstractAircraft, AbstractWing, AbstractFoil

export AbstractAircraft, AbstractWing, AbstractFoil

# Foil
import .AircraftGeometry: Foil, arc_length, kulfan_CST, naca4, camber_CST, make_panels, read_foil, leading_edge_index, upper_surface, lower_surface, split_surface, coordinates_to_camber_thickness, camber_thickness_to_coordinates, camber_thickness, camber_thickness_to_coordinates, cosine_interpolation, camber_thickness_to_CST, coordinates_to_CST, maximum_thickness_to_chord, translate, interpolate, rotate, affine, scale, reflect, camber, camber_line, control_surface

export Foil, arc_length, kulfan_CST, naca4, camber_CST, make_panels, read_foil, leading_edge_index, upper_surface, lower_surface, split_surface, coordinates_to_camber_thickness, camber_thickness_to_coordinates, camber_thickness, camber_thickness_to_coordinates, cosine_interpolation, camber_thickness_to_CST, coordinates_to_CST, maximum_thickness_to_chord, translate, interpolate, rotate, affine, scale, reflect, camber, camber_line, control_surface

# Fuselage
import .AircraftGeometry: Fuselage, projected_area, length, cosine_interpolation, volume

export Fuselage, projected_area, length, cosine_interpolation, volume

# Wing
import .AircraftGeometry: Wing, WingSection, affine_transformation, mean_aerodynamic_chord, span, aspect_ratio, projected_area, taper_ratio, leading_edge, trailing_edge, chop_leading_edge, chop_trailing_edge, chop_wing, chop_sections, chop_coordinates, chop_spanwise_sections, chop_chords, chop_spans, wing_bounds, make_panels, mesh_chords, mesh_wing, mesh_cambers, max_thickness_to_chord_ratio_sweeps, mean_aerodynamic_center, panel_wing, number_of_spanwise_panels, symmetric_spacing, coordinates, chord_coordinates, camber_coordinates, surface_coordinates, foils, chords, twists, spans, dihedrals, sweeps, position, orientation, WingMesh, chord_panels, camber_panels, surface_panels, AbstractSpacing, Sine, Cosine, Uniform, properties, wetted_area_ratio

# export Wing, WingSection

export Wing, WingSection, affine_transformation, mean_aerodynamic_chord, span, aspect_ratio, projected_area, taper_ratio, leading_edge, trailing_edge, chop_leading_edge, chop_trailing_edge, chop_wing, chop_sections, chop_coordinates, chop_spanwise_sections, chop_chords, chop_spans, wing_bounds, make_panels, mesh_chords, mesh_wing, mesh_cambers, max_thickness_to_chord_ratio_sweeps, mean_aerodynamic_center, panel_wing, number_of_spanwise_panels, symmetric_spacing, coordinates, chord_coordinates, camber_coordinates, surface_coordinates, foils, chords, twists, spans, dihedrals, sweeps, position, orientation, WingMesh, chord_panels, camber_panels, surface_panels, AbstractSpacing, Sine, Cosine, Uniform, properties, wetted_area_ratio

# Surfaces
# import .AircraftGeometry: HorizontalTail, VerticalTail

# export HorizontalTail, VerticalTail

# Controls
# import .AircraftGeometry: WingControlSurface

# export WingControlSurface

make_horseshoes(wing :: WingMesh) = Horseshoe.(chord_panels(wing), normal_vector.(camber_panels(wing)))

make_vortex_rings(wing :: WingMesh) = VortexRing.(camber_panels(wing))

export make_horseshoes, make_vortex_rings

## Aerodynamic analyses
#==========================================================================================#

## Doublet-source panel method

include("Aerodynamics/DoubletSource/DoubletSource.jl")
import .DoubletSource: doublet_matrix, source_matrix, boundary_vector, wake_panels, source_strengths, surface_velocities, lift_coefficient, quadrilateral_source_velocity, quadrilateral_source_velocity_farfield, quadrilateral_source_potential, quadrilateral_doublet_velocity, quadrilateral_doublet_velocity_farfield, quadrilateral_doublet_potential

export doublet_matrix, source_matrix, boundary_vector, wake_panels, source_strengths, surface_velocities, lift_coefficient, quadrilateral_source_velocity, quadrilateral_source_velocity_farfield, quadrilateral_source_potential, quadrilateral_doublet_velocity, quadrilateral_doublet_velocity_farfield, quadrilateral_doublet_potential

## Linear-strength source and vorticity panel method

include("Aerodynamics/LinearVortexSource/LinearVortexSource.jl")
import .LinearVortexSource: total_velocity, source_velocity, vortex_velocity, vortex_influence_matrix, source_influence_matrix, neumann_boundary_condition, kutta_condition, two_point_neumann_matrix, linear_source_neumann_matrix, linear_vortex_neumann_matrix, constant_source_matrix, constant_source_boundary_condition

export total_velocity, source_velocity, vortex_velocity, vortex_influence_matrix, source_influence_matrix, neumann_boundary_condition, kutta_condition, two_point_matrix, linear_source_matrix, linear_vortex_matrix, constant_source_matrix, constant_source_boundary_condition, constant_quadrilateral_source_velocity, constant_quadrilateral_source_velocity_farfield, constant_quadrilateral_doublet_velocity, constant_quadrilateral_doublet_velocity_farfield

## Vortex lattice

include("Aerodynamics/VortexLattice/VortexLattice.jl")
import .VortexLattice: Horseshoe, AbstractVortexLatticeSystem, VortexLatticeSystem, References, AbstractAxisSystem, Stability, Wind, Body, Geometry, streamlines, influence_coefficient, influence_matrix, boundary_condition, solve_system, transform, bound_leg, bound_leg_center, bound_leg_vector, r1, r2, control_point, Horseshoe, surface_velocity, surface_forces, surface_moments, nearfield_drag, geometry_to_wind_axes, geometry_to_stability_axes, stability_to_geometry_axes, wind_to_geometry_axes, rate_coefficient, nearfield, farfield, farfield_forces, surface_velocities, surface_forces, surface_dynamics, surface_coefficients, nearfield_coefficients, farfield_coefficients, VortexRing, Freestream, velocity, kinematic_viscosity, mach_number, stream_velocity, center_of_pressure, freestream_derivatives!, freestream_derivatives, print_coefficients, print_derivatives

export Horseshoe, VortexLatticeSystem, References, AbstractAxisSystem, Stability, Wind, Body, Geometry, streamlines, influence_coefficient, influence_matrix, boundary_condition, solve_system, transform, bound_leg, bound_leg_center, bound_leg_vector, r1, r2, control_point, Horseshoe, surface_velocity, surface_forces, surface_moments, nearfield_drag, geometry_to_wind_axes, geometry_to_stability_axes, stability_to_geometry_axes, wind_to_geometry_axes, rate_coefficient, nearfield, farfield, farfield_forces, surface_velocities, surface_forces, surface_dynamics, surface_coefficients, nearfield_coefficients, farfield_coefficients, VortexRing, Freestream, velocity, kinematic_viscosity, mach_number, stream_velocity, center_of_pressure, freestream_derivatives!, freestream_derivatives, print_coefficients, print_derivatives

## Profile drag estimation
include("Aerodynamics/profile_drag.jl")

export profile_drag_coefficient, wetted_area_drag, local_dissipation_drag, form_factor

## Cases

include("Aerodynamics/Cases/printing.jl")

export print_info

include("Aerodynamics/Cases/cases.jl")

export solve_case, spanwise_loading, triangle_connectivities, extrapolate_point_mesh

include("Aerodynamics/Cases/stability_cases.jl")

export longitudinal_stability_derivatives, longitudinal_stability_matrix, lateral_stability_derivatives, lateral_stability_matrix

include("Aerodynamics/Cases/foil_cases.jl")

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

export plot_panels, plot_streamlines, plot_planform, plot_surface

end