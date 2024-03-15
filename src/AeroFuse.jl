module AeroFuse

## Libraries
#==========================================================================================#

using LinearAlgebra
using StaticArrays
using CoordinateTransformations
using Rotations
using PrettyTables
using RecipesBase
using MacroTools

using Statistics: mean

using SplitApplyCombine: combinedimsview, combinedims, splitdimsview, splitdims
export combinedimsview, combinedims, splitdimsview, splitdims

using ComponentArrays
export ComponentVector, ComponentArray

using LabelledArrays

using Accessors
export @set

## Methods to be extended in submodules
#==========================================================================================#
function velocity end

function solve_linear end

function solve_nonlinear end

function solve_linear! end

function solve_nonlinear! end

function solve_system end

function surface_velocities end

function surface_coefficients end

export solve_linear, solve_nonlinear, solve_linear!, solve_nonlinear!

## Math tools
#==========================================================================================#

include("Tools/MathTools.jl")
# import .MathTools: weighted_vector, vectarray, slope, splitat, adj3, columns, extend_yz, reflect_mapper, cosine_interp, inverse_rotation, rotation, affine_2D, Point2D, Point3D, x, y, z, reshape_array, midpair_map, partition, uniform_spacing, linear_spacing, cosine_spacing, sine_spacing

# export forward_sum, forward_difference, forward_division, weighted_vector, vectarray, slope, splitat, adj3, columns, extend_yz, reflect_mapper, cosine_interp, structtolist, inverse_rotation, rotation, affine_2D, Point2D, Point3D, reshape_array, midpair_map, partition, uniform_spacing, linear_spacing, cosine_interpolation, sine_spacing

## Non-dimensionalization
#==========================================================================================#

include("Tools/NonDimensional.jl")
import .NonDimensional: dynamic_pressure, reynolds_number, force_coefficient, moment_coefficient, rate_coefficient, pressure_coefficient, aerodynamic_coefficients, force, moment

export dynamic_pressure, reynolds_number, force_coefficient, moment_coefficient, rate_coefficient, pressure_coefficient, force, moment

## Laplace
#==========================================================================================#

include("Tools/Laplace.jl")

import .Laplace: PointSingularity2D, PointSingularity3D, ConstantStrengthLineSingularity3D, potential, stream, velocity

export PointSingularity2D, PointSingularity3D, ConstantStrengthLineSingularity3D, potential, stream, velocity

# Traits
import .Laplace: Source2D, Doublet2D, Vortex2D, Uniform2D, Source3D, Doublet3D, SourceLine3D, DoubletLine3D

export Source2D, Doublet2D, Vortex2D, Uniform2D, Source3D, Doublet3D, SourceLine3D, DoubletLine3D

# Others
import .Laplace: Freestream, cartesian_to_freestream, freestream_to_cartesian

export Freestream, cartesian_to_freestream, freestream_to_cartesian

## Panels
#==========================================================================================#

include("Geometry/PanelGeometry/PanelGeometry.jl")
import .PanelGeometry: AbstractPanel, AbstractPanel2D, Panel2D, WakePanel2D, AbstractPanel3D, Panel3D, panel_length, transform_panel, transform_panel_points, panel_angle, tangent_vector, normal_vector, panel_location, panel_area, panel_coordinates, transform, midpoint, panel_points, wake_panel, wake_panels, reverse_panel, panel_velocity, panel_scalar, trailing_edge_panel, get_surface_values, panel_vector, distance, average_chord, average_width, wetted_area, make_panels, local_coordinate_system, get_transformation, trailing_edge_info, panel_coordinates, collocation_point

export AbstractPanel, AbstractPanel2D, Panel2D, WakePanel2D, AbstractPanel3D, Panel3D, transform, normal_vector, midpoint, panel_location, tangent_vector, panel_points, distance, wake_panel, wake_panels, panel_area, reverse_panel, panel_length, transform_panel, panel_angle, panel_vector, panel_velocity, panel_scalar, trailing_edge_panel, get_surface_values, average_chord, average_width, wetted_area, make_panels, local_coordinate_system, get_transformation, trailing_edge_info, panel_coordinates, collocation_point

## Aircraft geometry
#==========================================================================================#

include("Geometry/AircraftGeometry/AircraftGeometry.jl")

# Abstract types
import .AircraftGeometry: AbstractAircraft, AbstractWing, AbstractFoil, AbstractFuselage

export AbstractAircraft, AbstractWing, AbstractFoil, AbstractFuselage

# Foil
import .AircraftGeometry: Foil, arc_length, kulfan_CST, naca4, camber_CST, make_panels, read_foil, leading_edge_index, upper_surface, lower_surface, split_surface, coordinates_to_camber_thickness, camber_thickness_to_coordinates, camber_thickness, camber_thickness_to_coordinates, cosine_interpolation, camber_thickness_to_CST, coordinates_to_CST, maximum_thickness_to_chord, translate, interpolate, rotate, affine, scale, reflect, camber, camber_line, thickness_line, control_surface

export Foil, arc_length, kulfan_CST, naca4, camber_CST, make_panels, read_foil, leading_edge_index, upper_surface, lower_surface, split_surface, coordinates_to_camber_thickness, camber_thickness_to_coordinates, camber_thickness, camber_thickness_to_coordinates, cosine_interpolation, camber_thickness_to_CST, coordinates_to_CST, maximum_thickness_to_chord, translate, interpolate, rotate, affine, scale, reflect, camber, camber_line, thickness_line, control_surface

# Fuselage
import .AircraftGeometry: Fuselage, projected_area, length, cosine_interpolation, volume, HyperEllipseFuselage, curve

export Fuselage, projected_area, length, cosine_interpolation, volume, HyperEllipseFuselage, curve

# Wing
import .AircraftGeometry: Wing, WingSection, affine_transformation, mean_aerodynamic_chord, span, aspect_ratio, projected_area, taper_ratio, leading_edge, trailing_edge, chop_leading_edge, chop_trailing_edge, chop_wing, chop_sections, chop_coordinates, chop_spanwise_sections, chop_chords, chop_spans, make_panels, mesh_chords, mesh_wing, mesh_cambers, mean_aerodynamic_center, number_of_spanwise_panels, symmetric_spacing, coordinates, chord_coordinates, camber_coordinates, surface_coordinates, foils, chords, twists, spans, dihedrals, sweeps, position, orientation, WingMesh, chord_panels, camber_panels, surface_panels, AbstractSpacing, Sine, Cosine, Uniform, properties, wetted_area_ratio

export Wing, WingSection, affine_transformation, mean_aerodynamic_chord, span, aspect_ratio, projected_area, taper_ratio, leading_edge, trailing_edge, chop_leading_edge, chop_trailing_edge, chop_wing, chop_sections, chop_coordinates, chop_spanwise_sections, chop_chords, chop_spans, make_panels, mesh_chords, mesh_wing, mesh_cambers, mean_aerodynamic_center, number_of_spanwise_panels, symmetric_spacing, coordinates, chord_coordinates, camber_coordinates, surface_coordinates, foils, chords, twists, spans, dihedrals, sweeps, position, orientation, WingMesh, chord_panels, camber_panels, surface_panels, AbstractSpacing, Sine, Cosine, Uniform, properties, wetted_area_ratio

# Surfaces
# import .AircraftGeometry: HorizontalTail, VerticalTail

# export HorizontalTail, VerticalTail

# Controls
# import .AircraftGeometry: Flap, Aileron

# export Flap, Aileron

## Aerodynamic analyses
#==========================================================================================#

## Doublet-source panel method

include("Aerodynamics/DoubletSource/DoubletSource.jl")
import .DoubletSource: doublet_matrix, source_matrix, boundary_vector, wake_panels, source_strengths, surface_velocities, lift_coefficient, quadrilateral_source_potential, quadrilateral_doublet_potential

export doublet_matrix, source_matrix, boundary_vector, wake_panels, source_strengths, surface_velocities, lift_coefficient, quadrilateral_source_potential, quadrilateral_doublet_potential

## Linear-strength source and vorticity panel method

include("Aerodynamics/LinearVortexSource/LinearVortexSource.jl")
import .LinearVortexSource: total_velocity, source_velocity, vortex_velocity, vortex_influence_matrix, source_influence_matrix, neumann_boundary_condition, kutta_condition, two_point_neumann_matrix, linear_source_neumann_matrix, linear_vortex_neumann_matrix, constant_source_matrix, constant_source_boundary_condition

export total_velocity, source_velocity, vortex_velocity, vortex_influence_matrix, source_influence_matrix, neumann_boundary_condition, kutta_condition, two_point_matrix, linear_source_matrix, linear_vortex_matrix, constant_source_matrix, constant_source_boundary_condition, constant_quadrilateral_source_velocity, constant_quadrilateral_source_velocity_farfield, constant_quadrilateral_doublet_velocity, constant_quadrilateral_doublet_velocity_farfield

## Vortex lattice

include("Aerodynamics/VortexLattice/VortexLattice.jl")

# Vortex types
import .VortexLattice: Horseshoe, VortexRing, velocity, bound_leg_center, bound_leg_vector, control_point

export Horseshoe, VortexRing, velocity, bound_leg_center, bound_leg_vector, control_point

# Reference values
import .VortexLattice: References, kinematic_viscosity, mach_number

export References, kinematic_viscosity, mach_number

# Reference frames and traits
import .VortexLattice: AbstractAxisSystem, Stability, Wind, Body, Geometry, geometry_to_wind_axes, geometry_to_stability_axes, stability_to_geometry_axes, wind_to_geometry_axes, wind_to_body_axes

export AbstractAxisSystem, Stability, Wind, Body, Geometry, geometry_to_wind_axes, geometry_to_stability_axes, stability_to_geometry_axes, wind_to_geometry_axes, wind_to_body_axes

# System methods
import .VortexLattice: AbstractPotentialFlowSystem, VortexLatticeSystem, surface_velocity, surface_forces, surface_moments, nearfield_drag, rate_coefficient, nearfield, farfield, farfield_forces, surface_velocities, surface_forces, surface_dynamics, surface_coefficients, nearfield_coefficients, farfield_coefficients, center_of_pressure

export AbstractPotentialFlowSystem, VortexLatticeSystem, surface_velocity, surface_forces, surface_moments, nearfield_drag, rate_coefficient, nearfield, farfield, farfield_forces, surface_velocities, surface_forces, surface_dynamics, surface_coefficients, nearfield_coefficients, farfield_coefficients, center_of_pressure

# Derivatives
import .VortexLattice: freestream_derivatives

export freestream_derivatives

# Post-prrocessing
import .VortexLattice: print_coefficients, print_derivatives, streamlines

export print_coefficients, print_derivatives, streamlines

## Panel-VLM interface
include("Aerodynamics/vlm_interface.jl")

export make_horseshoes, make_vortex_rings

## Profile drag estimation
include("Aerodynamics/parasitic_drag.jl")

export form_factor, parasitic_drag_coefficient

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
import .Beams: Material, Tube, Beam, radii, area, moment_of_inertia, polar_moment_of_inertia, J_coeffs, Iyy_coeffs, Izz_coeffs, tube_stiffness_matrix, bending_stiffness_matrix, axial_stiffness_matrix, solve_cantilever_beam, elastic_modulus, shear_modulus, yield_stress, density, principal_stress, torsional_stress, von_mises_stress, beam_weight, structural_loads!, structural_loads

export Material, Tube, Beam, radii, area, moment_of_inertia, polar_moment_of_inertia, J_coeffs, Iyy_coeffs, Izz_coeffs, tube_stiffness_matrix, bending_stiffness_matrix, axial_stiffness_matrix, solve_cantilever_beam, elastic_modulus, shear_modulus, yield_stress, density, principal_stress, torsional_stress, von_mises_stress, beam_weight, structural_loads!, structural_loads

## Propulsion analyses
#==========================================================================================#

include("Propulsion/propulsion.jl")

# Actuator disc
import .Propulsion: ActuatorDisc, actuator_disc_induced_velocity

export ActuatorDisc, actuator_disc_induced_velocity

# Blade-element momentum theory
import .Propulsion: induced_velocity, induced_speed, inflow_angle, blade_solidity, slipstream_contraction

export induced_velocity, induced_speed, inflow_angle, blade_solidity, slipstream_contraction

## Post-processing
#==========================================================================================#

include("Tools/plot_tools.jl")

export plot_panel, plot_panels, plot_streamlines, plot_planform, plot_surface, plot_spanload

using PrecompileTools

@setup_workload begin
    # Putting some things in `@setup_workload` instead of `@compile_workload` can reduce the size of the
    # precompile file and potentially make loading faster.
    @compile_workload begin
        # all calls in this block will be precompiled, regardless of whether
        # they belong to your package or not (on Julia 1.8 and higher)
        include("precompile.jl")
    end
end


end