module PotentialFlow

using LinearAlgebra
using StaticArrays
using Rotations, CoordinateTransformations
using SplitApplyCombine
using TimerOutputs
using LabelledArrays
using Accessors
using PrettyTables
using ForwardDiff: jacobian!
using DiffResults: JacobianResult, jacobian, value

## Package imports
#==========================================================================================#

# Non-dimensionalization
import ..NonDimensional: dynamic_pressure, aerodynamic_coefficients, force_coefficient, moment_coefficient, rate_coefficient

# Some tools
import ..Laplace: AbstractFreestream, Freestream

import ..AeroFuse: ComponentArray, velocity, solve_linear, solve_nonlinear, solve_nonlinear!, surface_velocities, surface_coefficients

## Element types and interfaces
#==========================================================================================#

"""
    AbstractPotentialFlowElement

Supertype for populated elements in a coupled velocity-based potential-flow solve.
Elements provide linear-in-strength `velocity`, `control_point`, `normal_vector`,
and `transform` methods. Compressible analyses also require
`prandtl_glauert_scale_coordinates`. Non-standard boundary conditions override
`apply_bc_row!`; wake-producing elements opt in through `has_wake`.

Nearfield evaluation requires `bound_leg_center`, `bound_leg_vector`, and
`trailing_velocity` adapters and, when needed, a `block_force_override` method.
Each named component contains one concrete element type; different components
may use different types. Strength units depend on the element family.
"""
abstract type AbstractPotentialFlowElement end

# Legacy qualified imports remain valid.
const AbstractVortex = AbstractPotentialFlowElement

include("vortices.jl")

## Reference frames
#==========================================================================================#

## Axis transformations with respect to freestream
abstract type AbstractAxisSystem end

struct Geometry <: AbstractAxisSystem end
struct Body <: AbstractAxisSystem end
struct Stability <: AbstractAxisSystem end
struct Wind <: AbstractAxisSystem end

"""
    velocity(freestream :: Freestream, ::Geometry)

Compute the velocity of Freestream in the geometry axis system.
"""
velocity(fs::Freestream, ::Geometry) = velocity(fs)
# velocity(fs :: Freestream, ::Body) = flip_xz(velocity(fs, Geometry()))
# velocity(fs :: Freestream, ::Stability) = 
# velocity(fs :: Freestream, ::Wind) = 

include("reference_frames.jl")

geometry_to_wind_axes(xyz, fs::Freestream) = geometry_to_wind_axes(xyz, fs.alpha, fs.beta)
geometry_to_wind_axes(vor::AbstractPotentialFlowElement, fs::Freestream) =
    geometry_to_wind_axes(vor, fs.alpha, fs.beta)

function geometry_to_wind_axes(vortex::AbstractPotentialFlowElement, α, β)
    T = promote_type(eltype(α), eltype(β))
    return transform(vortex, LinearMap(RotZY{T}(β, α)))
end

function wind_to_geometry_axes(vor::AbstractPotentialFlowElement, α, β)
    T = promote_type(eltype(α), eltype(β))
    return transform(vor, LinearMap(RotYZ{T}(-α, -β)))
end

# Prandl-Glauert transformation
include("prandtl_glauert.jl")

## Matrix and residual setups
#==========================================================================================#

include("residuals.jl")

## Fuselage line singularity (slender-body coupling)
include("fuselage_line.jl")

## Body source panel (skin-panelled fuselage coupling)
include("body_panel.jl")

"""
    solve_linear(elements, U, Ω)

Solve for element strengths using the boundary-condition velocity `U` and
quasi-steady rotation vector `Ω`. Return strengths, the influence matrix, and
the boundary vector after applying element-specific boundary-row overrides.
"""
solve_linear(elements, U, Ω) = solve_linear(elements, U, map(el -> zero(control_point(el)), elements), Ω)

"""
    solve_linear(elements, U, Ups, Ω)

Variant of [`solve_linear`](@ref) that injects an extra per-collocation-point velocity
field ``U_{ps}`` into the boundary condition. Used to couple the prescribed fuselage source
(thickness) line into the system while the fuselage doublet strengths solve as unknowns in
the AIC. `Ups` must share the layout of `elements`.

After the generic system is assembled, per-element boundary-condition overrides are applied
(see [`apply_bc_rows!`]) so element types with a non-standard boundary condition (e.g. the
slender-body cylinder condition of a `FuselageLine`) rewrite their own rows in place.
"""
function solve_linear(elements, U, Ups, Ω)
    AIC = influence_matrix(elements)
    boco = boundary_condition(elements, U, Ups, Ω)
    apply_bc_rows!(AIC, boco, elements, U, Ups, Ω)
    strengths = AIC \ boco

    return strengths, AIC, boco
end

## Force evaluations
#==========================================================================================#

# Nearfield forces
include("nearfield.jl")

# Farfield forces
include("farfield.jl")

# System setups
#==========================================================================================#

# System
include("system.jl")

const VortexLatticeSystem = PotentialFlowSystem

# Propeller slipstream (blown lift)
include("slipstream.jl")

# Derivatives
include("stability.jl")

## Post-processing
#==========================================================================================#

# Pretty-printing
include("printing.jl")

# Streamlines
include("streamlines.jl")

end
