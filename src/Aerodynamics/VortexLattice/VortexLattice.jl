module VortexLattice

using LinearAlgebra
using StaticArrays
using Rotations
using OrderedCollections
using TimerOutputs
using ComponentArrays

## Package imports
#==========================================================================================#

# Math tools
import ..MathTools: weighted_vector, reshape_array

# Panel geometry
import ..PanelGeometry: Panel3D, panel_area, panel_coords, midpoint, panel_normal, transform, p1, p2, p3, p4, average_chord, average_width

# Non-dimensionalization
import ..NonDimensional: dynamic_pressure, aerodynamic_coefficients, force_coefficient, moment_coefficient, rate_coefficient

# Tools regarding solutions to Laplace's equation
import ..Laplace: cartesian_to_freestream, freestream_to_cartesian

## Horseshoe setup
#==========================================================================================#

include("horseshoes.jl")
include("finite_core.jl")

## Reference frames
#==========================================================================================#

include("reference_frames.jl")

## Influence matrix and solution of system
#==========================================================================================#

include("influences.jl")

"""
    make_horseshoes(horseshoe_panels)

Get the bound legs and collocation points of the horseshoes defined by horseshoe `Panel3D`s.
"""
make_horseshoes(horseshoe_panels) = @. Horseshoe(horseshoe_panels), horseshoe_point(horseshoe_panels)

quasi_steady_freestream(horseshoes, U, Ω) = map(hs -> U + Ω × horseshoe_point(hs), horseshoes)

"""
    solve_system(horseshoes, normals, U, Ω) 

Evaluate and return the vortex strengths ``\\Gamma``s given `Horseshoes`, their associated normal vectors (not necessarily the same as the panels' normals), the speed ``U`` and rotation vector ``\\Omega``.
"""
function solve_system(horseshoes, U, Ω, finite_core = false)
    AIC  = influence_matrix(horseshoes, -normalize(U), finite_core)
    boco = boundary_condition(quasi_steady_freestream(horseshoes, U, Ω), horseshoe_normal.(horseshoes))
    Γs   = AIC \ boco 
end

## Force evaluations
#==========================================================================================#

# Nearfield forces
include("nearfield.jl")

# Farfield forces
include("farfield.jl")


## Post-processing
#==========================================================================================#

"""
    evaluate_dynamics(Γ_comp, hs_comp, Γs, horseshoes, U, α, β, Ω, ρ, r_ref)

Compute the nearfield and farfield forces over a given component with associated vortex strengths `Γ_comp` and horseshoes `hs_comp`, using the vortex strengths ``\\Gamma``s and `horseshoes` of the entire configuration, the speed ``U``, angles of attack ``\\alpha`` and sideslip ``\\beta``, the rotation vector ``\\Omega``, the freestream density ``\\rho`` and the reference location for moments ``r_\\text{ref}``.
"""
function evaluate_dynamics(Γ_comp, hs_comp, Γs, horseshoes, U, α, β, Ω, ρ, r_ref)
    # Compute near-field dynamics
    surface_forces, surface_moments = nearfield_dynamics(Γ_comp, hs_comp, Γs, horseshoes, U, Ω, ρ, r_ref)

    # Compute farfield dynamics
    trefftz_force = trefftz_forces(Γ_comp, hs_comp, norm(U), α, β, ρ)

    surface_forces, surface_moments, trefftz_force
end

# In case of only one component
evaluate_dynamics(Γs, horseshoes, U, α, β, Ω, ρ, r_ref) = evaluate_dynamics(Γs, horseshoes, Γs, horseshoes, U, α, β, Ω, ρ, r_ref)
# evaluate_dynamics(Γs :: ComponentArray, horseshoes :: ComponentArray, U, α, β, Ω, ρ, r_ref) = evaluate_dynamics(Γs, horseshoes, Γs, horseshoes, U, α, β, Ω, ρ, r_ref)

function evaluate_coefficients(forces, moments, trefftz_force, U, α, β, ρ, S, c, b)
    V = norm(U)
    q = dynamic_pressure(ρ, V)

    # Computing summed coefficients
    force, moment = sum(forces), sum(moments)

    # Transform near-field dynamics to wind axes
    trans_force  = body_to_wind_axes(force, α, β)
    trans_force  = [ nearfield_drag(force, U); trans_force[2:end] ]
    trans_moment = body_to_wind_axes(stability_flip(moment), α, β)

    # Compute coefficients
    nearfield_coeffs = aerodynamic_coefficients(trans_force, trans_moment, V, S, b, c, ρ)
    farfield_coeffs  = force_coefficient(trefftz_force, q, S)

    # Non-dimensional panel coefficients
    CFs = force_coefficient.(forces, q, S)
    CMs = moment_coefficient.(moments, q, S, b, c)

    nearfield_coeffs, farfield_coeffs, CFs, CMs
end

# Streamlines
include("streamlines.jl")

# System setups
#==========================================================================================#

include("mutating_system.jl")
include("system.jl")

end