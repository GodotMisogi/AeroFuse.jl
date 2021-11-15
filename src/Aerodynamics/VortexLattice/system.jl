## Cases
#==========================================================================================#

"""
    solve_case(horseshoe_panels :: Matrix{Panel3D}, normals, U, α, β, Ω, rho_ref, r_ref, area_ref, chord_ref, span_ref)

Evaluate a vortex lattice case given an array of `Panel3D`s with associated normal vectors (not necessarily the same as the panels' normals), freestream speed `U`, angles of attack `α` and sideslip `β`,  reference density ``\\rho``, reference point ``r_\\text{ref}`` for moments, and reference values for area, chord, and span lengths.
"""
function evaluate_case(horseshoes :: Matrix{Horseshoe{T}}, U, α, β, Ω, rho_ref, r_ref, area_ref, chord_ref, span_ref) where T <: Real
    V = norm(U)

    # Solve system
    Γs = reshape(solve_system(horseshoes[:], U, Ω), size(horseshoes))

    # Compute forces and moments
    surface_forces  = nearfield_forces(Γs, horseshoes, U, Ω, rho_ref)
    surface_moments = nearfield_moments(horseshoes, surface_forces, r_ref)
    trefftz_force   = trefftz_forces(Γs, horseshoes, V, α, β, rho_ref)

    # Compute dynamic pressure
    q = dynamic_pressure(rho_ref, V)

    # Non-dimensional panel coefficients
    CFs = force_coefficient.(surface_forces, q, area_ref)
    CMs = moment_coefficient.(surface_moments, q, area_ref, span_ref, chord_ref)

    # Compute summed coefficients
    force, moment = sum(surface_forces), sum(surface_moments)

    # Transform near-field dynamics to wind axes
    trans_force  = body_to_wind_axes(force, α, β)
    trans_force  = [ nearfield_drag(force, U); trans_force[2:end] ]
    trans_moment = body_to_wind_axes(stability_flip(moment), α, β)

    # Compute coefficients
    nearfield_coeffs = aerodynamic_coefficients(trans_force, trans_moment, V, area_ref, span_ref, chord_ref, rho_ref)
    farfield_coeffs  = force_coefficient(trefftz_force, q, area_ref)

    nearfield_coeffs, farfield_coeffs, CFs, CMs, horseshoes, Γs
end

function evaluate_case(components, U, α, β, Ω, rho_ref, r_ref, area_ref, chord_ref, span_ref) 
    # Solve vortex lattice system.
    # Thanks to the magic of ComponentArrays, the components are automatically treated as one big vector.
    Γs = solve_system(components, U, Ω)

    # Compute nearfield forces and moments
    surface_forces  = nearfield_forces(Γs, components, U, Ω, rho_ref) 
    surface_moments = nearfield_moments(components, surface_forces, r_ref)

    # Transform to wind axes
    wind_forces  = body_to_wind_axes.(surface_forces, α, β)
    wind_moments = body_to_wind_axes.(stability_flip.(surface_moments), α, β)
    # drags        = nearfield_drag.(surface_forces, Ref(U))

    # Compute farfield forces in Trefftz plane
    trefftz  = map(comp -> trefftz_forces(Γs[comp], components[comp], norm(U), α, β, rho_ref), valkeys(components))
    
    # Surface force, moment, and farfield coefficients
    q        = dynamic_pressure(rho_ref, norm(U))
    CFs_comp = force_coefficient.(wind_forces, q, area_ref)
    CMs_comp = moment_coefficient.(wind_moments, q, area_ref, span_ref, chord_ref) 
    FF_comp  = force_coefficient.(trefftz, q, area_ref)

    # Collect data for each component
    data_comp  = map((ff, comp) -> (ff, CFs_comp[comp], CMs_comp[comp], components[comp], Γs[comp]), FF_comp, valkeys(components))

    # Set up named tuples (somewhat inelegant, but understandable)
    properties = (:farfield, :CFs, :CMs, :horseshoes, :circulations)
    tuple_comp = @views NamedTuple{properties}.(data_comp)
    
    NamedTuple{keys(components)}(tuple_comp)
end

nearfield(comp) = [ sum(comp.CFs); sum(comp.CMs) ]
farfield(comp) = comp.farfield

circulations(data) = @views reduce(vcat, map(comp -> vec(data[comp].circulations), keys(data)))
nearfield_coefficients(data) =  @views map(comp -> nearfield(data[comp]), keys(data))
farfield_coefficients(data) = @views map(comp -> farfield(data[comp]), keys(data))

function aircraft_data(data)
    # # Dictionary assembly
    # name_data = VLMResults(nf_coeffs, ff_coeffs, name_CFs, name_CMs, horseshoes, Γs)
end

# Components' non-dimensional forces and moments
# data = @views map((comp, tref) -> evaluate_coefficients(forces[comp], moments[comp], tref, U, α, β, rho_ref, area_ref, chord_ref, span_ref), keys(components), trefftz)
# nf_coeffs_comp = @views getindex.(data, 1) # Nearfield coefficients
# ff_coeffs_comp = @views getindex.(data, 2) # Farfield coefficients

# Aircraft's non-dimensional forces and moments
# nf_coeffs = reduce((x, y) -> x .+ y, nf_coeffs_comp) # Sum nearfield coefficients
# ff_coeffs = reduce((x, y) -> x .+ y, ff_coeffs_comp) # Sum farfield coefficients
# CFs       = reduce(vcat, vec.(CFs_comp))             # Collect surface force coefficients
# CMs       = reduce(vcat, vec.(CMs_comp))             # Collect surface moment coefficients