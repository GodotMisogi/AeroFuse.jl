## Cases
#==========================================================================================#

"""
    solve_case(horseshoe_panels :: Matrix{Panel3D}, normals, U, α, β, Ω, rho_ref, r_ref, area_ref, chord_ref, span_ref)

Evaluate a vortex lattice case given an array of `Panel3D`s with associated normal vectors (not necessarily the same as the panels' normals), freestream speed `U`, angles of attack `α` and sideslip `β`,  reference density ``\\rho``, reference point ``r_\\text{ref}`` for moments, and reference values for area, chord, and span lengths.
"""
function evaluate_case(horseshoes :: Matrix{Horseshoe{T}}, U, α, β, Ω, rho_ref, r_ref, area_ref, chord_ref, span_ref) where T <: Real
    # Solve system
    Γs = reshape(solve_system(horseshoes[:], U, Ω), size(horseshoes))

    # Compute forces and moments
    surface_forces, surface_moments, trefftz_force = evaluate_dynamics(Γs, horseshoes, U, α, β, Ω, rho_ref, r_ref)

    # Compute aerodynamic coefficients
    nearfield_coeffs, farfield_coeffs, CFs, CMs = evaluate_coefficients(surface_forces, surface_moments, trefftz_force, U, α, β, rho_ref, area_ref, chord_ref, span_ref)

    nearfield_coeffs, farfield_coeffs, CFs, CMs, horseshoes, Γs
end

function evaluate_case(components, U, α, β, Ω, rho_ref, r_ref, area_ref, chord_ref, span_ref) 
    # Get component names
    comp_names = keys(components)

    # Solve vortex lattice system.
    # Thanks to the magic of ComponentArrays, the components are automatically treated as one big vector.
    Γs = solve_system(components, U, Ω)

    hs_comp = getproperty.(Ref(components), comp_names)
    Γs_comp = getproperty.(Ref(Γs), comp_names)

    # Compute nearfield forces and moments
    forces  = nearfield_forces(Γs, components, Γs, components, U, Ω, rho_ref) 
    moments = nearfield_moments(components, forces, r_ref)

    # Transform to wind axes
    wind_forces  = body_to_wind_axes.(forces, α, β)
    wind_moments = body_to_wind_axes.(stability_flip.(moments), α, β)
    drags        = nearfield_drag.(forces, Ref(U))

    # Compute farfield forces in Trefftz plane
    trefftz  = trefftz_forces.(Γs_comp, hs_comp, norm(U), α, β, rho_ref)
    
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
    
    NamedTuple{comp_names}(tuple_comp)
end

nearfield(comp) = [ sum(comp.CFs); sum(comp.CMs) ]
farfield(comp) = comp.farfield

circulations(data) = reduce(vcat, map(comp -> vec(data[comp].circulations), keys(data)))
nearfield_coefficients(data) = map(comp -> nearfield(data[comp]), keys(data))
farfield_coefficients(data) = map(comp -> farfield(data[comp]), keys(data))

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