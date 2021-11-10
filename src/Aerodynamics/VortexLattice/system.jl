## Cases
#==========================================================================================#

struct VLMResults
    nearfield
    farfield
    surface_forces
    surface_moments
    horseshoes
    circulations
end

"""
    solve_case(horseshoe_panels :: Matrix{Panel3D}, normals, U, α, β, Ω, rho_ref, r_ref, area_ref, chord_ref, span_ref)

Evaluate a vortex lattice case given an array of `Panel3D`s with associated normal vectors (not necessarily the same as the panels' normals), freestream speed `U`, angles of attack `α` and sideslip `β`,  reference density ``\\rho``, reference point ``r_\\text{ref}`` for moments, and reference values for area, chord, and span lengths.
"""
function evaluate_case(horseshoes, U, α, β, Ω, rho_ref, r_ref, area_ref, chord_ref, span_ref)
    # Solve system
    Γs = reshape(solve_system(horseshoes[:], U, Ω), size(horseshoes))

    # Compute forces and moments
    surface_forces, surface_moments, trefftz_force = evaluate_dynamics(Γs, horseshoes, U, α, β, Ω, rho_ref, r_ref)

    # Compute aerodynamic coefficients
    nearfield_coeffs, farfield_coeffs, CFs, CMs = evaluate_coefficients(surface_forces, surface_moments, trefftz_force, U, α, β, rho_ref, area_ref, chord_ref, span_ref)

    nearfield_coeffs, farfield_coeffs, CFs, CMs, horseshoes, Γs
end

function evaluate_case(components, U, α, β, Ω, rho_ref, r_ref, area_ref, chord_ref, span_ref, name) 
    # Get dictionary keys and values, i.e. names and horseshoes
    # comp_names     = keys(components)
    # horseshoes_arr = values(components)
    # horseshoes     = @views reduce(vcat, vec.(horseshoes_arr)) # Flattening for VLM system solution

    # Solve vortex lattice system.
    # Thanks to the magic of ComponentArrays, the components are automatically treated as one big vector.
    Γs = solve_system(components, U, Ω)

    hs_comp = getproperty.(Ref(components), keys(components))
    Γs_comp = getproperty.(Ref(Γs), keys(components))

    # case_dynamics.(Γs_comp, hs_comp, Ref(Γs), Ref(horseshoes), Ref(U), α, β, Ref(Ω), rho_ref, Ref(r_ref))

    # Compute nearfield forces, moments and farfield forces in Trefftz plane
    forces = nearfield_forces(Γs, components, Γs, components, U, Ω, rho_ref) 
    # @views map(comp -> evaluate_dynamics(Γs[comp], components[comp], Γs, components, U, α, β, Ω, rho_ref, r_ref), keys(components))
    moments = nearfield_moments(components, forces, r_ref)
    trefftz = @views map(comp -> trefftz_forces(Γs[comp], components[comp], norm(U), α, β, rho_ref), keys(components))

    @show forces
    @show trefftz

    # Components' non-dimensional forces and moments
    # data = evaluate_coefficients.(forces, moments, trefftz, Ref(U), α, β, rho_ref, area_ref, chord_ref, span_ref)
    
    # nf_coeffs_comp = @views getindex.(data, 1) # Nearfield coefficients
    # ff_coeffs_comp = @views getindex.(data, 2) # Farfield coefficients
    # CFs_comp       = @views getindex.(data, 3) # Surface force coefficients
    # CMs_comp       = @views getindex.(data, 4) # Surface moment coefficients

    # # Aircraft's non-dimensional forces and moments
    # nf_coeffs = reduce((x, y) -> x .+ y, nf_coeffs_comp) # Sum nearfield coefficients
    # ff_coeffs = reduce((x, y) -> x .+ y, ff_coeffs_comp) # Sum farfield coefficients
    # name_CFs  = reduce(vcat, vec.(CFs))                  # Collect surface force coefficients
    # name_CMs  = reduce(vcat, vec.(CMs))                  # Collect surface moment coefficients

    # # Set up named tuples
    # properties = (:nearfield, :farfield, :CFs, :CMs, :horseshoes, :circulations)
    # full_data  = NamedTuple{properties}(nf_coeffs_comp, ff_coeffs_comp, CFs, CMs, components, Γs)
    # comp_data  = NamedTuple{properties}.(tuple.(nf_coeffs_comp, ff_coeffs_comp, CFs, CMs, hs_comp, Γs_comp))
    # results    = NamedTuple{(:full, keys(components)...)}(comp_data)


end

function aircraft_data(data)

    # Aircraft's non-dimensional forces and moments
    # nf_coeffs = reduce((x, y) -> x .+ y, nf_comp_coeffs) # Sum nearfield coefficients
    # ff_coeffs = reduce((x, y) -> x .+ y, ff_comp_coeffs) # Sum farfield coefficients
    # name_CFs  = reduce(vcat, vec.(CFs))                  # Collect surface force coefficients
    # name_CMs  = reduce(vcat, vec.(CMs))                  # Collect surface moment coefficients

    # # Dictionary assembly
    # name_data = VLMResults(nf_coeffs, ff_coeffs, name_CFs, name_CMs, horseshoes, Γs)
end