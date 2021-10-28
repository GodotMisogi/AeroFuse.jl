## Cases
#==========================================================================================#

"""
    solve_case(horseshoe_panels :: Matrix{Panel3D}, normals, U, α, β, Ω, rho_ref, r_ref, area_ref, chord_ref, span_ref)

Evaluate a vortex lattice case given an array of `Panel3D`s with associated normal vectors (not necessarily the same as the panels' normals), freestream speed `U`, angles of attack `α` and sideslip `β`,  reference density ``\\rho``, reference point ``r_\\text{ref}`` for moments, and reference values for area, chord, and span lengths.
"""
function evaluate_case(horseshoe_panels :: Array{<: Panel3D}, normals, U, α, β, Ω, rho_ref, r_ref, area_ref, chord_ref, span_ref)
    # Make horseshoes and collocation points
    horseshoes = Horseshoe.(horseshoe_panels, normals)

    # Solve system
    Γs = reshape(solve_system(horseshoes[:], U, Ω), size(horseshoe_panels))

    # Compute forces and moments
    surface_forces, surface_moments, trefftz_force = case_dynamics(Γs, horseshoes, U, α, β, Ω, rho_ref, r_ref)

    # Compute aerodynamic coefficients
    nearfield_coeffs, farfield_coeffs, CFs, CMs = evaluate_coefficients(surface_forces, surface_moments, trefftz_force, U, α, β, rho_ref, area_ref, chord_ref, span_ref)

    nearfield_coeffs, farfield_coeffs, CFs, CMs, horseshoes, Γs
end

function evaluate_case(components, U, α, β, Ω, rho_ref, r_ref, area_ref, chord_ref, span_ref, name) 
    # Get dictionary keys and values, i.e. names and horseshoes
    comp_names     = (collect ∘ keys)(components)
    horseshoes_arr = values(components)
    horseshoes     = @views reduce(vcat, vec.(horseshoes_arr)) # Flattening for VLM system solution

    # Solve system
    Γs = solve_system(horseshoes, U, Ω)

    # Reshaping
    panel_sizes = size.(horseshoes_arr)
    panel_inds  = [ 0; cumsum(prod.(panel_sizes)) ]
    Γs_arr      = reshape_array(Γs, panel_inds, panel_sizes)

    # Compute forces and moments
    results = case_dynamics.(Γs_arr, horseshoes_arr, Ref(Γs), Ref(horseshoes), Ref(U), α, β, Ref(Ω), rho_ref, Ref(r_ref))
    forces  = getindex.(results, 1)
    moments = getindex.(results, 2) 
    trefftz = getindex.(results, 3)

    # Components' non-dimensional forces and moments
    data = evaluate_coefficients.(forces, moments, trefftz, Ref(U), α, β, rho_ref, area_ref, chord_ref, span_ref)

    nf_comp_coeffs = getindex.(data, 1)
    ff_comp_coeffs = getindex.(data, 2)
    CFs            = getindex.(data, 3)
    CMs            = getindex.(data, 4)

    # Aircraft's non-dimensional forces and moments
    nf_coeffs = reduce((x, y) -> x .+ y, nf_comp_coeffs) # Sum nearfield coefficients
    ff_coeffs = reduce((x, y) -> x .+ y, ff_comp_coeffs) # Sum farfield coefficients
    name_CFs  = reduce(vcat, vec.(CFs))                  # Collect surface force coefficients
    name_CMs  = reduce(vcat, vec.(CMs))                  # Collect surface moment coefficients

    # Dictionary assembly
    name_data = (nf_coeffs, ff_coeffs, name_CFs, name_CMs, Γs)
    comp_data = tuple.(nf_comp_coeffs, ff_comp_coeffs, CFs, CMs, Γs_arr)

    names   = [ name     ; # Aircraft name
                comp_names ] # Component names
    data    = [ name_data  ; # Aircraft data
                comp_data  ] # Component data

    OrderedDict(names .=> data)
end