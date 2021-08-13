# Residual equation for linear systems
evaluate_linear_residual!(R, A, x, b) = R .= A * x - b

# Residual setup for stateful style
function solve_coupled_residual!(R, x, aero_system :: VLMSystem, aero_state :: VLMState, vlm_mesh, fem_mesh, stiffness_matrix, weight, load_factor)
    n = (length ∘ horseshoes)(aero_system) # VLMSystem size

    # Unpack aerodynamic and structural variables
    Γ = @view x[1:n] 
    δ = @view x[n+1:end-1]
    α = x[end]

    # Get residual vector views
    R_A = @view R[1:n]
    R_S = @view R[n+1:end-1]
    R_W = @view R[end]
    
    # Compute displacements
    δs  = @views reshape(δ[7:end], 6, length(fem_mesh))
    @timeit "Build Displacements" dxs = @views SVector.(δs[1,:], δs[2,:], δs[3,:])
    @timeit "Build Rotations" Ts  = @views rotation_matrix(δs[4:6,:])

    # New VLM variables
    @timeit "Transfer Displacements" new_vlm_mesh  = transfer_displacements(dxs, Ts, vlm_mesh, fem_mesh)
    @timeit "New Panels" new_panels    = make_panels(new_vlm_mesh)
    @timeit "New Horseshoes" aero_system.surfaces[1].horseshoes = horseshoe_line.(new_panels)
    # aero_system.surfaces[1].normals    = transfer_normals(Ts, reshape(aero_system.surfaces[1].normals, size(new_panels)))
    aero_state.alpha = α

    # Solve VLM system and update forces
    @timeit "Aerodynamics Setup" solve_aerodynamics!(Γ, aero_system, aero_state)

    # Compute loads
    @timeit "Surface Aerodynamic Centers" vlm_acs    = bound_leg_center.(horseshoes(aero_system.surfaces[1]))
    @timeit "Surface Aerodynamic Forces"  vlm_forces = surface_forces(aero_system.surfaces[1])
    
    # Build force vector with constraint
    @timeit "FEM Loads" fem_loads = fem_load_vector(vlm_acs, vlm_forces, fem_mesh)

    # Compute lift
    @timeit "Compute Total Lift" L = total_force(aero_system, aero_state)[3]

    # Aerodynamic residuals
    @timeit "Aerodynamic Residual" evaluate_linear_residual!(R_A, AIC(aero_system), Γ, RHS(aero_system))

    # Structural residuals
    @timeit "Structural Residual" evaluate_linear_residual!(R_S, stiffness_matrix, δ, fem_loads)

    # Weight residual
    @timeit "Load Factor Residual" evaluate_linear_residual!(R_W, weight, load_factor, L * cos(α))

    return R
end

## "Pure" versions - AD compatible (for now)
#==========================================================================================#

# Residual setup for single surface
function solve_coupled_residual!(R, x, speed, β, ρ, Ω, vlm_mesh, normies, fem_mesh, stiffness_matrix, weight, load_factor)
    n = prod(size(vlm_mesh) .- 1) # Get size of VLM system

    # Unpack aerodynamic and structural variables
    Γ = @view x[1:n] 
    δ = @view x[n+1:end-1]
    α = x[end]

    # Get residual vector views
    R_A = @view R[1:n]
    R_S = @view R[n+1:end-1]
    R_W = @view R[end]
    
    # Compute displacements
    δs  = @views reshape(δ[7:end], 6, length(fem_mesh))
    @timeit "Build Displacements" dxs = @views SVector.(δs[1,:], δs[2,:], δs[3,:])
    @timeit "Build Rotations" Ts  = @views rotation_matrix(δs[4:6,:])

    # New VLM variables
    @timeit "Transfer Displacements" new_vlm_mesh = transfer_displacements(dxs, Ts, vlm_mesh, fem_mesh)
    @timeit "New Panels" new_panels   = make_panels(new_vlm_mesh)
    @timeit "New Horseshoes" new_horsies  = horseshoe_line.(new_panels)

    # Set up VLM system
    U       = freestream_to_cartesian(-speed, α, β)
    @timeit "Influence Matrix" inf_mat = influence_matrix(new_horsies[:], horseshoe_point.(new_horsies[:]), normies[:], -normalize(U), false)
    @timeit "Boundary Condition" boco    = boundary_condition(quasi_steady_freestream(new_horsies[:], U, Ω), normies[:])

    # Compute loads
    @timeit "Surface Circulations" Γs         = reshape(Γ, size(new_horsies))
    @timeit "Surface Aerodynamic Centers" vlm_acs    = bound_leg_center.(new_horsies)
    @timeit "Surface Forces"  vlm_forces = nearfield_forces(Γs, new_horsies, Γs, new_horsies, U, Ω, ρ)
    
    # Compute lift for load factor residual
    @timeit "Compute Total Lift" L = sum(vlm_forces)[3]

    # Build force vector with constraint for structures
    @timeit "FEM Loads" fem_loads = fem_load_vector(vlm_acs, vlm_forces, fem_mesh)

    # Aerodynamic residuals
    @timeit "Aerodynamic Residual" evaluate_linear_residual!(R_A, inf_mat, Γ, boco)

    # Structural residuals
    @timeit "Structural Residual" evaluate_linear_residual!(R_S, stiffness_matrix, δ, fem_loads)

    # Weight residual
    @timeit "Load Factor Residual" evaluate_linear_residual!(R_W, weight, load_factor, L * cos(α))

    return R
end

# Residual setup for multiple surfaces
function solve_coupled_residual!(R, x, speed, β, ρ, Ω, vlm_mesh, other_panels, normies, fem_mesh, stiffness_matrix, weight, load_factor)
    n = length(normies) # Get size of VLM system

    # Unpack aerodynamic and structural variables
    Γ = @view x[1:n]
    δ = @view x[n+1:end-1]
    α = x[end]

    # Get residual vector views
    R_A = @view R[1:n]
    R_S = @view R[n+1:end-1]
    R_W = @view R[end]
    
    # Compute displacements
    δs  = @views reshape(δ[7:end], 6, length(fem_mesh))
    @timeit "Build Displacements" dxs = @views SVector.(δs[1,:], δs[2,:], δs[3,:])
    @timeit "Build Rotations" Ts  = @views rotation_matrix(δs[4:6,:])

    # New VLM variables
    @timeit "Transfer Displacements" new_vlm_mesh = transfer_displacements(dxs, Ts, vlm_mesh, fem_mesh)
    @timeit "New Panels" new_panels   = make_panels(new_vlm_mesh)
    @timeit "New Horseshoes" new_horsies  = horseshoe_line.(new_panels)
    @timeit "All Horseshoes" all_horsies  = [ new_horsies[:]; horseshoe_line.(other_panels) ]

    # Set up VLM system
    U       = freestream_to_cartesian(-speed, α, β)
    @timeit "Influence Matrix" inf_mat = influence_matrix(all_horsies, horseshoe_point.(all_horsies), normies[:], -normalize(U), false)
    @timeit "Boundary Condition" boco    = boundary_condition(quasi_steady_freestream(all_horsies, U, Ω), normies[:])

    # Compute forces for load factor residual
    @timeit "All Forces" all_forces = nearfield_forces(Γ, all_horsies, Γ, all_horsies, U, Ω, ρ)
    @timeit "Compute Total Lift" L = sum(all_forces)[3]

    # Compute component forces for structural residual
    @timeit "Surface Circulations" new_Γs     = @views reshape(Γ[1:length(new_horsies)], size(new_horsies))
    @timeit "Surface Aerodynamic Centers" new_acs    = bound_leg_center.(new_horsies)
    @timeit "Surface Forces" new_forces = nearfield_forces(new_Γs, new_horsies, Γ, all_horsies, U, Ω, ρ)

    # Build force vector with constraint for structures
    @timeit "FEM Loads" fem_loads = fem_load_vector(new_acs, new_forces, fem_mesh)

    # Aerodynamic residuals
    @timeit "Aerodynamic Residual" evaluate_linear_residual!(R_A, inf_mat, Γ, boco)

    # Structural residuals
    @timeit "Structural Residual" evaluate_linear_residual!(R_S, stiffness_matrix, δ, fem_loads)

    # Weight residual
    @timeit "Load Factor Residual" evaluate_linear_residual!(R_W, weight, load_factor, L * cos(α))

    return R
end
