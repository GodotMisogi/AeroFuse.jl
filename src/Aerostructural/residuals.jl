# Residual equation for linear systems
linear_residual!(R, A, x, b) = R .= A * x - b

# Residual setup for stateful style (seems pretty retarded right now...)
function solve_coupled_residual!(R, x, aero_system :: VLMSystem, aero_state :: VLMState, surf_index, vlm_mesh, cam_mesh, fem_mesh, stiffness_matrix, weight, load_factor)
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
    @timeit "Update VLM Mesh" new_vlm_mesh  = transfer_displacements(dxs, Ts, vlm_mesh, fem_mesh)
    @timeit "Update Camber Mesh" new_cam_mesh = transfer_displacements(dxs, Ts, cam_mesh, fem_mesh)

    # Update states
    @timeit "New Horseshoes" aero_system.surfaces[surf_index].horseshoes = horseshoe_line.(make_panels(new_vlm_mesh))
    @timeit "New Normals" aero_system.surfaces[surf_index].normals    = panel_normal.(make_panels(new_cam_mesh))
    @timeit "New Alpha" aero_state.alpha = α

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
    @timeit "Aerodynamic Residual" linear_residual!(R_A, AIC(aero_system), Γ, RHS(aero_system))

    # Structural residuals
    @timeit "Structural Residual" linear_residual!(R_S, stiffness_matrix, δ, fem_loads)

    # Weight residual
    @timeit "Load Factor Residual" linear_residual!(R_W, weight, load_factor, L * cos(α))

    return R
end

## "Pure" versions - AD compatible (for now)
#==========================================================================================#

# Residual setup for single surface
function solve_coupled_residual!(R, x, speed, β, ρ, Ω, vlm_mesh, cam_mesh, fem_mesh, stiffness_matrix, weight, load_factor)
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
    @timeit "Update VLM Mesh" new_vlm_mesh = transfer_displacements(dxs, Ts, vlm_mesh, fem_mesh)
    @timeit "Update Camber Mesh" new_cam_mesh = transfer_displacements(dxs, Ts, cam_mesh, fem_mesh)

    @timeit "New Horseshoes" new_horsies  = horseshoe_line.(make_panels(new_vlm_mesh))
    @timeit "New Normals" new_normies   = panel_normal.(make_panels(new_cam_mesh))

    # Set up VLM system
    U       = freestream_to_cartesian(-speed, α, β)

    # Compute loads
    @timeit "Surface Circulations" Γs         = reshape(Γ, size(new_horsies))
    @timeit "Surface Aerodynamic Centers" vlm_acs    = bound_leg_center.(new_horsies)
    @timeit "Surface Forces"  vlm_forces = nearfield_forces(Γs, new_horsies, Γs, new_horsies, U, Ω, ρ)
    
    # Compute lift for load factor residual
    @timeit "Compute Total Lift" L = sum(vlm_forces)[3]

    # Build force vector with constraint for structures
    @timeit "FEM Loads" fem_loads = fem_load_vector(vlm_acs, vlm_forces, fem_mesh)

    # Aerodynamic residuals
    @timeit "Aerodynamic Residual" aerodynamic_residual!(R_A, new_horsies[:], horseshoe_point.(new_horsies[:]), new_normies[:], Γ / speed, U / speed, Ω / speed)

    # Structural residuals
    @timeit "Structural Residual" linear_residual!(R_S, stiffness_matrix, δ, fem_loads)

    # Weight residual
    @timeit "Load Factor Residual" linear_residual!(R_W, weight, load_factor, L * cos(α))

    return R
end

# Residual setup for multiple surfaces
function solve_coupled_residual!(R, x, speed, β, ρ, Ω, vlm_mesh, cam_mesh, other_horsies, other_normies, fem_mesh, stiffness_matrix, weight, load_factor)
    n = length(x.aerodynamics) # Get size of VLM system (How do I do this without ComponentArrays now?)

    # Unpack aerodynamic and structural variables
    Γ = @views x[1:n]
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
    @timeit "Update VLM Mesh" new_vlm_mesh = transfer_displacements(dxs, Ts, vlm_mesh, fem_mesh)
    @timeit "New Horseshoes" new_horsies  = horseshoe_line.(make_panels(new_vlm_mesh))
    @timeit "All Horseshoes" all_horsies  = [ new_horsies[:]; other_horsies ]

    @timeit "Update Camber Mesh" new_cam_mesh = transfer_displacements(dxs, Ts, cam_mesh, fem_mesh)
    @timeit "New Normals" new_normies   = panel_normal.(make_panels(new_cam_mesh))
    @timeit "All Normals" all_normies   = [ new_normies[:]; other_normies ]

    # Compute velocity with new angle of attack
    U = freestream_to_cartesian(-speed, α, β)

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
    @timeit "Aerodynamic Residual" aerodynamic_residual!(R_A, all_horsies, horseshoe_point.(all_horsies), all_normies, Γ / speed, U / speed, Ω / speed)

    # Structural residuals
    @timeit "Structural Residual" linear_residual!(R_S, stiffness_matrix, δ, fem_loads)

    # Weight residual
    @timeit "Load Factor Residual" linear_residual!(R_W, weight, load_factor, L * cos(α))

    return R

end

# HEAVILY ALLOCATING AERODYNAMIC RESIDUAL COMPUTATIONS
# 
# Single surface:
# @timeit "Influence Matrix" inf_mat = influence_matrix(new_horsies[:], horseshoe_point.(new_horsies[:]), new_normies[:], -normalize(U), false)
# @timeit "Boundary Condition" boco    = boundary_condition(quasi_steady_freestream(new_horsies[:], U, Ω), new_normies[:])
# @timeit "Aerodynamic Residual" linear_residual!(R_A, inf_mat, Γ, boco)
#
# Multiple surfaces:
# @timeit "Normalizing Velocity" U_hat       = normalize(U)
# @timeit "Influence Matrix" inf_mat = influence_matrix(all_horsies, horseshoe_point.(all_horsies), normies[:], -U_hat, false)
# @timeit "Boundary Condition" boco    = boundary_condition(quasi_steady_freestream(all_horsies, U_hat, Ω / speed), normies[:])
# @timeit "Aerodynamic Residual" linear_residual!(R_A, inf_mat, Γ / speed, boco)