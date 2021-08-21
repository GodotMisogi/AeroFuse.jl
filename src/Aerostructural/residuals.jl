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
    dxs, Ts = translations_and_rotations(δs)

    # New VLM variables
    new_vlm_mesh  = transfer_displacements(dxs, Ts, vlm_mesh, fem_mesh)
    new_cam_mesh = transfer_displacements(dxs, Ts, cam_mesh, fem_mesh)

    # Update states
    @timeit "New Horseshoes" aero_system.surfaces[surf_index].horseshoes = Horseshoe.(make_panels(new_vlm_mesh), panel_normal.(make_panels(new_cam_mesh)))
    aero_state.alpha = α

    # Solve VLM system and update forces
    @timeit "Aerodynamics Setup" solve_aerodynamics!(Γ, aero_system, aero_state)

    # Compute loads
    vlm_acs    = bound_leg_center.(horseshoes(aero_system.surfaces[1]))
    @timeit "Surface Aerodynamic Forces"  vlm_forces = surface_forces(aero_system.surfaces[1])
    
    # Build force vector with constraint
    @timeit "FEM Loads" fem_loads = fem_load_vector(vlm_acs, vlm_forces, fem_mesh)

    # Compute lift
    L = total_force(aero_system, aero_state)[3]

    # Aerodynamic residuals
    horsies = horseshoes(aero_system)
    @timeit "Aerodynamic Residual" aerodynamic_residual!(R_A, horsies, horseshoe_point.(horsies), horseshoe_normal.(horsies), Γ, -aero_state.velocity, aero_state.omega)

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
    δs      = @views reshape(δ[7:end], 6, length(fem_mesh))
    dxs, Ts = translations_and_rotations(δs)

    # New VLM variables
    new_vlm_mesh = transfer_displacements(dxs, Ts, vlm_mesh, fem_mesh)
    new_cam_mesh = transfer_displacements(dxs, Ts, cam_mesh, fem_mesh)

    new_horsies  = Horseshoe.(make_panels(new_vlm_mesh), panel_normal.(make_panels(new_cam_mesh)))

    # Set up VLM system
    U = freestream_to_cartesian(-speed, α, β)

    # Compute loads
    Γs         = reshape(Γ, size(new_horsies))
    vlm_acs    = bound_leg_center.(new_horsies)
    @timeit "Surface Forces"  vlm_forces = nearfield_forces(Γs, new_horsies, Γs, new_horsies, U, Ω, ρ)
    
    # Compute lift for load factor residual
    L = sum(vlm_forces)[3]

    # Build force vector with constraint for structures
    @timeit "FEM Loads" fem_loads = fem_load_vector(vlm_acs, vlm_forces, fem_mesh)

    # Aerodynamic residuals
    @timeit "Aerodynamic Residual" aerodynamic_residual!(R_A, new_horsies[:], horseshoe_point.(new_horsies[:]), horseshoe_normal.(new_horsies[:]), Γ / speed, U / speed, Ω / speed)

    # Structural residuals
    @timeit "Structural Residual" linear_residual!(R_S, stiffness_matrix, δ, fem_loads)

    # Weight residual
    @timeit "Load Factor Residual" linear_residual!(R_W, weight, load_factor, L * cos(α))

    return R
end

# Residual setup for multiple aerodynamic surfaces
function solve_coupled_residual!(R, x, speed, β, ρ, Ω, vlm_mesh, cam_mesh, other_horsies, fem_mesh, stiffness_matrix, weight, load_factor)
    n = length(x.aerodynamics) # Get size of VLM system (How do I do this without ComponentArrays now?)

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
    dxs, Ts = translations_and_rotations(δs)

    # New VLM variables
    new_vlm_mesh = transfer_displacements(dxs, Ts, vlm_mesh, fem_mesh)
    new_cam_mesh = transfer_displacements(dxs, Ts, cam_mesh, fem_mesh)

    new_horsies  = Horseshoe.(make_panels(new_vlm_mesh), panel_normal.(make_panels(new_cam_mesh)))
    all_horsies  = [ new_horsies[:]; other_horsies ]

    # Compute velocity with new angle of attack
    U = freestream_to_cartesian(-speed, α, β)

    # Compute forces for load factor residual
    @timeit "All Forces" all_forces = nearfield_forces(Γ, all_horsies, Γ, all_horsies, U, Ω, ρ)
    L = sum(all_forces)[3]

    # Compute component forces for structural residual
    new_Γs     = @views reshape(Γ[1:length(new_horsies)], size(new_horsies))
    new_acs    = bound_leg_center.(new_horsies)
    @timeit "Surface Forces" new_forces = nearfield_forces(new_Γs, new_horsies, Γ, all_horsies, U, Ω, ρ)

    # Build force vector with constraint for structures
    @timeit "FEM Loads" fem_loads = fem_load_vector(new_acs, new_forces, fem_mesh)

    # Aerodynamic residuals
    @timeit "Aerodynamic Residual" aerodynamic_residual!(R_A, all_horsies, horseshoe_point.(all_horsies), horseshoe_normal.(all_horsies), Γ / speed, U / speed, Ω / speed)

    # Structural residuals
    @timeit "Structural Residual" linear_residual!(R_S, stiffness_matrix, δ, fem_loads)

    # Weight residual
    @timeit "Load Factor Residual" linear_residual!(R_W, weight, load_factor, L * cos(α))

    return R
end

# # Residual setup for multiple structural and aerodynamic surfaces
# function solve_coupled_residual!(R, x, speed, β, ρ, Ω, vlm_meshes :: Vector{Matrix{SVector,3}}, cam_meshes :: Vector{Matrix{SVector,3}}, other_horsies, other_normies, fem_mesh, stiffness_matrix, weight, load_factor)
#     n = length(x.aerodynamics) # Get size of VLM system (How do I do this without ComponentArrays now?)

#     # Unpack aerodynamic and structural variables
#     Γ = @view x[1:n]
#     δ = @view x[n+1:end-1]
#     α = x[end]

#     # Get residual vector views
#     R_A = @view R[1:n]
#     R_S = @view R[n+1:end-1]
#     R_W = @view R[end]
    
#     # Compute displacements
#     δs  = @views reshape(δ[7:end], 6, length(fem_mesh))
#     dxs, Ts = translations_and_rotations(δs)

#     # New VLM variables
#     new_vlm_mesh = transfer_displacements(dxs, Ts, vlm_mesh, fem_mesh)
#     new_horsies  = Horseshoe.(make_panels(new_vlm_mesh))
#     all_horsies  = [ new_horsies[:]; other_horsies ]

#     new_cam_mesh = transfer_displacements(dxs, Ts, cam_mesh, fem_mesh)
#     all_normies  = [ panel_normal.(make_panels.(new_cam_meshes))[:]; other_normies ]

#     # Compute velocity with new angle of attack
#     U = freestream_to_cartesian(-speed, α, β)

#     # Compute forces for load factor residual
#     @timeit "All Forces" all_forces = nearfield_forces(Γ, all_horsies, Γ, all_horsies, U, Ω, ρ)
#     L = sum(all_forces)[3]

#     # Compute component forces for structural residual
#     new_Γs     = @views reshape(Γ[1:length(new_horsies)], size(new_horsies))
#     new_acs    = bound_leg_center.(new_horsies)
#     @timeit "Surface Forces" new_forces = nearfield_forces(new_Γs, new_horsies, Γ, all_horsies, U, Ω, ρ)

#     # Build force vector with constraint for structures
#     @timeit "FEM Loads" fem_loads = fem_load_vector(new_acs, new_forces, fem_mesh)

#     # Aerodynamic residuals
#     @timeit "Aerodynamic Residual" aerodynamic_residual!(R_A, all_horsies, horseshoe_point.(all_horsies), all_normies, Γ / speed, U / speed, Ω / speed)

#     # Structural residuals
#     @timeit "Structural Residual" linear_residual!(R_S, stiffness_matrix, δ, fem_loads)

#     # Weight residual
#     @timeit "Load Factor Residual" linear_residual!(R_W, weight, load_factor, L * cos(α))

#     return R
# end