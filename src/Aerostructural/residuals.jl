# Residual equation for linear systems
linear_residual!(R, A, x, b) = R .= A * x - b

## Nonlinear block Gauss-Seidel
#==========================================================================================#

function gauss_seidel(x0, speed, β, ρ, Ω, vlm_mesh, cam_mesh, fem_mesh, stiffness_matrix, weight, load_factor; max_iters = 50, tol = 1e-9)
    x = copy(x0)

    for i = 1:max_iters
        Γ = x.aerodynamics
        δ = x.structures
        α = x.load_factor

        # Compute velocity with new angle of attack
        U = freestream_to_cartesian(-speed, α, β)

        # Compute displacements
        δs      = δ.displacement
        dxs, Ts = translations_and_rotations(δs)

        # New VLM variables
        new_horsies = new_horseshoes(dxs, Ts, vlm_mesh, cam_mesh, fem_mesh)

        all_forces = surface_forces(Γ, all_horsies, U, Ω, ρ)

        fem_loads = compute_loads(bound_leg_center.(new_horsies), vlm_forces, fem_mesh) 
        
        δ = stiffness_matrix \ fem_loads

        # Compute lift for load factor residual
        D, Y, L = body_to_wind_axes(sum(vlm_forces), α, β)

        # Weight residual
        α = acos(weight * load_factor - L)

    end

    return x
end

## Fully-simultaneous Newton method
#==========================================================================================#

# Residual setup for single aerostructural surface
function solve_coupled_residual!(R, x, speed, β, ρ, Ω, vlm_mesh, cam_mesh, fem_mesh, stiffness_matrix, weight, load_factor)
    # Unpack aerodynamic and structural variables
    Γ = x.aerodynamics
    δ = x.structures
    α = x.load_factor

    # Get residual vector views
    R_A = R.aerodynamics
    R_S = R.structures
    R_W = @view R[end]
    
    # Compute velocity with new angle of attack
    U = freestream_to_cartesian(-speed, α, β)

    # Compute displacements
    δs      = δ.displacement
    dxs, Ts = translations_and_rotations(δs)

    # New VLM variables
    new_horsies = new_horseshoes(dxs, Ts, vlm_mesh, cam_mesh, fem_mesh)

    # Compute aerodynamic residual and loads
    vlm_forces = evaluate_aerodynamics!(R_A, Γ, new_horsies, U, Ω, ρ, speed)

    # Compute structural residual and loads
    fem_loads = evaluate_structures!(R_S, δ, bound_leg_center.(new_horsies), vlm_forces, fem_mesh, stiffness_matrix)

    # Compute lift for load factor residual
    D, Y, L = body_to_wind_axes(sum(vlm_forces), α, β)

    # Weight residual
    linear_residual!(R_W, weight, load_factor, L * cos(α))

    return R
end

# Residual setup for single aerostructural surface and multiple aerodynamic surfaces
function solve_coupled_residual!(R, x, speed, β, ρ, Ω, vlm_mesh, cam_mesh, other_horsies, fem_mesh, stiffness_matrix, weight, load_factor)
    # Unpack aerodynamic and structural variables
    Γ = x.aerodynamics
    δ = x.structures
    α = x.load_factor

    # Get residual vector views
    R_A = R.aerodynamics
    R_S = R.structures
    R_W = @view R[end]

    # Compute velocity with new angle of attack
    U = freestream_to_cartesian(-speed, α, β)
    
    # Compute displacements
    δs      = @views reshape(δ[7:end], 6, length(fem_mesh))
    dxs, Ts = translations_and_rotations(δs)

    # New VLM variables
    new_horsies = new_horseshoes(dxs, Ts, vlm_mesh, cam_mesh, fem_mesh)
    all_horsies = [ new_horsies[:]; other_horsies ]

    @timeit "Surface Forces" all_forces = surface_forces(Γ, all_horsies, U, Ω, ρ)

    new_forces = @views reshape(all_forces[1:length(new_horsies)], size(new_horsies))

    # Compute lift
    L = sum(all_forces)[3]

    # Build force vector with constraint for structures
    fem_loads = fem_load_vector(bound_leg_center.(new_horsies), new_forces, fem_mesh)

    # Aerodynamic residuals
    @timeit "Aerodynamic Residual" aerodynamic_residual!(R_A, all_horsies, Γ / speed, U / speed, Ω / speed)

    # Structural residuals
    linear_residual!(R_S, stiffness_matrix, δ, fem_loads)

    # Weight residual
    linear_residual!(R_W, weight, load_factor, L * cos(α))

    return R
end

# Residual setup for multiple aerostructural surfaces and multiple aerodynamic surfaces
function solve_coupled_residual!(R, x, speed, β, ρ, Ω, syms :: Vector{Symbol}, vlm_meshes, cam_meshes, fem_meshes, other_horsies, stiffness_matrix, weight, load_factor)
    # Unpack aerodynamic and structural variables
    Γs = x.aerodynamics
    δs = x.structures
    α  = x.load_factor

    # Get residual vector views
    R_A = R.aerodynamics
    R_S = R.structures
    R_W = @view R[end]

    # Compute velocity with new angle of attack
    U = freestream_to_cartesian(-speed, α, β)
    
    # Compute displacements
    Δs    = map((key, n) -> reshape(δs[key][7:end], 6, n), valkeys(δs), length.(fem_meshes))
    dx_Ts = translations_and_rotations.(Δs)
    dxs   = getindex.(dx_Ts, 1)
    Ts    = getindex.(dx_Ts, 2)

    # New VLM variables
    new_horsies  = new_horseshoes.(dxs, Ts, vlm_meshes, cam_meshes, fem_meshes)
    all_horsies  = [ reduce(vcat, vec.(new_horsies)); other_horsies ]

    # Compute component forces for structural residual
    new_Γs       = getindex.(Ref(Γs), syms) 
    new_acs      = map(horsies -> bound_leg_center.(horsies), new_horsies)
    @timeit "New Forces" new_forces   = surface_forces.(new_Γs, new_horsies, Ref(Γs), Ref(all_horsies), Ref(U), Ref(Ω), Ref(ρ))
    
    # Compute other forces for load factor residual
    other_syms   = filter(sym -> !(sym ∈ syms), keys(Γs))
    other_Γs     = reduce(vcat, vec.(getindex.(Ref(Γs), other_syms)))
    @timeit "Other Forces" other_forces = surface_forces(other_Γs, other_horsies, Γs, all_horsies, U, Ω, ρ)[:]

    # Compute lift
    L = sum([ reduce(vcat, vec.(new_forces)); other_forces[:] ])[3]

    # Build force vector with constraint for structures
    @timeit "FEM Loads" fem_loads = reduce(vcat, fem_load_vector.(new_acs, new_forces, fem_meshes))

    # Aerodynamic residuals
    @timeit "Aerodynamic Residuals" aerodynamic_residual!(R_A, all_horsies, Γs / speed, U / speed, Ω / speed)

    # Structural residuals
    @timeit "Structural Residuals" linear_residual!(R_S, stiffness_matrix, δs, fem_loads)

    # Weight residual
    linear_residual!(R_W, weight, load_factor, L * cos(α))

    return R
end

# Residual setup for multiple aerostructural surfaces (NEED TO REDUCE REDUNDANCIES)
function solve_coupled_residual!(R, x, speed, β, ρ, Ω, syms :: Vector{Symbol}, vlm_meshes, cam_meshes, fem_meshes, stiffness_matrix, weight, load_factor)
    # Unpack aerodynamic and structural variables
    Γs = x.aerodynamics
    δs = x.structures
    α  = x.load_factor

    # Get residual vector views
    R_A = R.aerodynamics
    R_S = R.structures
    R_W = @view R[end]

    # Compute velocity with new angle of attack
    U = freestream_to_cartesian(-speed, α, β)
    
    # Compute displacements
    Δs    = map((key, n) -> reshape(δs[key][7:end], 6, n), valkeys(δs), length.(fem_meshes))
    dx_Ts = translations_and_rotations.(Δs)
    dxs   = getindex.(dx_Ts, 1)
    Ts    = getindex.(dx_Ts, 2)

    # New VLM variables
    new_horsies  = new_horseshoes.(dxs, Ts, vlm_meshes, cam_meshes, fem_meshes)
    all_horsies  = reduce(vcat, vec.(new_horsies))

    # Compute component forces for structural residual
    new_Γs       = getindex.(Ref(Γs), syms) 
    new_acs      = map(horsies -> bound_leg_center.(horsies), new_horsies)
    @timeit "New Forces" new_forces   = surface_forces.(new_Γs, new_horsies, Ref(Γs), Ref(all_horsies), Ref(U), Ref(Ω), Ref(ρ))

    # Compute lift
    L = sum(reduce(vcat, vec.(new_forces)))[3]

    # Build force vector with constraint for structures
    @timeit "FEM Loads" fem_loads = reduce(vcat, fem_load_vector.(new_acs, new_forces, fem_meshes))

    # Aerodynamic residuals
    @timeit "Aerodynamic Residuals" aerodynamic_residual!(R_A, all_horsies, Γs / speed, U / speed, Ω / speed)

    # Structural residuals
    @timeit "Structural Residuals" linear_residual!(R_S, stiffness_matrix, δs, fem_loads)

    # Weight residual
    linear_residual!(R_W, weight, load_factor, L * cos(α))

    return R
end