# Residual equation for linear systems
linear_residual!(R, A, x, b) = R .= A * x - b

## Fully-simultaneous Newton method
#==========================================================================================#

# Coupled aero-structural-load-factor residuals
function coupled_residuals!(R, all_horsies, Γs, U, Ω, speed, stiffness_matrix, δs, fem_loads, weight, load_factor, L)
    # Get residual vector views
    R_A = R.aerodynamics
    R_S = R.structures
    R_W = @view R[end]

    # Aerodynamic residuals
    @timeit "Aerodynamic Residuals" solve_nonlinear!(R_A, all_horsies, Γs / speed, U / speed, Ω / speed)

    # Structural residuals
    @timeit "Structural Residuals" linear_residual!(R_S, stiffness_matrix, δs, fem_loads)

    # Weight residual
    linear_residual!(R_W, weight, load_factor, L)

    nothing
end

# Residual setup for multiple aerostructural surfaces and multiple aerodynamic surfaces
function solve_coupled_residual!(R, x, speed, β, ρ, Ω, syms :: Vector{Symbol}, chord_meshes, camber_meshes, fem_meshes, other_horsies, stiffness_matrix, weight, load_factor)
    # Unpack aerodynamic and structural variables
    Γs = x.aerodynamics
    δs = x.structures
    α  = x.load_factor

    # Compute velocity with new angle of attack
    U = freestream_to_cartesian(-speed, α, β)
    
    # Compute displacements
    Δs    = map((key, fem_mesh) -> reshape(δs[key][7:end], 6, length(fem_mesh)), valkeys(δs), fem_meshes)
    # @timeit "Get Translations" dxs   = 
    # @timeit "Get Rotations" Ts    = 

    # New VLM variables
    @timeit "New Horseshoes" new_horsies = @. new_horseshoes(mesh_translation(Δs), mesh_rotation(Δs), chord_meshes, camber_meshes, fem_meshes)

    # Compute component forces for structural residual
    @timeit "Get Circulations" new_Γs       = getindex.(Ref(Γs), syms) 
    @timeit "Get Aerodynamic Centers" new_acs      = map(horsies -> bound_leg_center.(horsies), new_horsies)
    @timeit "All Horseshoes" all_horsies  = [ mapreduce(vec, vcat, new_horsies); vec(other_horsies) ]

    @timeit "New Forces" new_forces   = map((Γ_comp, hs_comp) -> surface_forces(hs_comp, Γ_comp, all_horsies, Γs, U, Ω, ρ), new_Γs, new_horsies)
    
    # Compute other forces for load factor residual
    @timeit "Get Symbols" other_syms   = filter(sym -> !(sym ∈ syms), keys(Γs))
    @timeit "Get Other Circulations" other_Γs     = mapreduce(sym -> vec(getindex(Γs, sym)), vcat, other_syms)
    @timeit "Other Forces" other_forces = surface_forces(vec(other_horsies), other_Γs, all_horsies, Γs, U, Ω, ρ)

    # Compute lift
    @timeit "All Forces" vlm_forces = [ mapreduce(vec, vcat, new_forces); vec(other_forces) ]
    @timeit "Transform Summed Forces" D, Y, L = geometry_to_wind_axes(sum(vlm_forces), α, β)

    # Build force vector with constraint for structures
    @timeit "FEM Loads" fem_loads = mapreduce(vec ∘ fem_load_vector, vcat, new_acs, new_forces, fem_meshes)

    # Compute residuals
    @timeit "Compute Residuals" coupled_residuals!(R, all_horsies, Γs, U, Ω, speed, stiffness_matrix, δs, fem_loads, weight, load_factor, L * cos(α))

    return R
end

# Residual setup for single aerostructural surface
function solve_coupled_residual!(R, x, speed, β, ρ, Ω, chord_mesh, camber_mesh, fem_mesh, stiffness_matrix, weight, load_factor)
    # Unpack aerodynamic and structural variables
    Γ = x.aerodynamics
    δ = x.structures
    α = x.load_factor
    
    # Compute velocity with new angle of attack
    U = freestream_to_cartesian(-speed, α, β)

    # Compute displacements
    Δs  = δ.displacement
    dxs = mesh_translation(Δs)
    Ts  = mesh_rotation(Δs)

    # New VLM variables
    @timeit "New Horseshoes" new_horsies = new_horseshoes(dxs, Ts, chord_mesh, camber_mesh, fem_mesh)

    # Compute aerodynamic forces
    @timeit "Surface Forces" vlm_forces = surface_forces(new_horsies, Γ, U, Ω, ρ)

    # Compute structural residual and loads
    @timeit "FEM Loads" fem_loads = fem_load_vector(bound_leg_center.(new_horsies), vlm_forces, fem_mesh)

    # Compute lift for load factor residual
    D, Y, L = geometry_to_wind_axes(sum(vlm_forces), α, β)

    # Compute residuals
    coupled_residuals!(R, new_horsies, Γ, U, Ω, speed, stiffness_matrix, δ, fem_loads, weight, load_factor, L * cos(α))

    return R
end

# # Residual setup for single aerostructural surface and multiple aerodynamic surfaces
# function solve_coupled_residual!(R, x, speed, β, ρ, Ω, chord_mesh, camber_mesh, other_horsies, fem_mesh, stiffness_matrix, weight, load_factor)
#     # Unpack aerodynamic and structural variables
#     Γ = x.aerodynamics
#     δ = x.structures
#     α = x.load_factor

#     # Get residual vector views
#     R_A = R.aerodynamics
#     R_S = R.structures
#     R_W = @view R[end]

#     # Compute velocity with new angle of attack
#     U = freestream_to_cartesian(-speed, α, β)
    
#     # Compute displacements
#     δs  = @views reshape(δ[7:end], 6, length(fem_mesh))
#     dxs = mesh_translation(Δs)
#     Ts  = mesh_rotation(Δs)

#     # New VLM variables
#     new_horsies = new_horseshoes(dxs, Ts, chord_mesh, camber_mesh, fem_mesh)
#     @timeit "Combine Horseshoes" all_horsies = [ vec(new_horsies); other_horsies ]

#     # Compute aerodynamic residual and loads
#     @timeit "Surface Forces" vlm_forces = surface_forces(Γ, all_horsies, U, Ω, ρ)

#     @timeit "FEM Forces" new_forces = @views reshape(vlm_forces[1:length(new_horsies)], size(new_horsies))

#     # Compute structural residual and loads
#     @timeit "FEM Loads" fem_loads = fem_load_vector(bound_leg_center.(new_horsies), new_forces, fem_mesh)

#     # Compute lift for load factor residual
#     D, Y, L = geometry_to_wind_axes(sum(vlm_forces), α, β)

#     # Compute residuals
#     coupled_residuals!(R, all_horsies, Γ, U, Ω, speed, stiffness_matrix, δ, fem_loads, weight, load_factor, L * cos(α))

#     return R
# end

# Residual setup for multiple aerostructural surfaces (NEED TO REDUCE REDUNDANCIES)
# function solve_coupled_residual!(R, x, speed, β, ρ, Ω, syms :: Vector{Symbol}, chord_meshes, camber_meshes, fem_meshes, stiffness_matrix, weight, load_factor)
#     # Unpack aerodynamic and structural variables
#     Γs = x.aerodynamics
#     δs = x.structures
#     α  = x.load_factor

#     # Compute velocity with new angle of attack
#     U = freestream_to_cartesian(-speed, α, β)
    
#     # Compute displacements
#     Δs    = map((key, n) -> reshape(δs[key][7:end], 6, n), valkeys(δs), length.(fem_meshes))
#     dxs   = mesh_translation.(Δs)
#     Ts    = mesh_rotation.(Δs)

#     # New VLM variables
#     new_horsies  = new_horseshoes.(dxs, Ts, chord_meshes, camber_meshes, fem_meshes)
#     all_horsies  = mapreduce(vec, vcat, new_horsies)

#     # Compute component forces for structural residual
#     new_Γs       = getindex.(Ref(Γs), syms) 
#     new_acs      = map(horsies -> bound_leg_center.(horsies), new_horsies)
#     @timeit "New Forces" new_forces = surface_forces.(new_Γs, new_horsies, Ref(Γs), Ref(all_horsies), Ref(U), Ref(Ω), Ref(ρ))

#     # Compute lift
#     D, Y, L = geometry_to_wind_axes(sum(mapreduce(vec, vcat, new_forces)), α, β)

#     # Build force vector with constraint for structures
#     @timeit "FEM Loads" fem_loads = mapreduce(fem_load_vector, vcat, new_acs, new_forces, fem_meshes)

#     # Compute residuals
#     coupled_residuals!(R, all_horsies, Γs, U, Ω, speed, stiffness_matrix, δs, fem_loads, weight, load_factor, L * cos(α))

#     return R
# end


## Nonlinear block Gauss-Seidel
#==========================================================================================#

function aerostruct_gauss_seidel(x0, speed, β, ρ, Ω, chord_mesh, camber_mesh, fem_mesh, stiffness_matrix, weight, load_factor; max_iters = 50, tol = 1e-9)
    x = deepcopy(x0)
    ε = 1e5
    for i = 1:max_iters
        xp = deepcopy(x)

        @show xp

        Γ = @views x.aerodynamics
        δ = @views x.structures
        α = @views x.load_factor

        # Compute velocity with new angle of attack
        U = freestream_to_cartesian(-speed, α, β)

        # Compute displacements
        δs  = δ.displacement
        dxs = mesh_translation(δs)
        Ts  = mesh_rotation(δs)

        # New geometric variables
        new_horsies = new_horseshoes(dxs, Ts, chord_mesh, camber_mesh, fem_mesh)

        # Solve circulations
        Γ = reshape(influence_matrix(vec(new_horsies), -U / speed) \ boundary_condition(quasi_steady_freestream(vec(new_horsies), U, Ω), horseshoe_normal.(vec(new_horsies))), size(new_horsies))

        x.aerodynamics = Γ

        @show x

        # Compute VLM forces
        vlm_forces = surface_forces(Γ, new_horsies, U, Ω, ρ)

        # Compute FEM loads
        fem_loads = fem_load_vector(bound_leg_center.(new_horsies), vlm_forces, fem_mesh) 
        
        # Solve displacements
        δ = reshape(stiffness_matrix \ vec(fem_loads), size(fem_loads))

        x.structures = δ

        @show x

        # Compute lift
        D, Y, L = geometry_to_wind_axes(sum(vlm_forces), α, β)

        # Weight residual
        α = acos(abs(weight * load_factor - L) / (weight * load_factor))

        x.load_factor = α

        @show L

        ε    = LinearAlgebra.norm(x - xp)
        @show (i, ε)

        if ε <= tol return x end # Needs NAN checks and everything like NLsolve

    end

    return x
end