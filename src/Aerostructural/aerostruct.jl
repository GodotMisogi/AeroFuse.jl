## Aerodynamic analysis
#==========================================================================================#

# lifting_line_forces(surf) = sum(surface_forces(surf), dims = 1)[:]
aerodynamic_forces(surfs) = reduce(vcat, (vec ∘ surface_forces).(surfs))

total_force(surfs, state) = body_to_wind_axes(sum(sum ∘ surface_forces, surfs), state.alpha, state.beta)

function solve_aerodynamics!(Γ, system :: VLMSystem, surfs :: Vector{<: VLMSurface}, state :: VLMState)
    # Update state velocity
    update_velocity!(state)

    # Assemble matrix system
    # compute_influence_matrix!(system, state.velocity)
    # compute_boundary_condition!(system, state.velocity, state.omega)
    @timeit "Generating AIC and RHS" generate_system!(system, state.velocity, state.omega) # Pre-allocated version for efficiency

    # Update circulations of system and surfaces
    @timeit "System Circulations"  system.circulations = Γ
    @timeit "Surface Circulations" update_circulations!(system)

    # Compute forces
    @timeit "Surface Forces"  map(surf -> compute_surface_forces!(surf, system, state.velocity, state.omega, state.rho_ref), surfs)
    @timeit "Surface Moments" map(surf -> compute_surface_moments!(surf, state.r_ref), surfs)
    @timeit "Farfield Forces" map(surf -> compute_farfield_forces!(surf, state.speed, state.alpha, state.beta, state.rho_ref), surfs)

    nothing
end

## Structural analysis
#==========================================================================================#

make_beam_mesh(vlm_mesh, fem_w) =  (1 - fem_w) * vlm_mesh[1,:] + fem_w * vlm_mesh[end,:]

## Axis transformation of stiffness matrices
function axis_transformation(fem_mesh, vlm_mesh)
    # Compute local beam axes
    ss = @. normalize(fem_mesh[2:end] - fem_mesh[1:end-1])
    ns = @. normalize((vlm_mesh[end,2:end] - vlm_mesh[1,1:end-1]) × (vlm_mesh[1,2:end] - vlm_mesh[end,1:end-1]))
    cs = @. ss × ns

    # WTF array of local coordinate systems - I can do better than this -_-
    wtf = zeros(3, 3, size(ss)...) # ((x,y,z), (c,s,n), chordwise, spanwise)
    for inds in CartesianIndices(wtf)
        l,k,i  = inds.I
        wtf[l,1,i] =  ns[i][l]
        wtf[l,2,i] = -cs[i][l]
        wtf[l,3,i] = -ss[i][l]
    end
    wtf
end

function transform_stiffy(perm_Ks, wtf)
    num = length(perm_Ks)
    K_tran = zeros(12, 12, num)
    perm_K = zeros(12, 12, num)
    for i in 1:num
        perm_K[:,:,i] .= perm_Ks[i]
        K_tran[:,:,i] .= kron(I(4), wtf[:,:,i])
    end

    ## Using Einstein summation convention for passive transformation of stiffness matrix: 
    #  D = Mᵗ(PK)M, where M = I₄ ⊗ T, and P is the permutation matrix acting on the original K
    @einsum D[j,k,i] := K_tran[l,j,i] * perm_K[l,m,i] * K_tran[m,k,i]
end

function permute_stiffy(Ks) 
    ## Simultaneous permutations of rows and columns:
    # 1. (Fx1, Fy1, Fz1, Mx1, My1, Mz1, Fx2, Fy2, Fz2, Mx2, My2, Mz2): [9,1,5,11,6,2,10,3,7,12,8,4]
    # 2. (Fx1, Fx2, Mx1, Mx2, Fy1, My1, Fy2, My2, Fz1, Mz1, Fz2, Mz2): [9,10,11,12,1,2,3,4,5,6,7,8]
    inds    = [9,1,5,11,6,2,10,3,7,12,8,4]
    perm_Ks = [ K[inds,:][:,inds] for K in Ks ] 
end

function build_big_stiffy(tubes, fem_mesh, vlm_mesh)
    # Building stiffness blocks
    Ks = [ tube_stiffness_matrix(tube) for tube in tubes ]

    perm_Ks = permute_stiffy(Ks)

    # Compute transformation matrix
    wtf = axis_transformation(fem_mesh, vlm_mesh)

    # Transform the stiffness matrix
    D = transform_stiffy(perm_Ks, wtf)    
end

solve_beam_residual!(R, K, δ, F) = R .= K * δ - F

## Load-displacement transfer mechanisms
#==========================================================================================#

# Sum adjacent values
adjacent_adder(x1, x2) = [ [ x1[1] ]; x1[2:end] .+ x2[1:end-1]; [ x2[end] ] ]

section_moments(vlm_acs, fem_pts, half_vlm_forces) = sum(x -> (x[1] .- fem_pts) .× x[2], zip(eachrow(vlm_acs), eachrow(half_vlm_forces)))

function compute_loads(vlm_acs, vlm_forces, fem_mesh)
    # Forces
    sec_forces   = sum(vlm_forces, dims = 1)[:] / 2
    beam_forces  = adjacent_adder(sec_forces / 2, sec_forces / 2)

    # Moments
    M_ins        = @views section_moments(vlm_acs, fem_mesh[1:end-1], vlm_forces / 2)
    M_outs       = @views section_moments(vlm_acs, fem_mesh[2:end],   vlm_forces / 2)
    beam_moments = adjacent_adder(M_ins, M_outs)

    # Concatenate forces and moments into loads array
    [ reduce(hcat, beam_forces); reduce(hcat, beam_moments) ]
end

# I guess JJ doesn't know differential forms.
rotation_matrix(Ωx, Ωy, Ωz) = @SMatrix [  0  -Ωz  Ωy ;
                                          Ωz  0  -Ωx ;
                                         -Ωy  Ωx  0  ]

rotation_matrix(θs) = rotation_matrix.(θs[1,:], θs[2,:], θs[3,:])

transfer_displacements(dxs, Ts, vlm_mesh, fem_mesh) = permutedims(reduce(hcat, map(xyz -> xyz + dxs + Ts .* (xyz - fem_mesh), eachrow(vlm_mesh))))

transfer_normals(Ts, normals) = permutedims(reduce(hcat, map(normie -> (Ts[1:end-1] + Ts[2:end]) / 2 .* normie, eachrow(normals))))

## Weights and fuel loads
#==========================================================================================#

load_factor_residual(L, W, n)     = L - n * W
load_factor_residual!(R, L, W, n) = R .= load_factor_residual(L, W, n)

## Coupled residual system
#==========================================================================================#

function solve_coupled_residual!(R, x, aero_system :: VLMSystem, aero_state :: VLMState, vlm_mesh, fem_mesh, stiffness_matrix, weight, load_factor)
    n = (prod ∘ size)(horseshoes(aero_system))   # Get panel size

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
    dxs = @views SVector.(δs[1,:], δs[2,:], δs[3,:])
    Ts  = @views rotation_matrix(δs[4:6,:])

    # New VLM variables
    new_vlm_mesh  = transfer_displacements(dxs, Ts, vlm_mesh, fem_mesh)
    new_panels    = make_panels(new_vlm_mesh)
    aero_system.surfaces[1].horseshoes = horseshoe_line.(new_panels)
    # aero_system.surfaces[1].normals    = transfer_normals(Ts, reshape(aero_system.surfaces[1].normals, size(new_panels)))
    aero_state.alpha = α

    # Solve VLM system and update forces
    @timeit "Aerodynamics Setup" solve_aerodynamics!(Γ, aero_system, aero_surfs, aero_state)

    # Compute loads
    @timeit "Get Aerodynamic Centers" vlm_acs    = bound_leg_center.(horseshoes(aero_system.surfaces[1]))
    @timeit "Get Aerodynamic Forces"  vlm_forces = surface_forces(aero_system.surfaces[1])
    
    # Build force vector with constraint
    @timeit "Computing Loads" fem_loads = [ zeros(6); compute_loads(vlm_acs, vlm_forces, fem_mesh)[:] ]

    # Compute lift
    @timeit "Computing Total Lift" L = total_force(aero_surfs, aero_state)[3]

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

evaluate_linear_residual!(R, A, x, b) = R .= A * x - b

# Single surface
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
    dxs = @views SVector.(δs[1,:], δs[2,:], δs[3,:])
    Ts  = @views rotation_matrix(δs[4:6,:])

    # New VLM variables
    new_vlm_mesh = transfer_displacements(dxs, Ts, vlm_mesh, fem_mesh)
    new_panels   = make_panels(new_vlm_mesh)
    new_horsies  = horseshoe_line.(new_panels)

    # Set up VLM system
    U       = freestream_to_cartesian(-speed, α, β)
    inf_mat = influence_matrix(new_horsies[:], horseshoe_point.(new_horsies[:]), normies[:], -normalize(U), false)
    boco    = boundary_condition(quasi_steady_freestream(new_horsies[:], U, Ω), normies[:])

    # Compute loads
    Γs         = reshape(Γ, size(new_horsies))
    vlm_acs    = bound_leg_center.(new_horsies)
    vlm_forces = nearfield_forces(Γs, new_horsies, Γs, new_horsies, U, Ω, ρ)
    
    # Build force vector with constraint for structures
    fem_loads = [ zeros(6); compute_loads(vlm_acs, vlm_forces, fem_mesh)[:] ]

    # Compute lift for load factor residual
    L = sum(vlm_forces)[3]

    # Aerodynamic residuals
    evaluate_linear_residual!(R_A, inf_mat, Γ, boco)

    # Structural residuals
    evaluate_linear_residual!(R_S, stiffness_matrix, δ, fem_loads)

    # Weight residual
    evaluate_linear_residual!(R_W, weight, load_factor, L * cos(α))

    return R
end

# Multiple surfaces
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
    dxs = @views SVector.(δs[1,:], δs[2,:], δs[3,:])
    Ts  = @views rotation_matrix(δs[4:6,:])

    # New VLM variables
    new_vlm_mesh = transfer_displacements(dxs, Ts, vlm_mesh, fem_mesh)
    new_panels   = make_panels(new_vlm_mesh)
    new_horsies  = horseshoe_line.(new_panels)
    all_horsies  = [ new_horsies[:]; horseshoe_line.(other_panels) ]

    # Set up VLM system
    U       = freestream_to_cartesian(-speed, α, β)
    inf_mat = influence_matrix(all_horsies, horseshoe_point.(all_horsies), normies[:], -normalize(U), false)
    boco    = boundary_condition(quasi_steady_freestream(all_horsies, U, Ω), normies[:])

    # Compute forces for load factor residual
    all_forces = nearfield_forces(Γ, all_horsies, Γ, all_horsies, U, Ω, ρ)
    L          = sum(all_forces)[3]

    # Compute component forces for structural residual
    new_Γs     = @views reshape(Γ[1:length(new_horsies)], size(new_horsies))
    new_acs    = bound_leg_center.(new_horsies)
    new_forces = nearfield_forces(new_Γs, new_horsies, Γ, all_horsies, U, Ω, ρ)

    # Build force vector with constraint for structures
    fem_loads = [ zeros(6); compute_loads(new_acs, new_forces, fem_mesh)[:] ]

    # Aerodynamic residuals
    evaluate_linear_residual!(R_A, inf_mat, Γ, boco)

    # Structural residuals
    evaluate_linear_residual!(R_S, stiffness_matrix, δ, fem_loads)

    # Weight residual
    evaluate_linear_residual!(R_W, weight, load_factor, L * cos(α))

    return R
end