mutable struct VLMState{T <: Real}
	U           :: T
    alpha       :: T
    beta        :: T
    omega       :: SVector{3,T}
	r_ref 		:: SVector{3,T}
	rho_ref 	:: T
	area_ref 	:: T
	chord_ref 	:: T
	span_ref 	:: T
	name 		:: String
end

VLMState(U :: T, α, β, Ω = zeros(3); rho_ref = 1.225, r_ref = zeros(3), area_ref = 1, chord_ref = 1, span_ref = 1, name = "Aircraft") where T <: Real = VLMState{T}(U, α, β, Ω, r_ref, rho_ref, area_ref, chord_ref, span_ref, name)

mutable struct VLMSurface{T <: Real}
    horseshoes         :: Matrix{Horseshoe{T}}
    collocation_points :: Matrix{SVector{3,T}}
    normals            :: Matrix{SVector{3,T}}
    surface_forces     :: Matrix{SVector{3,T}}
    surface_moments    :: Matrix{SVector{3,T}}
    circulations       :: Matrix{T}
    wake_vectors       :: Vector{SVector{3,T}}
    wake_AIC           :: Matrix{T}
    farfield_forces    :: MVector{3,T}
end

function VLMSurface(panels :: Matrix{Panel3D{T}}, normals :: Matrix{SVector{3,T}}) where T <: Real
    horseshoes         = horseshoe_line.(panels)
    collocation_points = horseshoe_point.(panels)
    
    m = size(panels)

    # Initialize
    wake_AIC        = Matrix{T}(undef, m)
    Γs              = Matrix{T}(undef, m)
    surface_forces  = Matrix{SVector{3,T}}(undef, m)
    surface_moments = Matrix{SVector{3,T}}(undef, m)
    wake_vectors    = Vector{SVector{3,T}}(undef, m[1])
    farfield_forces = MVector{3,T}(0., 0., 0.)

    VLMSurface{T}(horseshoes, collocation_points, normals, surface_forces, surface_moments, Γs, wake_vectors, wake_AIC, farfield_forces)
end

horseshoes(surf :: VLMSurface) = surf.horseshoes
collocation_points(surf :: VLMSurface) = surf.collocation_points
normals(surf :: VLMSurface) = surf.normals

function compute_wake_properties!(surface :: VLMSurface, α, β)
    # Reference velocity for broadcasting
    U_ref = (Ref ∘ SVector)(1, 0, 0)

    # Transform to wind axes
    wake_lines  = @. body_to_wind_axes(bound_leg(surface.horseshoes[end,:][:]), α, β)
    centers     = center.(wake_lines)
    wake_points = points(wake_lines)

    # Project trailing edge horseshoes' bound legs into Trefftz plane along wind axes
    surface.wake_vectors = @. project_vector(vector(wake_lines), U_ref)

    # Compute normal vectors
    wake_normals = @. normal(surface.wake_vectors, U_ref)

    surface.wake_AIC = trefftz_influence_matrix(centers, wake_normals, wake_points)
end

mutable struct VLMSystem{T <: Real}
    horseshoes         :: Vector{Horseshoe{T}}
    collocation_points :: Vector{SVector{3,T}}
    normals            :: Vector{SVector{3,T}}
    AIC                :: Matrix{T}
    RHS                :: Vector{T}
    circulations       :: Vector{T}
end

horseshoes(system :: VLMSystem) = system.horseshoes
collocation_points(system :: VLMSystem) = system.collocation_points
normals(system :: VLMSystem) = system.normals
AIC(system :: VLMSystem) = system.AIC
RHS(system :: VLMSystem) = system.RHS
circulations(system :: VLMSystem) = system.circulations

# Initialization
function VLMSystem{T}(horseshoes :: Vector{Horseshoe{T}}, collocation_points :: Vector{SVector{3,T}}, normals :: Vector{SVector{3,T}}) where T <: Real
    m   = length(horseshoes)
    AIC = Matrix{T}(undef, m, m)
    RHS = Vector{T}(undef, m)
    Γs  = Vector{T}(undef, m)
    VLMSystem{T}(horseshoes, collocation_points, normals, AIC, RHS, Γs)
end

function VLMSystem(surface :: VLMSurface{T}) where T <: Real
    horsies = horseshoes.(surface)[:]
    collies = collocation_points.(surface)[:]
    normies = normals.(surface)[:]
    VLMSystem{T}(horsies, collies, normies)
end

function VLMSystem(surfaces :: Vector{VLMSurface{T}}) where T <: Real
    # Flattening for system
    horsies = vcat((vec ∘ horseshoes).(surfaces)...)
    collies = vcat((vec ∘ collocation_points).(surfaces)...)
    normies = vcat((vec ∘ normals).(surfaces)...)
    VLMSystem{T}(horsies, collies, normies)
end

compute_horseshoes!(system :: VLMSystem, horseshoe_panels) = 
    system.horseshoes = horseshoe_line.(horseshoe_panels)

compute_collocation_points!(system :: VLMSystem, horseshoe_panels) = 
    system.collocation_points = horseshoe_point.(horseshoe_panels)

compute_influence_matrix!(system :: VLMSystem, V) = 
    system.AIC = influence_matrix(system.horseshoes, system.collocation_points, system.normals, -normalize(V))

compute_boundary_condition!(system :: VLMSystem, V, Ω) = 
    system.RHS = boundary_condition(map(r -> V + Ω × r, system.collocation_points), system.normals)

generate_system!(system :: VLMSystem, V, Ω) =
    matrix_assembly!(system.AIC, system.RHS, system.horseshoes, system.collocation_points, system.normals, V, Ω)

solve_system!(system :: VLMSystem) = 
    system.circulations = system.AIC \ system.RHS

solve_residual!(system :: VLMSystem, R, x) =
    R .= system.AIC * x - system.RHS

## Dynamics evaluations
compute_surface_forces!(surf :: VLMSurface, system :: VLMSystem, U, Ω, ρ) = 
    surf.surface_forces = nearfield_forces(surf.circulations, surf.horseshoes, system.circulations, system.horseshoes, U, Ω, ρ)

compute_surface_moments!(surf :: VLMSurface, r_ref) =
    surf.surface_moments = moments(surf.horseshoes, surf.surface_forces, r_ref)
    
function compute_farfield_forces!(surface :: VLMSurface, U, α, β, ρ);
    # Set up wake AIC and evaluate doublet normal derivatives
    Δφs = vec(sum(surface.circulations, dims = 1))
    compute_wake_properties!(surface, α, β)
    ∂φ_∂n = surface.wake_AIC * Δφs

    Δs = @. norm(surface.wake_vectors)
    θs = @. dihedral(surface.wake_vectors)

    surface.farfield_forces = trefftz_compute(Δφs, Δs, ∂φ_∂n, θs, U, ρ) 
end

## Evaluate system
function solve_case!(aircraft :: Dict{String, Tuple{Matrix{Panel3D{T}}, Matrix{SVector{3,T}}}}, state :: VLMState) where T <: Real
    # Collect surface data and freestream
    V        = freestream_to_cartesian(-state.U, state.alpha, state.beta)
    vals     = values(aircraft)
    horsies  = getindex.(vals, 1)
    normies  = getindex.(vals, 2)
    sizes    = size.(horsies)
    inds 	 = [ 0; cumsum(prod.(sizes)) ]

    # Build surfaces and systems
    surfs  = @. VLMSurface(horsies, normies)
    system = VLMSystem(surfs)

    # Assemble and solve matrix system
    # compute_influence_matrix!(system, V)
    # compute_boundary_condition!(system, V, state.omega)
    generate_system!(system, V, state.omega) # Pre-allocated version
    solve_system!(system)

    # Allocate surface circulations
    Γs = reshape_array(system.circulations, inds, sizes)
    for (surf, Γ_vec) in zip(surfs, Γs)
        surf.circulations = Γ_vec
    end

    # Evaluate forces
    compute_surface_forces!.(surfs, Ref(system), Ref(V), Ref(state.omega), state.rho_ref)
    compute_surface_moments!.(surfs, Ref(state.r_ref))
    compute_farfield_forces!.(surfs, state.U, state.alpha, state.beta, state.rho_ref)

    nf_coeffs = aerodynamic_coefficients.(surfs, Ref(state))
    ff_coeffs = farfield_force_coefficients.(surfs, Ref(state))
    
    total_nf_coeffs = sum(nf_coeffs)
    total_ff_coeffs = sum(ff_coeffs)

    # Generate dictionary
    results = Dict((collect ∘ keys)(aircraft) .=> zip(nf_coeffs, ff_coeffs, surfs))

    system, results, total_nf_coeffs, total_ff_coeffs
end

## Pure methods
surface_force_coefficients(surf :: VLMSurface, U, ρ, S) = 
    force_coefficient.(surf.surface_forces,  dynamic_pressure(ρ, U), S)

surface_moment_coefficients(surf :: VLMSurface, U, ρ, S, b, c) = 
    moment_coefficient.(surf.surface_moments, dynamic_pressure(ρ, U), S, b, c)

function aerodynamic_coefficients(surf :: VLMSurface, U, α, β, Ω, ρ, S, b, c) 
    CF_body = (sum ∘ surface_force_coefficients)(surf, U, ρ, S)
    CM_body = (sum ∘ surface_moment_coefficients)(surf, U, ρ, S, b, c)
    CF_wind = body_to_wind_axes(CF_body, α, β) # Consider axes specification in state
    CM_wind = body_to_stability_axes(stability_flip(CM_body), α)
    CR      = rate_coefficient(Ω, U, b, c)

    [ CF_wind; CM_wind; CR ]
end

surface_force_coefficients(surf :: VLMSurface, state :: VLMState) = surface_force_coefficients(surf, state.U, state.rho_ref, state.area_ref)
surface_moment_coefficients(surf :: VLMSurface, state :: VLMState) = surface_moment_coefficients(surf, state.U, state.rho_ref, state.area_ref, state.span_ref, state.chord_ref)
aerodynamic_coefficients(surf :: VLMSurface, state :: VLMState) = aerodynamic_coefficients(surf, state.U, state.alpha, state.beta, state.omega, state.rho_ref, state.area_ref, state.span_ref, state.chord_ref)

farfield_force_coefficients(surf :: VLMSurface, V, ρ, S) = force_coefficient(surf.farfield_forces, dynamic_pressure(ρ, V), S)
farfield_force_coefficients(surf :: VLMSurface, state :: VLMState) = farfield_force_coefficients(surf, state.U, state.rho_ref, state.area_ref)