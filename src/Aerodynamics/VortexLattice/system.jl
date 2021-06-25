mutable struct VLMState{T <: Real}
	U         :: T
    alpha     :: T
    beta      :: T
    V         :: SVector{3,T}
    omega     :: SVector{3,T}
	r_ref 	  :: SVector{3,T}
	rho_ref   :: T
	area_ref  :: T
	chord_ref :: T
	span_ref  :: T
	name 	  :: String
    VLMState(U, α, β, Ω :: AbstractVector{T} = zeros(3); r_ref = zeros(3), rho_ref = 1.225, area_ref = 1, chord_ref = 1, span_ref = 1, name = "Aircraft") where T <: Real = new{T}(U, α, β, freestream_to_cartesian(-U, α, β), SVector(Ω...), SVector(r_ref...), rho_ref, area_ref, chord_ref, span_ref, name)
end

rate_coefficient(state :: VLMState) = rate_coefficient(state.omega, state.U, state.span_ref, state.chord_ref)
update_velocity!(state :: VLMState) = 
    state.V = freestream_to_cartesian(-state.U, state.alpha, state.beta)

struct StabilityFrame{T <: Real}
    alpha :: T
end

struct WindFrame{T <: Real}
    alpha :: T
    beta  :: T
end

mutable struct VLMSurface{T <: Real}
    horseshoes      :: Matrix{Horseshoe{T}}
    normals         :: Matrix{SVector{3,T}}
    surface_forces  :: Matrix{SVector{3,T}}
    surface_moments :: Matrix{SVector{3,T}}
    circulations    :: Matrix{T}
    wake_vectors    :: Vector{SVector{3,T}}
    wake_AIC        :: Matrix{T}
    farfield_forces :: SVector{3,T}
end

function VLMSurface(panels :: Matrix{Panel3D{T}}, normals :: Matrix{SVector{3,T}}) where T <: Real
    # Initialize
    m               = size(panels)
    wake_AIC        = Matrix{T}(undef, m)
    Γs              = Matrix{T}(undef, m)
    surface_forces  = Matrix{SVector{3,T}}(undef, m)
    surface_moments = Matrix{SVector{3,T}}(undef, m)
    wake_vectors    = Vector{SVector{3,T}}(undef, m[1])
    farfield_forces = SVector{3,T}(0., 0., 0.)

    VLMSurface{T}(horseshoe_line.(panels), normals, surface_forces, surface_moments, Γs, wake_vectors, wake_AIC, farfield_forces)
end

horseshoes(surf :: VLMSurface)      = surf.horseshoes
normals(surf :: VLMSurface)         = surf.normals
surface_forces(surf :: VLMSurface)  = surf.surface_forces
surface_moments(surf :: VLMSurface) = surf.surface_moments
circulations(surf :: VLMSurface)    = surf.circulations

collocation_points(surf :: VLMSurface) = collocation_point.(horseshoes(surf))

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
    horseshoes   :: Vector{Horseshoe{T}}
    normals      :: Vector{SVector{3,T}}
    AIC          :: Matrix{T}
    RHS          :: Vector{T}
    circulations :: Vector{T}
end

horseshoes(system :: VLMSystem)   = system.horseshoes
normals(system :: VLMSystem)      = system.normals
AIC(system :: VLMSystem)          = system.AIC
RHS(system :: VLMSystem)          = system.RHS
circulations(system :: VLMSystem) = system.circulations

collocation_points(system :: VLMSystem) = collocation_point.(horseshoes(system))

# Initialization
function VLMSystem{T}(horseshoes :: Vector{Horseshoe{T}}, normals :: Vector{SVector{3,T}}) where T <: Real
    m   = length(horseshoes)
    AIC = Matrix{T}(undef, m, m)
    RHS = Vector{T}(undef, m)
    Γs  = Vector{T}(undef, m)
    VLMSystem{T}(horseshoes, normals, AIC, RHS, Γs)
end

function VLMSystem(surface :: VLMSurface{T}) where T <: Real
    horsies = horseshoes.(surface)[:]
    normies = normals.(surface)[:]
    VLMSystem{T}(horsies, normies)
end

function VLMSystem(surfaces :: Vector{VLMSurface{T}}) where T <: Real
    # Flattening for system
    horsies = reduce(vcat, (vec ∘ horseshoes).(surfaces))
    normies = reduce(vcat, (vec ∘ normals).(surfaces))
    VLMSystem{T}(horsies, normies)
end

function build_system(aircraft :: Dict{String, Tuple{Matrix{Panel3D{T}}, Matrix{SVector{3,T}}}}) where T <: Real
    # Get panels and normals
    vals   = values(aircraft)

    # Build surfaces and systems
    surfs  = @. VLMSurface(getindex.(vals, 1), getindex.(vals, 2))
    system = VLMSystem(surfs)

    system, surfs
end

compute_horseshoes!(system :: VLMSystem, horseshoe_panels) = 
    system.horseshoes = horseshoe_line.(horseshoe_panels)

compute_influence_matrix!(system, V) = 
    system.AIC = influence_matrix(horseshoes(system), collocation_points(system), normals(system), -normalize(V))

compute_boundary_condition!(system :: VLMSystem, V, Ω) = 
    system.RHS = boundary_condition(map(r -> V + Ω × r, collocation_points(system)), normals(system))

generate_system!(system :: VLMSystem, V, Ω) =
    matrix_assembly!(AIC(system), RHS(system), horseshoes(system), collocation_points(system), normals(system), V, Ω)

solve_system!(system :: VLMSystem) = 
    system.circulations = AIC(system) \ RHS(system)

function update_circulations!(Γ, surfs) 
    # Get sizes and indices for reshaping (UGLY AF)
    sizes = @. (size ∘ horseshoes)(surfs)
    inds  = [ 0; cumsum(prod.(sizes)) ]
    Γs    = reshape_array(Γ, inds, sizes)

    # Allocate surface circulations
    for (surf, Γ_vec) in zip(surfs, Γs)
        surf.circulations = Γ_vec
    end
end

## Dynamics evaluations
compute_surface_forces!(surf :: VLMSurface, system :: VLMSystem, U, Ω, ρ) = 
    surf.surface_forces = nearfield_forces(circulations(surf), horseshoes(surf), circulations(system), horseshoes(system), U, Ω, ρ)

compute_surface_moments!(surf :: VLMSurface, r_ref) =
    surf.surface_moments = moments(horseshoes(surf), surface_forces(surf), r_ref)
    
function compute_farfield_forces!(surface :: VLMSurface, U, α, β, ρ);
    # Set up wake AIC and evaluate doublet normal derivatives
    Δφs = vec(sum(circulations(surface), dims = 1))
    compute_wake_properties!(surface, α, β)
    ∂φ_∂n = surface.wake_AIC * Δφs

    Δs = @. norm(surface.wake_vectors)
    θs = @. dihedral(surface.wake_vectors)

    surface.farfield_forces = trefftz_compute(Δφs, Δs, ∂φ_∂n, θs, U, ρ) 
end

function solve_case!(system :: VLMSystem, surfs :: Vector{<: VLMSurface}, state :: VLMState)
    # Update state velocity
    update_velocity!(state)
    
    # Assemble and solve matrix system
    # compute_influence_matrix!(system, state.V)
    # compute_boundary_condition!(system, state.V, state.omega)
    generate_system!(system, state.V, state.omega) # Pre-allocated version for efficiency
    solve_system!(system)
    update_circulations!(circulations(system), surfs)

    # Evaluate forces
    compute_surface_forces!.(surfs, Ref(system), Ref(state.V), Ref(state.omega), state.rho_ref)
    compute_surface_moments!.(surfs, Ref(state.r_ref))
    compute_farfield_forces!.(surfs, state.U, state.alpha, state.beta, state.rho_ref)

    nothing
end

function solve_case!(aircraft :: Dict{String, Tuple{Matrix{Panel3D{T}}, Matrix{SVector{3,T}}}}, state :: VLMState) where T <: Real
    # Build surfaces and systems
    system, surfs = build_system(aircraft)

    # Solve case
    solve_case!(system, surfs, state)

    # Generate dictionaries
    results = Dict(keys(aircraft) .=> surfs)

    system, results
end

## Pure methods
surface_force_coefficients(surf :: VLMSurface, U, ρ, S) = 
    force_coefficient.(surf.surface_forces,  dynamic_pressure(ρ, U), S)

surface_moment_coefficients(surf :: VLMSurface, U, ρ, S, b, c) = 
    moment_coefficient.(surf.surface_moments, dynamic_pressure(ρ, U), S, b, c)

function nearfield_coefficients(surf :: VLMSurface, U, α, β, ρ, S, b, c) 
    CF_body = (sum ∘ surface_force_coefficients)(surf, U, ρ, S)
    CM_body = (sum ∘ surface_moment_coefficients)(surf, U, ρ, S, b, c)
    CF_wind = body_to_wind_axes(CF_body, α, β) # Consider axes specification in state
    CM_wind = body_to_wind_axes(stability_flip(CM_body), α, β)

    [ CF_wind; CM_wind ]
end

farfield_coefficients(surf :: VLMSurface, V, ρ, S) = force_coefficient(surf.farfield_forces, dynamic_pressure(ρ, V), S)

aerodynamic_coefficients(surf :: VLMSurface, state :: VLMState) = nearfield_coefficients(surf, state), farfield_coefficients(surf, state)

function aerodynamic_coefficients(surfs :: Vector{<: VLMSurface}, state :: VLMState) 
    coeffs    = aerodynamic_coefficients.(surfs, Ref(state))
    nf_coeffs = reduce(hcat, first.(coeffs))
    ff_coeffs = reduce(hcat, last.(coeffs))
    
    all_nf_coeffs = [ sum(nf_coeffs, dims = 2) nf_coeffs ]
    all_ff_coeffs = [ sum(ff_coeffs, dims = 2) ff_coeffs ]
    
    all_nf_coeffs, all_ff_coeffs
end

function aerodynamic_coefficients(surfs :: Dict{String, VLMSurface{T}}, state :: VLMState) where T <: Real 
    nf, ff    = aerodynamic_coefficients((collect ∘ values)(surfs), state)
    all_names = [ state.name; (collect ∘ keys)(surfs) ]
    Dict( all_names .=> zip(eachcol(nf), eachcol(ff)) )
end

# State versions
surface_force_coefficients(surf :: VLMSurface, state :: VLMState)  = surface_force_coefficients(surf, state.U, state.rho_ref, state.area_ref)

surface_moment_coefficients(surf :: VLMSurface, state :: VLMState) = surface_moment_coefficients(surf, state.U, state.rho_ref, state.area_ref, state.span_ref, state.chord_ref)

nearfield_coefficients(surf :: VLMSurface, state :: VLMState) = nearfield_coefficients(surf, state.U, state.alpha, state.beta, state.rho_ref, state.area_ref, state.span_ref, state.chord_ref)

farfield_coefficients(surf :: VLMSurface, state :: VLMState)  = farfield_coefficients(surf, state.U, state.rho_ref, state.area_ref)

## Residual setup
#=============================================#

evaluate_residual!(R, Γ, system :: VLMSystem) =
    R .= AIC(system) * Γ - RHS(system)
