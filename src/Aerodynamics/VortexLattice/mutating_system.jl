mutable struct VLMState{T <: Real}
    speed     :: T
    alpha     :: T
    beta      :: T
    velocity  :: SVector{3,T}
    omega     :: SVector{3,T}
    r_ref     :: SVector{3,T}
    rho_ref   :: T
    area_ref  :: T
    chord_ref :: T
    span_ref  :: T
    name      :: String
    VLMState(U, α, β, Ω :: AbstractVector{T} = zeros(3); r_ref = zeros(3), rho_ref = 1.225, area_ref = 1, chord_ref = 1, span_ref = 1, name = "Aircraft") where T <: Real = new{T}(U, α, β, freestream_to_cartesian(-U, α, β), SVector(Ω...), SVector(r_ref...), rho_ref, area_ref, chord_ref, span_ref, name)
end

rate_coefficient(state :: VLMState) = rate_coefficient(state.omega, state.speed, state.span_ref, state.chord_ref)
update_velocity!(state :: VLMState) = 
    state.velocity = freestream_to_cartesian(-state.speed, state.alpha, state.beta)

name(state :: VLMState) = state.name

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
    name            :: String
end

function VLMSurface(panels :: Matrix{Panel3D{T}}, normals :: Matrix{SVector{3,T}}, name = "Wing") where T <: Real
    # Initialize
    m               = size(panels)
    wake_AIC        = Matrix{T}(undef, m)
    Γs              = Matrix{T}(undef, m)
    surface_forces  = Matrix{SVector{3,T}}(undef, m)
    surface_moments = Matrix{SVector{3,T}}(undef, m)
    wake_vectors    = Vector{SVector{3,T}}(undef, m[1])
    farfield_forces = SVector{3,T}(0., 0., 0.)

    VLMSurface{T}(Horseshoe.(panels, normals), normals, surface_forces, surface_moments, Γs, wake_vectors, wake_AIC, farfield_forces, name)
end

horseshoes(surf :: VLMSurface)      = surf.horseshoes
normals(surf :: VLMSurface)         = surf.normals
surface_forces(surf :: VLMSurface)  = surf.surface_forces
surface_moments(surf :: VLMSurface) = surf.surface_moments
farfield_forces(surf :: VLMSurface) = surf.farfield_forces
circulations(surf :: VLMSurface)    = surf.circulations
name(surf :: VLMSurface)            = surf.name

collocation_points(surf :: VLMSurface) = horseshoe_point.(horseshoes(surf))

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
    surfaces     :: Vector{VLMSurface{T}}
    AIC          :: Matrix{T}
    RHS          :: Vector{T}
    circulations :: Vector{T}
end

surfaces(system :: VLMSystem)     = system.surfaces
AIC(system :: VLMSystem)          = system.AIC
RHS(system :: VLMSystem)          = system.RHS
circulations(system :: VLMSystem) = system.circulations
horseshoes(system :: VLMSystem)   = reduce(vcat, (vec ∘ horseshoes).(surfaces(system)))
normals(system :: VLMSystem)      = reduce(vcat, (vec ∘ normals).(surfaces(system)))

collocation_points(system :: VLMSystem) = horseshoe_point.(horseshoes(system))

# Initialization
function VLMSystem(surfaces :: AbstractVector{VLMSurface{T}}) where T <: Real
    m   = sum(length ∘ horseshoes, surfaces)
    AIC = Matrix{T}(undef, m, m)
    RHS = Vector{T}(undef, m)
    Γs  = Vector{T}(undef, m)
    VLMSystem{T}(surfaces, AIC, RHS, Γs)
end

function VLMSystem(surface :: VLMSurface{T}) where T <: Real
    horsies = horseshoes.(surface)[:]
    normies = normals.(surface)[:]
    VLMSystem{T}(horsies, normies)
end

# function VLMSystem(surfaces :: AbstractVector{VLMSurface{T}}) where T <: Real
#     # Flattening for system
#     horsies = reduce(vcat, (vec ∘ horseshoes).(surfaces))
#     normies = reduce(vcat, (vec ∘ normals).(surfaces))
#     VLMSystem{T}(horsies, normies)
# end

function build_system(aircraft :: Dict{String, Tuple{Matrix{Panel3D{T}}, Matrix{SVector{3,T}}}}) where T <: Real
    # Get panels and normals
    names  = keys(aircraft)
    vals   = values(aircraft)

    # Build surfaces and systems
    surfs  = @. VLMSurface(getindex.(vals, 1), getindex.(vals, 2), names)
    system = VLMSystem(surfs)

    system
end

compute_horseshoes!(system :: VLMSystem, horseshoe_panels) = 
    system.horseshoes = Horseshoe.(horseshoe_panels)

compute_influence_matrix!(system, V) = 
    system.AIC = influence_matrix(horseshoes(system), -normalize(V))

compute_boundary_condition!(system :: VLMSystem, V, Ω) = 
    system.RHS = boundary_condition(map(r -> V + Ω × r, collocation_points(system)), normals(system))

generate_system!(system :: VLMSystem, V, Ω) =
    matrix_assembly!(AIC(system), RHS(system), horseshoes(system), collocation_points(system), normals(system), V, Ω)

solve_system!(system :: VLMSystem) = 
    system.circulations = AIC(system) \ RHS(system)

function update_circulations!(system) 
    # Get sizes and indices for reshaping (UGLY AF)
    surfs = surfaces(system)
    sizes = @. (size ∘ horseshoes)(surfs)
    inds  = [ 0; cumsum(prod.(sizes)) ]
    Γs    = reshape_array(circulations(system), inds, sizes)

    # Allocate surface circulations
    for (surf, Γ_vec) in zip(surfs, Γs)
        surf.circulations = Γ_vec
    end
end


function evaluate_case!(system :: VLMSystem, state :: VLMState)
    # Update state velocity
    update_velocity!(state)

    # Get surfaces
    surfs = surfaces(system)
    
    # Assemble and solve matrix system
    compute_influence_matrix!(system, state.velocity)
    compute_boundary_condition!(system, state.velocity, state.omega)
    # generate_system!(system, state.velocity, state.omega) # Pre-allocated version for efficiency
    solve_system!(system)
    update_circulations!(system)

    # Evaluate forces
    compute_surface_forces!.(surfs, Ref(system), Ref(state.velocity), Ref(state.omega), state.rho_ref)
    compute_surface_moments!.(surfs, Ref(state.r_ref))
    compute_farfield_forces!.(surfs, state.speed, state.alpha, state.beta, state.rho_ref)

    system, state
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