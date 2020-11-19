module VortexLattice

include("MathTools.jl")

using LinearAlgebra
using StaticArrays
using Rotations
using TimerOutputs
using .MathTools: accumap, structtolist

export Laplace, Uniform3D, Uniform, velocity, 
Panel, Panel3D, area, panel_coords, midpoint, panel_normal,
stability_axes, wind_axes,
solve_horseshoes, nearfield_dynamics, farfield_dynamics, nearfield_drag, 
streamlines, aerodynamic_coefficients, print_dynamics, horseshoe_lines

"""
Solutions to Laplace's equation.
"""
abstract type Laplace end

"""
A Uniform3D type expressing a uniform flow in spherical polar coordinates with angles α, β and magnitude mag.
"""
struct Uniform3D <: Laplace
    mag :: Real
    α :: Real 
    β :: Real
end

"""
Uniform3D constructor in degrees.
"""
Uniform(mag, α_deg, β_deg) = Uniform3D(-mag, deg2rad(α_deg), deg2rad(β_deg))

"""
Converts uniform flow coordinates to Cartesian coordinates.
"""
freestream_to_cartesian(r, θ, φ) = r .* SVector(cos(θ) * cos(φ), -sin(φ), sin(θ) * cos(φ))

"""
Computes the velocity of Uniform3D.
"""
velocity(uni :: Uniform3D) = freestream_to_cartesian(uni.mag, uni.α, uni.β)

# Axis transformations
#==========================================================================================#

"""
Converts coordinates into stability axes.
"""
stability_axes(coords, uni :: Uniform3D) =  [cos(uni.α) 0 sin(uni.α); 
                                                  0     1     0     ;
                                             sin(uni.α) 0 sin(uni.α)] * coords
                                            # RotY{Float64}(uni.α) * coords

"""
Converts coordinates into wind axes. This is supposedly not used and is just ceremonially implemented.
"""
wind_axes(coords, uni :: Uniform3D) = RotZY{Float64}(uni.β, uni.α) * coords


"""
Flips the y axis of a given vector.
"""
xz_flip(vector :: SVector{3, Float64}) = SVector(vector[1], -vector[2], vector[3])

"""
Flips the x and z axes of a given vector.
"""
stab_flip(vector :: SVector{3, Float64}) = SVector(-vector[1], vector[2], -vector[3])

"""
Transforms forces and moments into stability axes.
"""
stability_axes(force, moment, Ω, uniform :: Uniform3D) = stability_axes(force, uniform), stability_axes(stab_flip(moment), uniform), stability_axes(stab_flip(Ω), uniform)

"""
Transforms forces and moments into wind axes.
"""
wind_axes(force, moment, uniform :: Uniform3D) = wind_axes(force, uniform), wind_axes([ -moment[1], moment[2], -moment[3] ], uniform)


# Panel setup
#==========================================================================================#

"""
Placeholder. Panels should be an abstract type as they have some common methods, at least in 2D and 3D. 
"""
abstract type Panel end

"""
A composite type consisting of 4 coordinates. The following ASCII art depicts the order:

z → y
↓
x
        p1 —→— p4
        |       |
        ↓       ↓
        |       |
        p2 —→— p3
"""
struct Panel3D <: Panel
    p1 :: SVector{3,Real}
    p2 :: SVector{3,Real}
    p3 :: SVector{3,Real}
    p4 :: SVector{3,Real}
end

"""
Computes the area of Panel3D.
"""
area(panel :: Panel3D) = (abs ∘ norm)((panel.p2 - panel.p1) .* (panel.p3 - panel.p2))

"""
Computes the coordinates of a Panel3D for plotting purposes.
"""
panel_coords(panel :: Panel3D) = structtolist(panel)


"""
Computes the midpoint of Panel3D.
"""
midpoint(panel :: Panel3D) = (panel.p1 .+ panel.p2 .+ panel.p3 .+ panel.p4) / 4

"""
Computes the normal vector of Panel3D.
"""
panel_normal(panel :: Panel3D) = let p31 = panel.p3 .- panel.p1, 
                                     p42 = panel.p4 .- panel.p2, 
                                     p31_x_p42 = p31 × p42;
                                     p31_x_p42 / norm(p31_x_p42) end

"""
Computes the weighted value between two points.
"""
weighted(x1, x2, wx) = (1 - wx) * x1 + wx * x2

"""
Computes the weighted point between two points in a given direction.
"""
weighted_point((x1, y1, z1), (x2, y2, z2), wx, wy, wz) = SVector(weighted(x1, x2, wx), weighted(y1, y2, wy), weighted(z1, z2, wz))

"""
Computes the quarter point between two points in the x-z plane.
"""
quarter_point(p1, p2) = weighted_point(p1, p2, 1/4, 0, 1/4)

"""
Computes the 3-quarter point between two points in the x-z plane.
"""
three_quarter_point(p1, p2) = weighted_point(p1, p2, 3/4, 0, 3/4)

"""
Helper function to compute the bound leg of Panel3D for horseshoes/vortex rings, which is the quarter point on each side in the x-z plane.
"""
bound_leg(p1, p2, p3, p4) = [ quarter_point(p1, p2), quarter_point(p4, p3) ]

"""
Helper function to compute the collocation point of Panel3D for horseshoes/vortex rings, which is the 3-quarter point on each side in the x-z plane.
"""
collocation_point(p1, p2, p3, p4) = ( three_quarter_point(p1, p2) .+ three_quarter_point(p4, p3) ) ./ 2

"""
A composite type consisting of two vectors to define a line.
"""
struct Line
    r1 :: SVector{3,Real}
    r2 :: SVector{3,Real}
end

"""
Helper function to compute the velocity induced by a bound vortex leg.
"""
bound_leg_velocity(a, b, Γ) = Γ/4π * (1/norm(a) + 1/norm(b)) * a × b / (norm(a) * norm(b) + dot(a, b))

"""
Placeholder.
"""
bound_leg_velocity(r, line :: Line, Γ, ε = 1e-8) = let   a = r .- line.r1, b = r .- line.r2; on_line(a, b, ε) ? zeros(3) : bound_leg_velocity(a, b, Γ) end


"""
Helper function to compute the velocity induced by trailing vortex legs.
"""
trailing_legs_velocities(a, b, Γ, direction) = Γ/4π * ((a × direction) / (norm(a) - dot(a, direction)) / norm(a) - (b × direction) / (norm(b) - dot(b, direction)) / norm(b))

"""
Placeholder.
"""
trailing_legs_velocities(r, line :: Line, Γ, direction) = @timeit "Trailing Leg" let a = r .- line.r1, b = r .- line.r2; trailing_legs_velocities(a, b, Γ, direction) end

"""
Helper function to check if any point is on the line.
"""
on_line(a, b, ε) = any(<(ε), norm.([a, b, a × b]))


"""
Computes the velocity induced at a point `r` by a vortex Line with uniform strength Γ. Checks if `r` lies on the line and sets it to `(0, 0, 0)` if so, as the velocity is singular there.
"""
function horseshoe_velocity(r, line :: Line, Γ, ε = 1e-8; direction = [1, 0, 0])
    a = r .- line.r1
    b = r .- line.r2

    # Compute bound leg velocity
    @timeit "Bound Leg" bound_velocity = on_line(a, b, ε) ? zeros(3) : bound_leg_velocity(a, b, Γ)
    
    # Compute velocities of trailing legs
    @timeit "Trailing Leg" trailing_velocity = trailing_legs_velocities(a, b, Γ, direction)

    # Sums bound and trailing legs' velocities
    @timeit "Sum Legs" (bound_velocity .+ trailing_velocity)
end

"""
Placeholder. Vortex rings and horseshoes basically have the same methods, and are arrays of vortex lines.
"""
abstract type AbstractVortexArray end

"""
A horseshoe type consisting of a bound leg of type Line represening a vortex line.
"""
struct Horseshoe <: AbstractVortexArray
    bound_leg :: Line
end

"""
Computes the bound leg for a Panel3D.
"""
bound_leg(panel :: Panel3D) = bound_leg(panel.p1, panel.p2, panel.p3, panel.p4)

"""
Computes the collocation point of a Panel3D.
"""
collocation_point(panel :: Panel3D) = collocation_point(panel.p1, panel.p2, panel.p3, panel.p4)

"""
Returns a Horseshoe on a Panel3D.
"""
horseshoe_lines(panel :: Panel3D) = (Horseshoe ∘ Line)(bound_leg(panel)...)

"""
Computes the induced velocities at a point `r` of a Horseshoe with uniform strength Γ and trailing legs pointing in a given direction.
"""
velocity(r, horseshoe :: Horseshoe, Γ, û) = horseshoe_velocity(r, horseshoe.bound_leg, Γ, direction = û)

trailing_velocity(r, horseshoe :: Horseshoe, Γ, û) = trailing_legs_velocities(r, horseshoe.bound_leg, Γ, û)


"""
A vortex ring type consisting of vortex lines. TODO: Consider if better fields are "bound_vortices" and "trailing_vortices"
"""
struct VortexRing <: AbstractVortexArray
    left_leg :: Line
    bound_leg :: Line
    back_leg :: Line
    right_leg :: Line
end

transform(line :: Line, rotation, translation) = Line(rotation * line.r1 + translation, rotation * line.r2 + translation)

wind_axes(line :: Line, uniform :: Uniform3D) = Line(wind_axes(line.r1, uniform), wind_axes(line.r2, uniform)) 

"""
Helper function to compute the vortex ring given four points following Panel3D ordering.
"""
function vortex_ring(p1, p2, p3, p4)
    v1, v4 = quarter_point(p1, p2), quarter_point(p4, p3)
    v2, v3 = v1 .+ p2 .- p1, v4 .+ p3 .- p4
    [ v1, v2, v3, v4 ]
end

"""
Computes the vortex rings on a Panel3D.
"""
vortex_ring(panel :: Panel3D) = vortex_ring(panel.p1, panel.p2, panel.p3, panel.p4)

"""
Constructor for vortex rings on a Panel3D using Lines. The following convention is adopted:

``` 
    p1 —bound_leg→ p4
    |               |
left_line       right_line
    ↓               ↓
    p2 —back_line→ p3
```
"""
function VortexRing(panel :: Panel3D)
    v1, v2, v3, v4 = vortex_ring(panel)

    bound_leg = Line(v1, v4)
    left_leg = Line(v2, v1)
    right_leg = Line(v3, v2)
    back_leg = Line(v4, v3)

    VortexRing(left_leg, bound_leg, back_leg, right_leg)
end

"""
Computes the midpoint of the bound leg of a Horseshoe or Vortex Ring.
"""
bound_leg_center(vortex :: AbstractVortexArray) = (vortex.bound_leg.r1 .+ vortex.bound_leg.r2) / 2

"""
Computes the direction vector of the bound leg of a Horseshoe or Vortex Ring.
"""
bound_leg_vector(vortex :: AbstractVortexArray) = vortex.bound_leg.r2 .- vortex.bound_leg.r1

"""
Sums the velocities evaluated at a point `r` of vortex lines with uniform strength Γ.
"""
sum_vortices(r, vortex_lines :: Array{Line}, Γ) = sum(velocity(r, line, Γ) for line ∈ vortex_lines)

"""
Computes the induced velocities at a point `r` of a Vortex Ring with uniform strength Γ.
"""
velocity(r, vortex_ring :: VortexRing, Γ) = sum_vortices(r, structtolist(vortex_ring), Γ)

# Matrix setup
#==========================================================================================#

"""
Computes the velocity using the method of images for a symmetric case in the x-z plane.
"""
function mirror_velocity(collocation_point, horseshoe, Γ, û)
    mirror_point = xz_flip(collocation_point)
    mir_vel = (xz_flip ∘ velocity)(mirror_point, horseshoe, Γ, û)
end

"""
Computes the influence coefficient of the velocity of a vortex line at a collocation point projected to a normal vector.
"""
function influence_coefficient(collocation_point, horseshoe :: Horseshoe, panel_normal, û, symmetry = false) 
    if symmetry
        col_vel = velocity(collocation_point, horseshoe, 1., û)
        mir_vel = mirror_velocity(collocation_point, horseshoe, 1., û)
        return dot(col_vel .+ mir_vel, panel_normal)
    else 
        return dot(velocity(collocation_point, horseshoe, 1., û), panel_normal)
    end
end

"""
Assembles the Aerodynamic Influence Coefficient (AIC) matrix given vortex lines, collocation points and associated normal vectors.
"""
influence_matrix(colpoints, normals, vortex_lines :: Array{Horseshoe}, û, symmetry = false) = @timeit "Matrix Construction" [ @timeit "Influence Coefficient" influence_coefficient(colpoint_i, horsie_j, normal_i, û, symmetry) for (colpoint_i, normal_i) ∈ zip(colpoints, normals), horsie_j ∈ vortex_lines ]

"""
Computes the projection of a velocity vector with respect to normal vectors of panels. Corresponds to construction of the boundary condition for the RHS of the AIC system.
"""
boundary_condition(velocities, normals) = @timeit "Dotting" dot.(velocities, normals)

"""
Generates a wake Panel3D aligned with a normalised direction vector.
"""
function wake_panel(panel :: Panel3D, uniform :: Uniform3D, bound = 1.)
    wake = bound .* velocity(uniform) / uniform.mag
    Panel3D(panel.p2, panel.p2 .+ wake, panel.p3 .+ wake, panel.p3)
end

"""
Generates wake Panel3Ds given the Panel3Ds of a wing corresponding to the trailing edge and a Uniform3D. The wake length and number of wake panels must be provided. 
"""
make_wake(last_panels :: Array{Panel3D}, uniform :: Uniform3D, wake_length, wake_num) = (permutedims ∘ accumap)(panel -> wake_panel(panel, uniform, wake_length / wake_num), wake_num, last_panels)

"""
Solves the AIC matrix with the boundary condition given Panel3Ds and a freestream velocity.
"""
function solve_horseshoes(horseshoe_panels :: Array{Panel3D}, camber_panels :: Array{Panel3D}, û, Ω̄, symmetry = false) 

    @timeit "Horseshoes" horseshoes = horseshoe_lines.(horseshoe_panels)
    @timeit "Collocation Points" colpoints = collocation_point.(horseshoe_panels)
    @timeit "Normals" normals = panel_normal.(camber_panels)
    @timeit "Total Velocity" total_vel = [ û .+ (Ω̄ × rc_i) for rc_i ∈ colpoints ]

    @timeit "AIC" AIC = influence_matrix(colpoints[:], normals[:], horseshoes[:], û, symmetry) 
    @timeit "RHS" boco = boundary_condition(normals[:], total_vel[:])
    @timeit "Solve AIC" Γs = AIC \ boco

    @timeit "Reshape" output = reshape(Γs, size(horseshoes)...), horseshoes, AIC

    output
end

# Force evaluations
#==========================================================================================#

# Trefftz plane evaluations
trefftz_potential(r_i, r_j, Γ_j, û) = let r = r_i .- r_j; Γ_j/2π * û × r / dot(r, r) end

trefftz_matrix(horseshoes, normals, û) = [ dot(trefftz_potential(bound_leg_center(horseshoe_i), horseshoe_j.bound_leg.r1, 1., û), n̂_i) for (horseshoe_i, n̂_i) in zip(horseshoes, normals), horseshoe_j in horseshoes ]

"""
Computes the aerodynamic forces in the Trefftz plane.
"""
function trefftz_forces(Γs, horseshoes, freestream, ρ)
    # Transform horseshoes into Trefftz plane along freestream
    lines = bound_leg_vector.(horseshoes[end,:])
    vels = repeat([freestream], length(lines))
    projs = dot.(vels, lines)
    projected_vecs = lines .- (projs ./ norm.(projs) .* vels)
    normals = vels .× projected_vecs
    normals = normals ./ norm.(normals)
    # lines = [ wind_axes(horseshoe.bound_leg, uniform) for horseshoe in horseshoes ]

    # Compute matrices
    # dihedrals = trefftz_dihedral.(horseshoes)[:]
    @timeit "Dihedrals" dihedrals = [ atan(v[3], v[2]) for v in projected_vecs ]
    @timeit "Projected Leg Norms" Δs = norm.(projected_vecs)
    @timeit "Sum Γs" Δφs = sum(Γs, dims = 1)
    @timeit "Trefftz AIC" AIC = trefftz_matrix(horseshoes[end,:], normals, freestream / norm(freestream))

    # Solve system
    @timeit "∂φ/∂n" ∂φ_∂n = AIC * Δφs[:]

    # Compute forces
    pots_lens = Δφs .* Δs
    D_i = -0.5 * ρ * sum(pots_lens .* ∂φ_∂n)
    Y = - ρ * sum(pots_lens .* sin.(dihedrals))
    L = ρ * sum(pots_lens .* cos.(dihedrals))

    SVector(L, D_i, Y)
end 

"""
Evaluates the induced velocity by the trailing legs at the midpoint of a given Horseshoe `r`, by summing over the velocities of Horseshoes with vortex strengths `Γ`s, rotation rates `Ω`, and freestream flow vector `freestream` in the aircraft reference frame.
"""
midpoint_velocity(r, Ω :: SVector{3, Float64}, horseshoes :: Array{Horseshoe}, Γs :: Array{<: Real}, freestream) = sum(trailing_velocity(r, horseshoe, Γ, freestream / norm(freestream)) for (horseshoe, Γ) ∈ zip(horseshoes, Γs)) .- freestream .- Ω × r


"""
Evaluates the total induced velocity at a point `r` given Horseshoes, vortex strengths `Γ`s, rotation rates `Ω`, and freestream flow vector `freestream` in the global reference frame.
"""
stream_velocity(r, Ω :: SVector{3, Float64}, horseshoes :: Array{Horseshoe}, Γs :: Array{<: Real}, freestream) = sum(velocity(r, horseshoe, Γ, freestream / norm(freestream)) for (horseshoe, Γ) ∈ zip(horseshoes, Γs)) .+ freestream .+ Ω × r

"""
Computes the nearfield forces via the local Kutta-Jowkowski given an array of horseshoes, their associated vortex strengths Γs, a Uniform3D, and a density ρ. The velocities are evaluated at the midpoint of the bound leg of each horseshoe.
"""
function nearfield_forces(Γs :: Array{<: Real}, horseshoes :: Array{Horseshoe}, freestream, Ω, ρ)
    Γ_shoes = zip(Γs, horseshoes)
    @timeit "Summing Forces" [ 
      let r_i = bound_leg_center(horseshoe_i), l_i = bound_leg_vector(horseshoe_i); 
      ρ * midpoint_velocity(r_i, Ω, horseshoes, Γs, freestream) × l_i * Γ_focus end 
      for (Γ_focus, horseshoe_i) ∈ Γ_shoes 
    ]
end

"""
Computes the dynamic pressure given density ρ and speed V.
"""
dynamic_pressure(ρ, V) = 0.5 * ρ * V^2

"""
Placeholder. Unsure whether to change this to a generic moment computation function.
"""
moments(horseshoes :: Array{Horseshoe}, forces, r_ref) = [ (bound_leg_center(vortex_ring) .- r_ref) × force for (force, vortex_ring) ∈ zip(forces, horseshoes) ]

"""
Computes the non-dimensional force coefficient corresponding to standard aerodynamics.
"""
force_coefficient(force, q, S) = force / (q * S)

"""
Computes the non-dimensional moment coefficient corresponding to standard flight dynamics.
"""
moment_coefficient(moment, q, S, ref_length) = force_coefficient(moment, q, S) / ref_length

"""
Computes the non-dimensional angular velocity coefficient corresponding to standard flight dynamics.
"""
rate_coefficient(angular_speed, V, ref_length) = angular_speed * ref_length / 2V

"""
Computes the Kutta-Jowkowski forces and associated moments for horseshoes.
"""
function nearfield_dynamics(Γs :: Array{<: Real}, horseshoes :: Array{<: AbstractVortexArray}, freestream, Ω, r_ref, ρ = 1.225)
    @timeit "Forces" geom_forces = nearfield_forces(Γs, horseshoes, freestream, Ω, ρ)

    @timeit "Moments" geom_moments = moments(horseshoes, geom_forces, r_ref)

    geom_forces, geom_moments
    # , trefftz_force
end

function farfield_dynamics(Γs :: Array{<: Real}, horseshoes :: Array{<: AbstractVortexArray}, freestream, r_ref, ρ = 1.225)
        
    @timeit "Trefftz Force" trefftz_force = trefftz_forces(Γs, horseshoes, freestream, ρ)

    @timeit "Trefftz Moment" trefftz_moment = r_ref × trefftz_force

    trefftz_force, trefftz_moment
end

"""
Computes the near-field drag given the sum of the local Kutta-Jowkowski forces and a Uniform3D.
"""
nearfield_drag(force, uniform :: Uniform3D) = - dot(force, velocity(uniform) / uniform.mag)


# Post-processing
#==========================================================================================#

"""
Computes the streamlines from a given starting point, a Uniform3D, Horseshoes and their associated strengths Γs. The length of the streamline and the number of evaluation points must also be specified.
"""
function streamlines(point, uniform :: Uniform3D, Ω, horseshoes, Γs, length, num_steps)
    streamlines = [point]
    vel = velocity(uniform)
    @timeit "Iterating" for i ∈ 1:num_steps
        @timeit "Updating Velocity" update = stream_velocity(streamlines[end], Ω, horseshoes, Γs, -vel)
        @timeit "Adding Streamline" streamlines = [ streamlines..., streamlines[end] .+ (update / norm(update) * length / num_steps) ]
    end
    streamlines
end

"""
Computes the streamlines from the collocation points of given Horseshoes with the relevant previous inputs.
"""
streamlines(uniform :: Uniform3D, Ω, horseshoe_panels :: Array{<: Panel}, horseshoes :: Array{Horseshoe}, Γs :: Array{<: Real}, length :: Real, num_steps :: Integer) = [ streamlines(SVector(hs), uniform, Ω, horseshoes, Γs, length, num_steps) for hs ∈ collocation_point.(horseshoe_panels)[:] ]

"""
Prints the relevant aerodynamic/flight dynamics information.
"""
function aerodynamic_coefficients(force, moment, Ω, V, S, b, c, ρ)
    q = dynamic_pressure(ρ, V)
    CDi = force_coefficient(force[1], q, S)
    CY = force_coefficient(force[2], q, S)
    CL = force_coefficient(force[3], q, S)

    Cl = moment_coefficient(moment[1], q, S, b)
    Cm = moment_coefficient(moment[2], q, S, c)
    Cn = moment_coefficient(moment[3], q, S, b)

    p̄ = rate_coefficient(Ω[1], V, b)
    q̄ = rate_coefficient(Ω[2], V, c)
    r̄ = rate_coefficient(Ω[3], V, b)

    CL, CDi, CY, Cl, Cm, Cn, p̄, q̄, r̄
end

"""
Prints the relevant aerodynamic/flight dynamics information with a separate entry for drag.
"""
function aerodynamic_coefficients(force, moment, drag, Ω, V, S, b, c, ρ)
    q = dynamic_pressure(ρ, V)
    CDi = force_coefficient(drag, q, S)
    CY = force_coefficient(force[2], q, S)
    CL = force_coefficient(force[3], q, S)

    Cl = moment_coefficient(moment[1], q, S, b)
    Cm = moment_coefficient(moment[2], q, S, c)
    Cn = moment_coefficient(moment[3], q, S, b)

    p̄ = rate_coefficient(Ω[1], V, b)
    q̄ = rate_coefficient(Ω[2], V, c)
    r̄ = rate_coefficient(Ω[3], V, b)

    CL, CDi, CY, Cl, Cm, Cn, p̄, q̄, r̄
end

function print_dynamics(CL, CDi, CY, Cl, Cm, Cn, p̄, q̄, r̄)
    # println("Lift Coefficient")
    println("CL: $CL")
    # println("Drag Coefficient")
    println("CDi: $CDi")
    # println("Side Force Coefficient")
    println("CY: $CY")
    # println("Lift-to-Drag Ratio")
    println("L/D: $(CL/CDi)")
    # println("Rolling Moment Coefficient")
    println("Cl: $Cl")
    # println("Pitching Moment Coefficient")
    println("Cm: $Cm")
    # println("Yawing Moment Coefficient")
    println("Cn: $Cn")
    # println("Rolling Rate Coefficient")
    println("p̄: $p̄")
    # println("Pitching Rate Coefficient")
    println("q̄: $q̄")
    # println("Yawing Rate Coefficient")
    println("r̄: $r̄")
end

end