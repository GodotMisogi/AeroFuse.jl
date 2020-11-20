module VortexLattice

include("MathTools.jl")

using LinearAlgebra
using StaticArrays
using Rotations
using TimerOutputs
using .MathTools: accumap, structtolist

export Laplace, Uniform3D, Freestream, velocity, aircraft_velocity,
Panel, Panel3D, area, panel_coords, midpoint, panel_normal,
stability_axes, wind_axes,
solve_horseshoes, nearfield_dynamics, farfield_dynamics, nearfield_drag, 
streamlines, aerodynamic_coefficients, print_dynamics, horseshoe_lines

"""
Solutions to Laplace's equation.
"""
abstract type Laplace end

"""
A Uniform3D type expressing a freestream flow in spherical polar coordinates with angles α, β and magnitude mag.
"""
struct Uniform3D <: Laplace
    mag :: Real
    α :: Real 
    β :: Real
end

"""
Uniform3D constructor in degrees.
"""
Freestream(mag, α_deg, β_deg) = Uniform3D(mag, deg2rad(α_deg), deg2rad(β_deg))

"""
Converts freestream flow coordinates to Cartesian coordinates.
"""
freestream_to_cartesian(r, θ, φ) = r .* SVector(cos(θ) * cos(φ), -sin(φ), sin(θ) * cos(φ))

"""
Computes the velocity of Uniform3D.
"""
velocity(uni :: Uniform3D) = freestream_to_cartesian(uni.mag, uni.α, uni.β)

"""
Computes the velocity of Uniform3D in the aircraft reference frame.
"""
aircraft_velocity(uni :: Uniform3D) = -velocity(uni)


# Axis transformations
#==========================================================================================#

"""
Converts coordinates into stability axes.
"""
stability_axes(coords, uni :: Uniform3D) =  
                                            RotY{Float64}(uni.α) * coords
                                            # [cos(uni.α) 0 -sin(uni.α); 
                                            #       0     1     0     ;
                                            #  sin(uni.α) 0 sin(uni.α)] * coords
                                            

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
stability_axes(force, moment, Ω, freestream :: Uniform3D) = stability_axes(force, freestream), stability_axes(stab_flip(moment), freestream), stability_axes(stab_flip(Ω), freestream)

"""
Transforms forces and moments into wind axes.
"""
wind_axes(force, moment, freestream :: Uniform3D) = wind_axes(force, freestream), wind_axes([ -moment[1], moment[2], -moment[3] ], freestream)


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
    p1 :: SVector{3, Float64}
    p2 :: SVector{3, Float64}
    p3 :: SVector{3, Float64}
    p4 :: SVector{3, Float64}
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
                                     p31_x_p42 end

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
    r1 :: SVector{3, Float64}
    r2 :: SVector{3, Float64}
end

vector(line :: Line) = line.r2 - line.r1
center(line :: Line) = (line.r1 .+ line.r2) / 2

transform(line :: Line, rotation, translation) = Line(rotation * line.r1 + translation, rotation * line.r2 + translation)

"""
Helper function to compute the velocity induced by a bound vortex leg.
"""
bound_leg_velocity(a, b, Γ) = Γ/4π * (1/norm(a) + 1/norm(b)) * a × b / (norm(a) * norm(b) + dot(a, b))

"""
Helper function to compute the velocity induced by trailing vortex legs.
"""
trailing_legs_velocities(a, b, Γ, û) = Γ/4π * ((a × û) / (norm(a) - dot(a, û)) / norm(a) - (b × û) / (norm(b) - dot(b, û)) / norm(b))

"""
Helper function to check if any point is on the line.
"""
on_line(a, b, ε) = any(<(ε), norm.([a, b, a × b]))

"""
Computes the velocity induced at a point `r` by a vortex Line with constant strength Γ. Checks if `r` lies on the line and sets it to `(0, 0, 0)` if so, as the velocity is singular there.
"""
function horseshoe_velocity(r, line :: Line, Γ, ε = 1e-8; direction = SVector(1, 0, 0))
    a, b = r .- line.r1, r .- line.r2

    # Compute bound leg velocity
    @timeit "Bound Leg" bound_velocity = on_line(a, b, ε) ? zeros(3) : bound_leg_velocity(a, b, Γ)
    
    # Compute velocities of trailing legs
    @timeit "Trailing Leg" trailing_velocity = trailing_legs_velocities(a, b, Γ, direction)

    # Sums bound and trailing legs' velocities
    @timeit "Sum Legs" bound_velocity .+ trailing_velocity
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
Computes the midpoint of the bound leg of a Horseshoe or Vortex Ring.
"""
bound_leg_center(vortex :: AbstractVortexArray) = center(vortex.bound_leg)

"""
Computes the direction vector of the bound leg of a Horseshoe or Vortex Ring.
"""
bound_leg_vector(vortex :: AbstractVortexArray) = vector(vortex.bound_leg)


# Matrix setup
#==========================================================================================#

"""
Computes the velocity using the method of images for a symmetric case in the x-z plane.
"""
function mirror_velocity(collocation_point, horseshoe, Γ, Û)
    mirror_point = xz_flip(collocation_point)
    mir_vel = (xz_flip ∘ velocity)(mirror_point, horseshoe, Γ, Û)
end

"""
Computes the induced velocities at a point `r` of a Horseshoe with constant strength Γ and trailing legs pointing in a given direction.
"""
velocity(r, horseshoe :: Horseshoe, Γ, V̂) = horseshoe_velocity(r, horseshoe.bound_leg, Γ, direction = V̂)

"""
Computes the influence coefficient of the velocity of a vortex line at a collocation point projected to a normal vector.
"""
function influence_coefficient(collocation_point, horseshoe :: Horseshoe, panel_normal, V̂, symmetry = false) 
    if symmetry
        col_vel = velocity(collocation_point, horseshoe, 1., V̂)
        mir_vel = mirror_velocity(collocation_point, horseshoe, 1., V̂)
        return dot(col_vel .+ mir_vel, panel_normal)
    else 
        return dot(velocity(collocation_point, horseshoe, 1., V̂), panel_normal)
    end
end

"""
Assembles the Aerodynamic Influence Coefficient (AIC) matrix given vortex lines, collocation points and associated normal vectors.
"""
influence_matrix(colpoints, normals, vortex_lines :: Array{Horseshoe}, V̂, symmetry = false) = @timeit "Matrix Construction" [ @timeit "Influence Coefficient" influence_coefficient(colpoint_i, horsie_j, normal_i, V̂, symmetry) for (colpoint_i, normal_i) ∈ zip(colpoints, normals), horsie_j ∈ vortex_lines ]

"""
Computes the projection of a velocity vector with respect to normal vectors of panels. Corresponds to construction of the boundary condition for the RHS of the AIC system.
"""
boundary_condition(velocities, normals) = @timeit "Dotting" dot.(velocities, normals)

"""
Generates a wake Panel3D aligned with a normalised direction vector.
"""
function wake_panel(panel :: Panel3D, freestream :: Uniform3D, bound = 1.)
    wake = bound .* (normalize ∘ velocity)(freestream)
    Panel3D(panel.p2, panel.p2 .+ wake, panel.p3 .+ wake, panel.p3)
end

"""
Generates wake Panel3Ds given the Panel3Ds of a wing corresponding to the trailing edge and a Uniform3D. The wake length and number of wake panels must be provided.
"""
make_wake(last_panels :: Array{Panel3D}, freestream :: Uniform3D, wake_length, wake_num) = (permutedims ∘ accumap)(panel -> wake_panel(panel, freestream, wake_length / wake_num), wake_num, last_panels)

"""
Solves the AIC matrix with the boundary condition given Panel3Ds and a freestream velocity.
"""
function solve_horseshoes(horseshoe_panels :: Array{Panel3D}, camber_panels :: Array{Panel3D}, uniform :: Uniform3D, Ω, symmetry = false) 
    U = aircraft_velocity(uniform)
    @timeit "Horseshoes" horseshoes = horseshoe_lines.(horseshoe_panels)
    @timeit "Collocation Points" colpoints = collocation_point.(horseshoe_panels)
    @timeit "Normals" normals = panel_normal.(camber_panels)
    @timeit "Total Velocity" total_vel = [ U .+ (Ω × rc_i) for rc_i ∈ colpoints ]

    @timeit "AIC" AIC = influence_matrix(colpoints[:], normals[:], horseshoes[:], -normalize(U), symmetry) 
    @timeit "RHS" boco = boundary_condition(normals[:], total_vel[:])
    @timeit "Solve AIC" Γs = AIC \ boco

    @timeit "Reshape" output = reshape(Γs, size(horseshoes)...), horseshoes

    output
end


# Force evaluations
#==========================================================================================#

"""
Placeholder. Unsure whether to change this to a generic moment computation function.
"""
moments(horseshoes :: Array{Horseshoe}, forces, r_ref) = [ (bound_leg_center(vortex_ring) .- r_ref) × force for (force, vortex_ring) ∈ zip(forces, horseshoes) ]


#-----------------------------Nearfield evaluations-----------------------#

"""
Placeholder.
"""
bound_leg_velocity(r, line :: Line, Γ, ε = 1e-8) = let a = r .- line.r1, b = r .- line.r2; on_line(a, b, ε) ? zeros(3) : bound_leg_velocity(a, b, Γ) end

"""
Placeholder.
"""
trailing_legs_velocities(r, line :: Line, Γ; direction = SVector(1., 0., 0.)) = @timeit "Trailing Leg" let a = r .- line.r1, b = r .- line.r2; trailing_legs_velocities(a, b, Γ, direction) end

"""
Placeholder.
"""
trailing_velocity(r, horseshoe :: Horseshoe, Γ, V) = trailing_legs_velocities(r, horseshoe.bound_leg, Γ, direction = V)

"""
Evaluates the induced velocity by the trailing legs at the midpoint of a given Horseshoe `r`, by summing over the velocities of Horseshoes with vortex strengths `Γ`s, rotation rates `Ω`, and freestream flow vector `freestream` in the aircraft reference frame.
"""
midpoint_velocity(r :: SVector{3, Float64}, Ω :: SVector{3, Float64}, horseshoes :: Array{Horseshoe}, Γs :: Array{<: Real}, U) = sum(trailing_velocity(r, horseshoe, Γ, U) for (horseshoe, Γ) ∈ zip(horseshoes, Γs)) .- U .- Ω × r

"""
Computes the nearfield forces via the local Kutta-Jowkowski given an array of horseshoes, their associated vortex strengths Γs, a Uniform3D, and a density ρ. The velocities are evaluated at the midpoint of the bound leg of each horseshoe.
"""
function nearfield_forces(Γs :: Array{<: Real}, horseshoes :: Array{Horseshoe}, U, Ω, ρ)
    Γ_shoes = zip(Γs, horseshoes)
    @timeit "Summing Forces" [ 
      let r_i = bound_leg_center(horseshoe_i), l_i = bound_leg_vector(horseshoe_i); 
      ρ * midpoint_velocity(r_i, Ω, horseshoes, Γs, U) × l_i * Γ_i end 
      for (Γ_i, horseshoe_i) ∈ Γ_shoes 
    ]
end

"""
Computes the near-field drag given the sum of the local Kutta-Jowkowski forces and a Uniform3D.
"""
nearfield_drag(force, V :: Uniform3D) = dot(force, velocity(V) / V.mag)

"""
Computes the nearfield forces and associated moments.
"""
function nearfield_dynamics(Γs :: Array{<: Real}, horseshoes :: Array{Horseshoe}, uniform :: Uniform3D, Ω, r_ref, ρ = 1.225)
    U = aircraft_velocity(uniform)
    @timeit "Forces" geom_forces = nearfield_forces(Γs, horseshoes, U, Ω, ρ)
    @timeit "Moments" geom_moments = moments(horseshoes, geom_forces, r_ref)

    geom_forces, geom_moments
end

#---------------------------Farfield evaluations---------------------------#

trefftz_potential(r_i, r_j, Γ_j, Û) = let r = r_i .- r_j; Γ_j/2π * Û × r / dot(r, r) end

trefftz_matrix(trefftz_lines, normals, Û) = [ dot(trefftz_potential(center(tline_i), tline_j.r1, 1., Û), n̂_i) for (tline_i, n̂_i) in zip(trefftz_lines, normals), tline_j in trefftz_lines ]

trefftz_projection(line) = Line(SVector(0, line.r1[2], line.r1[3]), SVector(0, line.r2[2], line.r2[3]))

wind_axes(line :: Line, freestream :: Uniform3D) = Line(wind_axes(line.r1, freestream), wind_axes(line.r2, freestream)) 


"""
Computes the aerodynamic forces in the Trefftz plane normal to the freestream.
"""
function trefftz_forces(Γs, horseshoes :: Array{Horseshoe}, uniform :: Uniform3D, ρ)

    # lines = bound_leg_vector.(horseshoes[end,:])

    # projs = dot.(vels, lines)
    # projected_vecs = lines .- (projs ./ norm.(projs) .* vels

    # Project horseshoes' bound legs into Trefftz plane along freestream
    U = velocity(uniform)
    lines = [ wind_axes(horseshoe.bound_leg, uniform) for horseshoe in horseshoes[end,:] ][:]
    trefftz_lines = trefftz_projection.(lines)
    trefftz_vectors = vector.(trefftz_lines)

    Us = repeat([U], length(trefftz_lines))
    normals = Us .× trefftz_vectors
    normals = normals ./ norm.(normals)

    # Compute matrices
    @timeit "Dihedrals" dihedrals = [ atan(vec[3], vec[2]) for vec in trefftz_vectors ]
    @timeit "Projected Leg Norms" Δs = norm.(trefftz_vectors)
    @timeit "Sum Γs" Δφs = sum(Γs, dims = 1)
    @timeit "Trefftz AIC" AIC = trefftz_matrix(trefftz_lines, normals, normalize(U))

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
Compute farfield forces and moments.
"""
function farfield_dynamics(Γs :: Array{<: Real}, horseshoes :: Array{Horseshoe}, uniform :: Uniform3D, r_ref, ρ = 1.225)
    @timeit "Trefftz Force" trefftz_force = trefftz_forces(Γs, horseshoes, uniform, ρ)
    @timeit "Trefftz Moment" trefftz_moment = r_ref × trefftz_force

    trefftz_force, trefftz_moment
end


# Post-processing
#==========================================================================================#

"""
Evaluates the total induced velocity at a point `r` given Horseshoes, vortex strengths `Γ`s, rotation rates `Ω`, and freestream flow vector `freestream` in the global reference frame.
"""
stream_velocity(r, Ω :: SVector{3, Float64}, horseshoes :: Array{Horseshoe}, Γs :: Array{<: Real}, V) = sum(velocity(r, horseshoe, Γ, V / norm(V)) for (horseshoe, Γ) ∈ zip(horseshoes, Γs)) .+ freestream .+ Ω × r

"""
Computes the streamlines from a given starting point, a Uniform3D, Horseshoes and their associated strengths Γs. The length of the streamline and the number of evaluation points must also be specified.
"""
function streamlines(point, uniform :: Uniform3D, Ω, horseshoes, Γs, length, num_steps)
    streamlines = [point]
    V = velocity(uniform)
    @timeit "Iterating" for i ∈ 1:num_steps
        @timeit "Updating Velocity" update = stream_velocity(streamlines[end], Ω, horseshoes, Γs, V)
        @timeit "Adding Streamline" streamlines = [ streamlines..., streamlines[end] .+ (update / norm(update) * length / num_steps) ]
    end
    streamlines
end

"""
Computes the streamlines from the collocation points of given Horseshoes with the relevant previous inputs.
"""
streamlines(uniform :: Uniform3D, Ω, horseshoe_panels :: Array{<: Panel}, horseshoes :: Array{Horseshoe}, Γs :: Array{<: Real}, length :: Real, num_steps :: Integer) = [ streamlines(SVector(hs), uniform, Ω, horseshoes, Γs, length, num_steps) for hs ∈ collocation_point.(horseshoe_panels)[:] ]

end


# """
# A vortex ring type consisting of vortex lines. TODO: Consider if better fields are "bound_vortices" and "trailing_vortices"
# """
# struct VortexRing <: AbstractVortexArray
#     left_leg :: Line
#     bound_leg :: Line
#     back_leg :: Line
#     right_leg :: Line
# end

# """
# Helper function to compute the vortex ring given four points following Panel3D ordering.
# """
# function vortex_ring(p1, p2, p3, p4)
#     v1, v4 = quarter_point(p1, p2), quarter_point(p4, p3)
#     v2, v3 = v1 .+ p2 .- p1, v4 .+ p3 .- p4
#     [ v1, v2, v3, v4 ]
# end

# """
# Computes the vortex rings on a Panel3D.
# """
# vortex_ring(panel :: Panel3D) = vortex_ring(panel.p1, panel.p2, panel.p3, panel.p4)

# """
# Constructor for vortex rings on a Panel3D using Lines. The following convention is adopted:

# ``` 
#     p1 —bound_leg→ p4
#     |               |
# left_line       right_line
#     ↓               ↓
#     p2 —back_line→ p3
# ```
# """
# function VortexRing(panel :: Panel3D)
#     v1, v2, v3, v4 = vortex_ring(panel)

#     bound_leg = Line(v1, v4)
#     left_leg = Line(v2, v1)
#     right_leg = Line(v3, v2)
#     back_leg = Line(v4, v3)

#     VortexRing(left_leg, bound_leg, back_leg, right_leg)
# end

# """
# Computes the induced velocities at a point `r` of a Vortex Ring with constant strength Γ.
# """
# velocity(r, vortex_ring :: VortexRing, Γ) = sum_vortices(r, structtolist(vortex_ring), Γ)

# """
# Sums the velocities evaluated at a point `r` of vortex lines with constant strength Γ.
# """
# sum_vortices(r, vortex_lines :: Array{Line}, Γ) = sum(velocity(r, line, Γ) for line ∈ vortex_lines)