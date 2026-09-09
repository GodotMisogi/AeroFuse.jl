# Propeller slipstream (blown lift)
#==========================================================================================#
#
# A prescribed actuator-disk slipstream used to model blown/powered lift. The disk induces a
# uniform, fully-developed slipstream *tube* aligned with its thrust axis:
#
#   * Axial jump. From actuator-disk momentum theory the fully-developed axial speed inside
#     the tube is `V_s = √(V∞² + 2T/(ρ A))`, i.e. an increment `Δu = V_s - V∞` over the
#     freestream. The tube contracts from the disk radius `R` to the developed radius
#     `R_s = R √(V_disk / V_s)` with the disk-plane speed `V_disk = ½(V∞ + V_s)` (mass
#     continuity). The developed values are applied uniformly for all stations downstream of
#     the disk; the near-disk development region is not resolved.
#
#   * Swirl. The disk torque `Q` spins the slipstream. Conserving angular momentum in the
#     developed tube with a solid-body swirl `V_θ(r) = ω_s r` gives
#     `ω_s = 2Q / (π ρ V_s R_s⁴)`. `sense = ±1` sets the rotation direction (right-hand rule
#     about the thrust axis); `Q = 0` gives an axial-only slipstream.
#
#   * Turning (deflected slipstream / jet flap). A wing-flap immersed in the slipstream turns
#     the high-momentum jet. This is modelled by rotating the developed axial jet about a hinge
#     axis by up to `turn`, ramped smoothly (tanh) across the axial turning station `x_turn`
#     (e.g. the flap hinge): axial ahead of the station, turned behind it. The redirected axial
#     momentum in the lift direction is the powered-flap augmentation, which a purely axial
#     slipstream cannot produce. The tube centreline itself is left straight (only the velocity
#     vector is turned) and the swirl is unaffected.
#
# Unlike the fuselage line, the slipstream is a *prescribed* field (not an unknown): it enters
# the boundary condition at every collocation point (changing the circulations) and the local
# Kutta–Joukowsky velocity at every bound-leg centre (the dynamic-pressure boost inside the
# slipstream). The tube has a hard edge, so the field is discontinuous across `R_s`.

"""
    PropellerDisk(center, axis, radius, thrust, torque, sense;
                  hinge, turn, x_turn, turn_length)
    PropellerDisk(;
        center, axis,
        radius, thrust,
        torque = 0, sense = 1,
        hinge = [0, 1, 0], turn = 0, x_turn = 0, turn_length = radius,
    )

An actuator-disk propeller producing a prescribed slipstream for blown-lift analyses, defined
in **geometry axes**.

# Arguments
- `center :: SVector{3}`: Disk-centre location (m).
- `axis   :: SVector{3}`: Thrust direction; normalized on construction. The slipstream washes
  the region downstream of the disk along this axis.
- `radius :: Real`: Disk radius `R` (m).
- `thrust :: Real`: Disk thrust `T` (N), sets the axial slipstream increment.
- `torque :: Real`: Disk torque `Q` (N·m), sets the solid-body swirl; `0` for no swirl.
- `sense  :: Real`: Rotation sense `±1` about `axis` (right-hand rule).

Deflected-slipstream (jet-flap) turning, all defaulting to no turning:
- `hinge       :: SVector{3}`: Axis the jet turns about; normalized on construction. With the
  default spanwise `ŷ = [0,1,0]` a positive `turn` deflects a streamwise jet downward (lift).
- `turn        :: Real`: Maximum jet turning angle (radians), typically the flap deflection.
- `x_turn      :: Real`: Axial station of the turn (distance downstream of the disk along
  `axis`), e.g. the flap hinge location.
- `turn_length :: Real`: Streamwise ramp length over which the turn develops.
"""
struct PropellerDisk{T}
    center      :: SVector{3,T}
    axis        :: SVector{3,T}
    radius      :: T
    thrust      :: T
    torque      :: T
    sense       :: T
    hinge       :: SVector{3,T}
    turn        :: T
    x_turn      :: T
    turn_length :: T
end

# Raw field-wise constructor (all fields, in struct order): promotes to a common element type
# and stores verbatim. This is the constructor `setproperties`/`ConstructionBase` calls, so it
# must NOT normalize — the axis-transform helpers rely on it preserving the rotated vectors.
function PropellerDisk(center, axis, radius, thrust, torque, sense, hinge, turn, x_turn, turn_length)
    T = promote_type(eltype(center), eltype(axis), typeof(radius), typeof(thrust), typeof(torque), typeof(sense), eltype(hinge), typeof(turn), typeof(x_turn), typeof(turn_length))
    PropellerDisk{T}(SVector{3,T}(center), SVector{3,T}(axis), T(radius), T(thrust), T(torque), T(sense), SVector{3,T}(hinge), T(turn), T(x_turn), T(turn_length))
end

# Convenience constructor: normalizes the thrust and hinge directions, turning optional.
function PropellerDisk(center, axis, radius, thrust, torque = zero(radius), sense = one(radius);
                       hinge = SVector(0.0, 1.0, 0.0), turn = zero(radius), x_turn = zero(radius), turn_length = radius)
    PropellerDisk(center, axis / norm(axis), radius, thrust, torque, sense, hinge / norm(hinge), turn, x_turn, turn_length)
end

PropellerDisk(; center, axis, radius, thrust, torque = 0.0, sense = 1.0, hinge = SVector(0.0, 1.0, 0.0), turn = 0.0, x_turn = 0.0, turn_length = radius) =
    PropellerDisk(center, axis, radius, thrust, torque, sense; hinge, turn, x_turn, turn_length)

Base.broadcastable(disk :: PropellerDisk) = Ref(disk)

"""
    slipstream_velocity(r, disk :: PropellerDisk, refs :: References)

Induced velocity of the propeller slipstream at a point `r`, using the fully-developed
actuator-disk tube (uniform axial increment plus solid-body swirl). Returns the zero vector
ahead of the disk and outside the developed tube. See [`PropellerDisk`](@ref) for the model.
"""
function slipstream_velocity(r, disk :: PropellerDisk, refs :: References)
    T = promote_type(eltype(r), eltype(disk.center))
    z = zero(SVector{3,T})

    V∞ = refs.speed
    ρ  = refs.density
    A  = π * disk.radius^2

    d = SVector{3,T}(r) - disk.center
    x = dot(d, disk.axis)             # Axial distance downstream of the disk
    x <= 0 && return z                # No slipstream ahead of the disk

    rad = d - x * disk.axis           # Radial offset from the axis
    rr  = norm(rad)

    Vs = sqrt(V∞^2 + 2 * disk.thrust / (ρ * A)) # Fully-developed axial speed
    Δu = Vs - V∞
    Rs = disk.radius * sqrt((V∞ + Vs) / (2 * Vs)) # Developed tube radius (mass continuity)
    rr > Rs && return z               # Outside the slipstream tube

    # Deflected-slipstream turning: rotate the developed axial jet about the hinge by up to
    # `turn`, ramped smoothly across the turning station `x_turn`. The rotation preserves the
    # jet momentum magnitude and redirects it, so the lift-direction component is the
    # powered-flap augmentation. `AngleAxis` normalizes the hinge internally.
    θ       = disk.turn * (1 + tanh((x - disk.x_turn) / disk.turn_length)) / 2
    ax_turn = AngleAxis(θ, disk.hinge[1], disk.hinge[2], disk.hinge[3]) * disk.axis
    v_axial = Δu * ax_turn

    # Solid-body swirl from angular-momentum conservation; tangential direction by right-hand
    # rule about the thrust axis. Guard the axis (rr = 0) where swirl vanishes anyway.
    ω_s     = 2 * disk.torque / (π * ρ * Vs * Rs^4)
    t_hat   = rr > 0 ? disk.sense * cross(disk.axis, rad) / rr : z
    v_swirl = ω_s * rr * t_hat

    return v_axial + v_swirl
end

# Sum the slipstream of every disk in a collection.
slipstream_velocity(r, disks, refs :: References) = sum(disk -> slipstream_velocity(r, disk, refs), disks)

## Auto-derived turning from the immersed wing/flap
#==========================================================================================#

"""
    auto_turn(disk :: PropellerDisk, panels, refs :: References)

Return a copy of `disk` with its deflected-slipstream turning (`turn`, `x_turn`,
`turn_length`) derived from the wing `panels` immersed in the disk's slipstream tube, instead
of prescribed. The jet turning is the mean chordwise flow deflection of the trailing-edge
panels lying in the tube — i.e. the wing/flap turns the jet by the angle its trailing edge
deflects the flow, read from the panel normals — and it is ramped across the rear of the
immersed chord so the jet is fully turned by the trailing edge.

`panels` is a chordwise×spanwise array of `AbstractVortex` in geometry axes (e.g. from
`make_vortex_rings`/`make_horseshoes`), with the first chordwise index at the leading edge. If
no trailing-edge panel is immersed, the disk is returned unchanged.
"""
function auto_turn(disk :: PropellerDisk, panels, refs :: References)
    ĥ  = disk.hinge / norm(disk.hinge)
    up = cross(disk.axis, ĥ)
    up = up / norm(up) # Lift-direction reference: the normal of an undeflected panel

    # Developed tube radius (identical model to `slipstream_velocity`)
    V∞, ρ = refs.speed, refs.density
    Vs = sqrt(V∞^2 + 2 * disk.thrust / (ρ * π * disk.radius^2))
    Rs = disk.radius * sqrt((V∞ + Vs) / (2 * Vs))

    axial(r)  = dot(r - disk.center, disk.axis)
    radial(r) = norm((r - disk.center) - axial(r) * disk.axis)
    # Chordwise flow deflection encoded in a panel's normal (0 for an undeflected panel).
    deflection(vor) = atan(dot(normal_vector(vor), disk.axis), dot(normal_vector(vor), up))

    nc     = size(panels, 1)
    thresh = deg2rad(2.0) # Deflection rise (over the leading-edge baseline) that marks the hinge
    Σδ      = 0.0
    Σx_te   = 0.0
    Σx_hinge = 0.0
    count   = 0
    for j in axes(panels, 2)
        r_te = bound_leg_center(panels[nc, j])
        x    = axial(r_te)
        (x > 0 && radial(r_te) <= Rs) || continue # Trailing-edge panel must be in the tube

        # Hinge = first chordwise station deflected past the leading-edge baseline (the flap);
        # falls back to the trailing edge for an unflapped strip.
        base   = deflection(panels[1, j])
        i_h    = findfirst(i -> deflection(panels[i, j]) - base > thresh, 1:nc)
        Σx_hinge += axial(bound_leg_center(panels[isnothing(i_h) ? nc : i_h, j]))
        Σδ       += deflection(panels[nc, j]) # Turning set by the trailing-edge flow angle
        Σx_te    += x
        count    += 1
    end
    count == 0 && return disk

    turn        = Σδ / count
    x_te        = Σx_te / count
    x_hinge     = Σx_hinge / count
    flap_chord  = x_te - x_hinge
    turn_length = flap_chord > 0 ? flap_chord / 2 : disk.radius / 4
    x_turn      = (x_hinge + x_te) / 2 # Centre the ramp on the flap so the jet turns across it

    return setproperties(disk, turn = turn, x_turn = x_turn, turn_length = turn_length, hinge = ĥ)
end

## Axis transforms (wind-axis rotation and Prandtl-Glauert scaling), mirroring the aircraft
## so the slipstream can be evaluated in the frame the linear system is solved in.
#==========================================================================================#

transform(disk :: PropellerDisk, T :: LinearMap) = setproperties(disk,
    center = T(disk.center),
    axis   = T(disk.axis),
    hinge  = T(disk.hinge),
)

geometry_to_wind_axes(disk :: PropellerDisk, α, β) = transform(disk, LinearMap(RotZY(β, α)))

prandtl_glauert_scale_coordinates(disk :: PropellerDisk, β) = setproperties(disk,
    center = prandtl_glauert_scale_coordinates(disk.center, β),
    axis   = prandtl_glauert_scale_normal(disk.axis, β),
    hinge  = prandtl_glauert_scale_normal(disk.hinge, β),
)
