function quarter_point(p1, p2) 
    μ = SVector(1/4, 0, 1/4)
    return @. (1 - μ) * p1 + μ * p2
end 

function three_quarter_point(p1, p2) 
    μ = SVector(3/4, 0, 3/4)
    return @. (1 - μ) * p1 + μ * p2
end 

control_point(p1, p2, p3, p4) = (three_quarter_point(p1, p2) + three_quarter_point(p4, p3)) / 2
bound_leg(p1, p2, p3, p4) = (quarter_point(p1, p2), quarter_point(p4, p3))

"""
    bound_leg(panel :: Panel3D)

Compute the bound leg for a `Panel3D`, for horseshoes/vortex rings, which quarter point on each side of the trailing legs.
"""
bound_leg(panel :: Panel3D) = bound_leg(panel.p1, panel.p2, panel.p3, panel.p4)

"""
    control_point(panel :: Panel3D)

Compute the control point of a `Panel3D` for horseshoes/vortex rings, which is the average of the 3-quarter point on each side of the trailing legs.
"""
control_point(panel :: Panel3D) = control_point(panel.p1, panel.p2, panel.p3, panel.p4)


"""
    HorseshoeVortex(panel :: Panel3D, normal, drift = zeros(3))

Generate a `HorseshoeVortex` corresponding to a `Panel3D`, an associated normal vector, and a "drift velocity".
"""
function HorseshoeVortex(panel :: Panel3D, normal, drift = @SVector zeros(3); core_size = 0.)
    r1, r2 = bound_leg(panel)
    rc = control_point(panel) + drift
    HorseshoeVortex(r1, r2, rc, normal, core_size)
end

"""
Constructor for making a `RingVortex` with a `Panel3D`. The following convention is adopted:

```
    p1 —front leg→ p4
    |               |
left leg       right leg
    ↓               ↓
    p2 —back leg-→ p3
```
"""
function RingVortex(panel :: Panel3D{T}, rc, normal, trailing = false; core_size = zero(T)) where T <: Number
    # r1 = quarter_point(panel.p1, panel.p2)
    # r4 = quarter_point(panel.p4, panel.p3)
    # r2 = normalize(panel.p2 - panel.p1) * 0.25 + r1
    # r3 = normalize(panel.p3 - panel.p4) * 0.25 + r2
    # rc = control_point(panel) # (r1 + r2 + r3 + r4) / 4
    RingVortex{T}(panel.p1, panel.p2, panel.p3, panel.p4, rc, normal, trailing, core_size)
end

## Generative element builder
#==========================================================================================#

"""
    elements(geometry, model :: AbstractElementModel)

Build the populated singularity elements for a component from its `geometry` and an
element-model specification. This is the generic entry point for assembling a component:
pass a mesh (or panels) together with the element type to use.

Supported combinations:
- `elements(wing :: WingMesh, :: Horseshoe)` → `HorseshoeVortex` elements from the chord
  panels and camber-surface normals.
- `elements(wing :: WingMesh, :: VortexRing)` → `RingVortex` lattice from the camber
  coordinates, with trailing-edge wake identification.

See [`Horseshoe`](@ref), [`VortexRing`](@ref) and [`SourcePanel`](@ref) for the models.
"""
elements(wing :: WingMesh, model :: Horseshoe) =
    map((cho, cam) -> HorseshoeVortex(cho, normal_vector(cam); core_size = model.core_size), chord_panels(wing), camber_panels(wing))

elements(wing :: WingMesh, model :: VortexRing) = make_vortex_rings(camber_coordinates(wing); core_size = model.core_size)

"""
    make_horseshoes(wing :: WingMesh)

Generate an array of `HorseshoeVortex` elements defined by the chord coordinates and normal
vectors of the camber distribution of a `WingMesh`. Convenience wrapper for
`elements(wing, Horseshoe())`.
"""
make_horseshoes(wing :: WingMesh) = elements(wing, Horseshoe())

"""
    make_vortex_rings(wing :: WingMesh)

Generate an array of `RingVortex` elements defined by the camber coordinates and normal
vectors of a `WingMesh`. Convenience wrapper for `elements(wing, VortexRing())`.
"""
make_vortex_rings(wing_mesh :: WingMesh) = elements(wing_mesh, VortexRing())

@views function make_vortex_rings(cam_coo; core_size = 0.)
    # Generate vortex ring mesh
    cams = combinedimsview(cam_coo, (1,2))
    vor_cams = similar(cams)
    vor_cams[1:end-1,:,:] = 0.75 * cams[1:end-1,:,:] + 0.25 * cams[2:end,:,:]
    vor_cams[end,:,:] = cams[end,:,:]

    # Construct vortex rings with trailing edge identification for boundary condition
    cam_pan = make_panels(cam_coo)
    vor_pans = make_panels(splitdimsview(vor_cams, (1,2)))
    rings = map(CartesianIndices(vor_pans)) do ind
        i, j = ind.I
        if i == size(vor_pans, 2)
            RingVortex(vor_pans[i,j], control_point(cam_pan[i,j]), normal_vector(cam_pan[i,j]); core_size = core_size)
        else
            RingVortex(vor_pans[i,j], control_point(cam_pan[i,j]), normal_vector(cam_pan[i,j]), true; core_size = core_size)
        end
    end

    return rings
end

"""
    deflect_normals(
        rings :: AbstractMatrix{<:RingVortex}, δ;
        hinge = 0.75, axis = SVector(0., 1., 0.), sense = _ -> 1,
    )

Emulate a trailing-edge control-surface deflection on a vortex-ring mesh by rotating the
normal vectors of the rings behind the normalized chordwise `hinge` location by the angle
`δ` (radians) about `axis`. Rings ahead of the hinge are returned unchanged.

`sense(j)` is an optional multiplier evaluated per spanwise column `j`, e.g. antisymmetric
`±1` for an aileron or `0` to exclude a station; it defaults to a symmetric deflection
(elevator/flap). The mesh is assumed to be indexed `(chordwise, spanwise)`.

Each affected ring is rebuilt through the promoting `RingVortex` constructor, so the mesh
adopts the element type of `δ` and the deflection is differentiable with respect to `δ`
(compatible with ForwardDiff and other AD backends).
"""
function deflect_normals(rings :: AbstractMatrix{<:RingVortex}, δ; hinge = 0.75, axis = SVector(0., 1., 0.), sense = _ -> 1)
    nc = size(rings, 1)
    i_hinge = floor(Int, hinge * nc)   # chordwise rows with index > i_hinge lie aft of the hinge
    return map(CartesianIndices(rings)) do idx
        i, j = Tuple(idx)
        ring = rings[idx]
        i > i_hinge || return ring
        n̂ = normalize(ring.normal + sense(j) * δ * cross(axis, ring.normal))
        RingVortex(ring.r1, ring.r2, ring.r3, ring.r4, ring.rc, n̂, ring.trailing, ring.core)
    end
end

reynolds_number(refs :: References) = refs.density * refs.speed * refs.chord / refs.viscosity

"""
    make_fuselage_line(fuse :: HyperEllipseFuselage; n = 20)

Generate an array of `FuselageLine` singularity segments modelling a `HyperEllipseFuselage`
as a slender body along its axis, for coupling into a `VortexLatticeSystem`. Each segment
carries a prescribed source (thickness) strength `ΔS = π ΔR²` from the cross-sectional area
distribution and an unknown doublet (cross-flow lift) whose strength solves in the AIC. `n`
sets the number of stations per section (nose, cabin, rear).

Assemble the returned vector into the aircraft as a `:fuse` block, e.g.
`ComponentVector(wing = make_horseshoes(mesh), fuse = make_fuselage_line(fuse))`.
"""
function make_fuselage_line(fuse :: HyperEllipseFuselage; n = 20)
    ts = LinRange(0., 1., n)

    # Radius distribution R(x) and axis x-stations in the local (pre-affine) frame
    xn, xc, xr, Rn, Rc, Rr = undrooped_curve(fuse, ts)
    xs = [ xn; xc; xr ]   # x-stations
    Rs = [ Rn; Rc; Rr ]   # radii R(x)

    # Centerline droop = top-surface profile (R + droop, from `curve`) − R
    z_cen = @views curve(fuse, ts)[:,2] .- Rs

    # Local direction → world direction (the affine translation cancels in the difference)
    aff = fuse.affine
    to_world_dir(d) = aff(SVector(d...)) - aff(SVector(0., 0., 0.))
    normal = normalize(to_world_dir(SVector(0., 0., 1.))) # Vertical (cross-flow / doublet) axis

    # Skip the zero-length segments at the nose/cabin and cabin/rear junctions (shared points)
    N    = length(xs)
    segs = [ i for i in 1:N-1 if xs[i+1] - xs[i] > 1e-9 ]

    return map(segs) do i
        # Axis segment endpoints (centerline) and midpoint, mapped to world coordinates
        r1 = aff(SVector(xs[i],   0., z_cen[i]))
        r2 = aff(SVector(xs[i+1], 0., z_cen[i+1]))
        rc = (r1 + r2) / 2

        # Prescribed source (thickness) strength ΔS = π ΔR², per unit freestream speed
        sigma = π * (Rs[i+1]^2 - Rs[i]^2)

        # Local body radius for the 2-D cross-flow cylinder condition
        radius = (Rs[i] + Rs[i+1]) / 2

        FuselageLine(r1, r2, rc, normal, sigma, radius)
    end
end

"""
    make_fuselage_panels(fuse :: HyperEllipseFuselage; n_secs = 20, n_circ = 20)

Panel the skin of a `HyperEllipseFuselage` into a `Matrix{Panel3D}` for a 3-D body-panel
analysis. `n_secs` sets the number of axial stations per section (nose, cabin, rear) and
`n_circ` the number of circumferential points per ring. The panel winding is oriented so each
panel's `normal_vector` points outward from the body axis (needed for a consistent Neumann
sign); near-degenerate nose/tail rings are kept and should be filtered by area downstream.

`n_circ` is rounded up to the next odd value so the circumferential distribution
`LinRange(0, 2π, n_circ)` is mirror-symmetric about the x-z plane. Without this the ±y panels
are not exact mirror pairs, and a sharp asymmetric pressure feature (e.g. a wing root near the
body) integrates to a spurious side force on an otherwise symmetric configuration.
"""
function make_fuselage_panels(fuse :: HyperEllipseFuselage; n_secs = 20, n_circ = 20)
    n_circ = isodd(n_circ) ? n_circ : n_circ + 1   # Enforce x-z mirror-symmetric rings
    ts   = LinRange(0., 1., n_secs)
    coo  = coordinates(fuse, ts, n_circ)   # (n_circ, 3·n_secs, 3) numeric skin grid
    xyzs = splitdimsview(coo, (1,2))       # (n_circ, 3·n_secs) grid of SVector{3}
    panels = make_panels(xyzs)

    # Ring centres per axial station, used to test the outward-normal orientation.
    ncirc, nsec = size(xyzs)
    centres = [ sum(@view xyzs[:,j]) / ncirc for j in 1:nsec ]

    # Test a mid-body panel; the mesh has uniform winding, so flip the whole mesh if inward.
    i0, j0  = max(size(panels, 1) ÷ 2, 1), max(size(panels, 2) ÷ 2, 1)
    pan     = panels[i0, j0]
    axis_pt = (centres[j0] + centres[j0 + 1]) / 2
    inward  = dot(normal_vector(pan), midpoint(pan) - axis_pt) < 0

    return inward ? map(p -> Panel3D(p.p1, p.p4, p.p3, p.p2), panels) : panels
end

"""
    elements(panels :: AbstractArray{<:Panel3D}, model :: SourcePanel)

Build a vector of `SourcePanel3D` constant-strength source panels from an existing `Panel3D`
skin mesh. Panels with area below `model.min_area` (the collapsed nose/tail rings) are
dropped. The panel winding is assumed already oriented outward (see [`make_fuselage_panels`]).
"""
function elements(panels :: AbstractArray{<:Panel3D}, model :: SourcePanel)
    good = filter(p -> panel_area(p) > model.min_area, vec(panels))

    return map(good) do pan
        SourcePanel3D(pan.p1, pan.p2, pan.p3, pan.p4, midpoint(pan), normalize(normal_vector(pan)), panel_area(pan))
    end
end

"""
    elements(fuse :: HyperEllipseFuselage, model :: SourcePanel; n_secs = 20, n_circ = 20)

Panel the skin of a `HyperEllipseFuselage` and build `SourcePanel3D` source panels from it for
monolithic coupling into a `VortexLatticeSystem`. `n_secs`/`n_circ` set the axial/circumferential
panel density (see [`make_fuselage_panels`]).
"""
elements(fuse :: HyperEllipseFuselage, model :: SourcePanel; n_secs = 20, n_circ = 20) =
    elements(make_fuselage_panels(fuse; n_secs, n_circ), model)

"""
    make_body_panels(fuse :: HyperEllipseFuselage; n_secs = 20, n_circ = 20, min_area = 1e-8)

Generate a vector of `SourcePanel3D` constant-strength source panels modelling the skin of a
`HyperEllipseFuselage`, for monolithic coupling into a `VortexLatticeSystem`. Convenience
wrapper for `elements(fuse, SourcePanel(; min_area); n_secs, n_circ)`.

Assemble the returned vector into the aircraft as a `:body` block, e.g.
`ComponentVector(wing = make_horseshoes(mesh), body = make_body_panels(fuse))`.
"""
make_body_panels(fuse :: HyperEllipseFuselage; n_secs = 20, n_circ = 20, min_area = 1e-8) =
    elements(fuse, SourcePanel(; min_area); n_secs, n_circ)