# Body source panel (skin-panelled fuselage coupling)
#==========================================================================================#
#
# A constant-strength source panel with a Neumann (no-penetration) boundary condition — the
# classic Hess–Smith non-lifting body element. Unlike the slender-body `FuselageLine`, this is
# a real 3-D skin panel: its unknown is the panel source strength `σ`, and its row in the
# shared AIC is the standard `velocity·normal` no-penetration row, so it couples monolithically
# with the lifting-surface vortices with no bespoke solve and no row overwrite.
#
# The kernel is the constant-source quadrilateral velocity of Hess–Smith (Katz & Plotkin,
# Low-Speed Aerodynamics §10.4.1). It is defined locally here (self-contained, mirroring
# `fuselage_line.jl`'s point kernels) rather than reusing the DoubletSource *potential* kernel.

"""
    SourcePanel(; min_area = 1e-8)

Element-model specification selecting constant-strength source-panel (Hess–Smith non-lifting
body) discretization. Pass it to [`elements`](@ref) to build `SourcePanel3D` elements from a
`Panel3D` skin mesh or a `HyperEllipseFuselage`, e.g. `elements(fuse, SourcePanel())`. Panels
with area below `min_area` (the collapsed nose/tail rings) are dropped.
"""
struct SourcePanel <: AbstractElementModel
    min_area::Float64
end

SourcePanel(; min_area = 1e-8) = SourcePanel(min_area)

"""
    SourcePanel3D(p1, p2, p3, p4, rc, normal, area, core)

A constant-strength source panel modelling one quadrilateral of a fuselage skin, coupled into
a `VortexLatticeSystem`. `p1…p4` are the panel corners (wound so `cross(p3-p1, p4-p2)` points
outward), `rc` is the centroid (control point), `normal` the **unit outward** normal, `area`
the panel area, and `core` a finite-core size (unused, kept for interface symmetry). The
unknown solved in the AIC is the panel source strength.

Assemble the returned vector into the aircraft as a `:body` block, e.g.
`ComponentVector(wing = make_horseshoes(mesh), body = make_body_panels(fuse))`.
"""
struct SourcePanel3D{T} <: AbstractVortex
    p1     :: SVector{3,T}
    p2     :: SVector{3,T}
    p3     :: SVector{3,T}
    p4     :: SVector{3,T}
    rc     :: SVector{3,T}   # Centroid (control point)
    normal :: SVector{3,T}   # Unit outward normal
    area   :: T
    core   :: T
end

function SourcePanel3D(p1, p2, p3, p4, rc, normal, area, core = zero(eltype(p1)))
    T = promote_type(eltype(p1), eltype(p2), eltype(p3), eltype(p4), eltype(rc), eltype(normal), typeof(area), typeof(core))
    SourcePanel3D{T}(p1, p2, p3, p4, rc, normal, area, core)
end

Base.length(::SourcePanel3D) = 1

# `control_point`/`normal_vector` are inherited from the AbstractVortex accessors (fields rc/normal).

## Constant-source quadrilateral velocity kernel (Hess–Smith)
#==========================================================================================#

# Panel-local orthonormal frame `[ŝ l̂ n̂]` (columns) with `n̂` the stored outward normal and the
# origin at the centroid. The in-plane axis `ŝ` is a corner-based tangent projected onto the
# panel plane, so the frame stays consistent with the stored normal after wind-axis rotation and
# Prandtl–Glauert scaling (which transform corners and normal by different maps).
function _source_panel_frame(el::SourcePanel3D)
    n  = el.normal
    s0 = (el.p2 + el.p3) / 2 - el.rc
    s  = normalize(s0 - dot(s0, n) * n)
    l  = cross(n, s)
    return el.rc, hcat(s, l, n)
end

"""
    quadrilateral_source_velocity(σ, el :: SourcePanel3D, r)

Induced velocity at `r` of a constant-strength source panel `el` of strength `σ`
(Hess–Smith constant-source quadrilateral, Katz & Plotkin §10.4.1). The point is transformed
to the panel-local frame, the four edge log/atan velocity terms are summed, and the result is
rotated back to global coordinates. At the panel's own centroid the field reduces to the
self-induced normal jump `σ/2 · n̂` (outward), giving a diagonal influence coefficient of `½`.
"""
function quadrilateral_source_velocity(σ, el::SourcePanel3D, r)
    c, R = _source_panel_frame(el)

    # Corners and evaluation point in the panel-local frame (origin at centroid, z ∥ n̂)
    q1, q2, q3, q4 = R' * (el.p1 - c), R' * (el.p2 - c), R' * (el.p3 - c), R' * (el.p4 - c)
    corners = (q1, q2, q3, q4)
    p = R' * (r - c)
    x, y, z = p

    T = eltype(p)
    ε = 1e-10

    # Self term: at the centroid the panel induces only the ±½ normal jump on the outward side.
    if abs(z) <= ε && abs(x) <= ε && abs(y) <= ε
        return σ / 2 * el.normal
    end

    u = zero(T); v = zero(T); w = zero(T)
    onplane = abs(z) <= ε # Coplanar off-panel points feel no normal (w) velocity.

    for i in 1:4
        j = i % 4 + 1
        xi, yi, _ = corners[i]
        xj, yj, _ = corners[j]

        ri  = sqrt((x - xi)^2 + (y - yi)^2 + z^2)
        rj  = sqrt((x - xj)^2 + (y - yj)^2 + z^2)
        dij = sqrt((xj - xi)^2 + (yj - yi)^2)

        # A collapsed edge (dij ≈ 0, e.g. the nose/tail cap triangles) contributes nothing to
        # either term; skipping it also avoids the 0/0 slope `mij` in the atan term below.
        dij > ε || continue

        # Tangential (in-plane) log terms
        lij = log((ri + rj - dij) / (ri + rj + dij))
        u += (yj - yi) / dij * lij
        v += (xi - xj) / dij * lij

        # Normal (out-of-plane) atan terms. The textbook per-edge form
        # `atan((mᵢⱼeᵢ-hᵢ)/(z rᵢ)) - atan((mᵢⱼeⱼ-hⱼ)/(z rⱼ))` divides by the local edge slope
        # `mᵢⱼ=(yⱼ-yᵢ)/(xⱼ-xᵢ)`, which is finite in value (atan saturates) but NaN under AD for a
        # near-vertical edge. The arctan-subtraction identity `atan p - atan q = atan2(p-q, 1+pq)`
        # gives the identical angle as a single division-free `atan2` of polynomials, so the slope
        # cancels and the whole kernel is differentiable.
        if !onplane
            dx = xj - xi
            dy = yj - yi
            ei = (x - xi)^2 + z^2
            hi = (x - xi) * (y - yi)
            ej = (x - xj)^2 + z^2
            hj = (x - xj) * (y - yj)
            Ni = dy * ei - dx * hi
            Nj = dy * ej - dx * hj
            Di = dx * z * ri
            Dj = dx * z * rj
            w += atan(Ni * Dj - Nj * Di, Di * Dj + Ni * Nj)
        end
    end

    # Negated so a positive source blows fluid outward (+n̂ side), matching the +σ/2 self term.
    return -R * (σ / 4π * SVector(u, v, w))
end

"""
    velocity(r, el :: SourcePanel3D, σ, V̂ = x̂)

Induced velocity at `r` of the body source panel `el` of strength `σ`. The trailing-direction
argument `V̂` is ignored (a source panel has no wake).
"""
velocity(r, el::SourcePanel3D, σ, V_hat = SVector{3,promote_type(eltype(r), typeof(σ))}(1, 0, 0)) = quadrilateral_source_velocity(σ, el, r)

## Nearfield force hooks (no-ops)
#==========================================================================================#
# A source panel carries no bound-leg circulation, so the Kutta–Joukowsky nearfield loop must
# ignore it: a zero bound-leg vector makes `kutta_joukowsky` vanish, and it induces no trailing
# (wake) velocity. Its pressure force is integrated separately (see `body_forces`).

bound_leg_center(el::SourcePanel3D) = el.rc
bound_leg_vector(el::SourcePanel3D{T}) where T = zero(SVector{3,T})
trailing_velocity(r, el::SourcePanel3D, σ, V) = zero(SVector{3, promote_type(eltype(r), typeof(σ))})

# No trailing wake, so it contributes nothing to the Trefftz-plane farfield.
has_wake(::SourcePanel3D) = false

## Axis transforms (wind-axis rotation and Prandtl-Glauert scaling)
#==========================================================================================#

transform(el::SourcePanel3D, T::LinearMap) = setproperties(el,
    p1     = T(el.p1),
    p2     = T(el.p2),
    p3     = T(el.p3),
    p4     = T(el.p4),
    rc     = T(el.rc),
    normal = normalize(T(el.normal)),
)

prandtl_glauert_scale_coordinates(el::SourcePanel3D, β) = setproperties(el,
    p1     = prandtl_glauert_scale_coordinates(el.p1, β),
    p2     = prandtl_glauert_scale_coordinates(el.p2, β),
    p3     = prandtl_glauert_scale_coordinates(el.p3, β),
    p4     = prandtl_glauert_scale_coordinates(el.p4, β),
    rc     = prandtl_glauert_scale_coordinates(el.rc, β),
    normal = normalize(prandtl_glauert_scale_normal(el.normal, β)),
)
