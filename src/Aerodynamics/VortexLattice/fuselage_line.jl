# Fuselage line singularity (slender-body coupling)
#==========================================================================================#

# A slender-body model of the fuselage as a line of singularities along its axis, coupled
# into the same linear system as the lifting-surface vortices. The physics splits in two:
#
#   * Thickness (displacement) is a *source* line whose strength is fixed by the
#     cross-sectional area distribution, `σ(x) = V∞ · dS/dx`. Because it is prescribed it is
#     NOT an unknown — it enters the boundary condition of every collocation point (wing,
#     tail and fuselage) as a known induced-velocity field (see `source_line_velocity`).
#
#   * Cross-flow lift is a *doublet* line whose strengths `λ` ARE the unknowns. Following
#     slender-body theory, each axial station behaves as a 2-D circular cross-section in the
#     local cross-flow `W`, so its doublet obeys the cylinder relation `λ = -2π R² W`, with
#     `W` the total vertical velocity there (freestream + wing-induced). That relation couples
#     the fuselage to the wing (wing → fuselage) while the doublet's 3-D field couples back
#     onto the wing (fuselage → wing). Cross-planes are otherwise independent, so the
#     fuselage-fuselage block of the system is the identity. This yields zero net body lift
#     but the destabilizing Munk pitching moment, exactly as slender-body theory predicts.
#
# Each element lumps one axis segment: a point source of strength `sigma` and a lumped 3-D
# cross-flow doublet (the unknown) at the segment midpoint on the axis.

"""
    FuselageLine(r1, r2, rc, normal, sigma, radius, core)

A slender-body fuselage segment. `r1`, `r2` are the axis segment endpoints and `rc` their
midpoint (the singularity location and cross-flow evaluation point); `normal` is the vertical
(lift / cross-flow) direction `ẑ`, which is also the doublet axis; `sigma` is the prescribed
source (thickness) strength `ΔS` per unit freestream speed; `radius` is the local body radius
`R` used in the 2-D cross-flow cylinder condition.
"""
struct FuselageLine{T} <: AbstractVortex
    r1     :: SVector{3,T}   # Axis segment start
    r2     :: SVector{3,T}   # Axis segment end
    rc     :: SVector{3,T}   # Axis midpoint (singularity + cross-flow evaluation point)
    normal :: SVector{3,T}   # Vertical (cross-flow / doublet) axis ẑ
    sigma  :: T              # Prescribed source (thickness) strength ΔS, per unit V∞
    radius :: T              # Local body radius R
    core   :: T              # Finite-core size (unused, kept for interface symmetry)
end

function FuselageLine(r1, r2, rc, normal, sigma, radius, core = zero(eltype(r1)))
    T = promote_type(eltype(r1), eltype(r2), eltype(rc), eltype(normal), typeof(sigma), typeof(radius), typeof(core))
    FuselageLine{T}(r1, r2, rc, normal, sigma, radius, core)
end

Base.length(::FuselageLine) = 1

# `control_point`/`normal_vector` are inherited from the AbstractVortex accessors (fields rc/normal).

# Axis segment length, used to lump the cross-flow doublet-line density into a point doublet.
segment_length(el::FuselageLine) = norm(el.r2 - el.r1)

## Velocity kernels
#==========================================================================================#

# Point-doublet induced velocity per unit moment `m` with axis `p̂` at `d = r - c`:
#   v = m/4π [ 3(p̂·d) d / ρ⁵ - p̂ / ρ³ ],  ρ = √(|d|² + ε²)
# Derived as ∇φ of the doublet potential φ = -m/4π (p̂·d)/ρ³. The finite core `ε` (set to the
# local body radius) regularizes the field: the singularity represents a body of radius `ε`,
# whose field is physically smooth outside `ε`. Without it, wing panels embedded in the
# fuselage (root chord inside the body radius) see an unphysically large 1/r³ near-field.
# Implemented locally rather than reusing `Laplace.Doublet3D`, whose kernel is unreliable.
function point_doublet_velocity(m, p, d, ε = 0.)
    ρ = sqrt(dot(d, d) + ε^2)
    return m / 4π * (3 * dot(p, d) * d / ρ^5 - p / ρ^3)
end

point_source_velocity(σ, d, ε = 0.) = σ / 4π * d / sqrt(dot(d, d) + ε^2)^3

"""
    velocity(r, el :: FuselageLine, λ, V̂ = x̂)

Induced velocity at `r` of the fuselage segment's cross-flow doublet (the lifting unknown)
with strength `λ`, using the 3-D doublet field of the lumped segment (moment `λ · Δx`). The
prescribed source part is excluded here so the AIC column stays linear in `λ`; the source
enters the boundary condition via [`source_line_velocity`].
"""
velocity(r, el::FuselageLine, λ, V_hat = SVector{3,promote_type(eltype(r), typeof(λ))}(1, 0, 0)) = point_doublet_velocity(λ * segment_length(el), el.normal, r - el.rc, el.radius)

"""
    source_line_velocity(r, fuse_elems, U)

Induced velocity at `r` of the prescribed fuselage *source* (thickness) line — the sum of
the lumped point sources of every `FuselageLine` in `fuse_elems`, scaled by the freestream
speed `norm(U)`. Used to build the extra boundary-condition velocity `Ups` injected at
every collocation point, which is how the fuselage displacement upwash reaches the wing.
"""
function source_line_velocity(r, fuse_elems, U)
    v = zero(SVector{3, promote_type(eltype(r), eltype(U))})
    for el in fuse_elems
        v += point_source_velocity(el.sigma, r - el.rc, el.radius)
    end
    return norm(U) * v
end

## Coupled solve with the 2-D cross-flow cylinder condition
#==========================================================================================#

"""
    solve_linear_fuselage(vortices, U, Ups, Ω)

Solve the coupled lifting-surface + fuselage system. The generic AIC and boundary condition
are assembled first (this gives the correct wing rows, including the fuselage doublet
columns). The fuselage rows are then overwritten with the slender-body 2-D cross-flow
cylinder condition `λ_i = -2π R_i² W_i`: each fuselage row equals `-2π R_i²` times its
generic (`velocity · ẑ`) row, and the fuselage-fuselage block is set to the identity because
slender-body cross-planes are independent.
"""
function solve_linear_fuselage(vortices, U, Ups, Ω)
    AIC  = influence_matrix(vortices)
    boco = boundary_condition(vortices, U, Ups, Ω)

    N        = length(vortices)
    fuse_col = [ vortices[j] isa FuselageLine for j in 1:N ]

    # Overwrite each fuselage row with the cylinder condition  λ_i = -2π R_i² W_i,
    # where W_i = (V∞ + Ω×r_c + source)·ẑ + (wing-induced)·ẑ is the external cross-flow.
    # Note V∞ = -U (the boundary condition stores U = -freestream), and the wing-induced
    # normal velocity is exactly the generic influence AIC[i, wing_j], so the wing columns are
    # scaled by +2π R_i² and moved to the left-hand side.
    @views for i in 1:N
        el = vortices[i]
        el isa FuselageLine || continue
        s  = 2π * el.radius^2
        rc = control_point(el)
        ni = normal_vector(el)
        for j in 1:N
            if fuse_col[j]
                AIC[i, j] = ifelse(i == j, one(eltype(AIC)), zero(eltype(AIC)))
            else
                AIC[i, j] *= s
            end
        end
        boco[i] = -s * dot(-U + Ω × rc + Ups[i], ni)
    end

    Γs = AIC \ boco

    return Γs, AIC, boco
end

## Nearfield force hooks (no-ops)
#==========================================================================================#
# The nearfield Kutta–Joukowsky loop maps over every element. A source/doublet line carries
# no bound-leg circulation, so it contributes zero force: a zero bound-leg vector makes
# `kutta_joukowsky` vanish, and it induces no trailing (wake) velocity on the surfaces.
# Proper fuselage forces are computed separately by slender-body integration.

bound_leg_center(el::FuselageLine) = el.rc
bound_leg_vector(el::FuselageLine{T}) where T = zero(SVector{3,T})
trailing_velocity(r, el::FuselageLine, Γ, V) = zero(SVector{3, promote_type(eltype(r), typeof(Γ))})

## Axis transforms (wind-axis rotation and Prandtl-Glauert scaling)
#==========================================================================================#

transform(el::FuselageLine, T::LinearMap) = setproperties(el,
    r1     = T(el.r1),
    r2     = T(el.r2),
    rc     = T(el.rc),
    normal = T(el.normal),
)

prandtl_glauert_scale_coordinates(el::FuselageLine, β) = setproperties(el,
    r1     = prandtl_glauert_scale_coordinates(el.r1, β),
    r2     = prandtl_glauert_scale_coordinates(el.r2, β),
    rc     = prandtl_glauert_scale_coordinates(el.rc, β),
    normal = prandtl_glauert_scale_normal(el.normal, β),
)
