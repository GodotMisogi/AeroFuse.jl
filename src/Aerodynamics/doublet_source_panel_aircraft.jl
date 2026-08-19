# Generic multi-component 3D doublet-source (Morino) panel solver
#==========================================================================================#
#
# Assembles and solves an aircraft made of several lifting/closed surfaces (matrices of
# `Panel3D`) and, optionally, a slender-body fuselage line (`make_fuselage_line`), in a
# single coupled linear system — mirroring the `ComponentVector`-of-components style of the
# `VortexLatticeSystem`. Each surface contributes constant-strength doublet panels with a
# Morino/Kutta wake; the fuselage contributes a prescribed source (thickness) field plus
# unknown cross-flow doublets closed by the slender-body 2-D cylinder condition.
#
# The aircraft is a `NamedTuple` (surfaces are `Matrix{<:Panel3D}`, an optional `fuse` key
# is a `Vector{<:FuselageLine}`); the solution is packed into `ComponentVector`s keyed by
# component, so `system.doublets.wing`, `system.wake_strengths.htail`, etc. are accessible.

import ForwardDiff

# ---- point-singularity potentials (fuselage line coupling) ------------------------------
# The panel doublet kernel `quadrilateral_doublet_potential` is the NEGATIVE of the standard
# point-doublet potential (built so a closed body's influence rows sum to +1), so the
# fuselage doublet column is negated to match.
_point_doublet_potential(m, phat, d, ε = 0.0) = (ρ = sqrt(dot(d, d) + ε^2); -m / 4π * dot(phat, d) / ρ^3)
_point_source_potential(σ, d, ε = 0.0)        = (ρ = sqrt(dot(d, d) + ε^2); -σ / (4π * ρ))

# Solid-angle doublet potential (Van Oosterom–Strackee), used for the α-dependent wake block.
# The panel-local `quadrilateral_doublet_potential` is exact but not ForwardDiff-safe: an
# axis-aligned edge in the local frame drives `mij = (yj-yi)/(xj-xi) → Inf`, so `atan(Inf)`
# returns a correct finite `±π/2` but a `NaN` derivative. The solid-angle form (a smooth
# `atan2` with a strictly positive denominator) matches it to ~1e-15 and differentiates
# cleanly. A constant-strength doublet panel's potential is exactly `-μ/4π · Ω`, with `Ω` the
# signed solid angle it subtends at the field point.
function _tri_solid_angle(a, b, c)
    na, nb, nc = norm(a), norm(b), norm(c)
    num = dot(a, cross(b, c))
    den = na*nb*nc + dot(a, b)*nc + dot(a, c)*nb + dot(b, c)*na
    2 * atan(num, den)
end

function solidangle_doublet_potential(μ, pan :: AbstractPanel3D, r)
    a1, a2, a3, a4 = p1(pan), p2(pan), p3(pan), p4(pan)
    Ω = _tri_solid_angle(a1 - r, a2 - r, a3 - r) + _tri_solid_angle(a1 - r, a3 - r, a4 - r)
    -μ / 4π * Ω
end

# A constant-strength doublet panel is equivalent to a vortex ring of circulation Γ = μ
# around its four edges; its induced velocity is the Biot–Savart sum of the edges.
function panel_ring_velocity(pan :: AbstractPanel3D, r)
    a1, a2, a3, a4 = p1(pan), p2(pan), p3(pan), p4(pan)
    bound_leg_velocity(r - a1, r - a2, 1.0) + bound_leg_velocity(r - a2, r - a3, 1.0) +
    bound_leg_velocity(r - a3, r - a4, 1.0) + bound_leg_velocity(r - a4, r - a1, 1.0)
end

# ---- system type ------------------------------------------------------------------------
"""
    DoubletSourcePanelSystem

Result of a coupled 3D doublet-source panel analysis. Accessible fields:
- `surfaces`:       `NamedTuple` of the surface panel matrices (`Matrix{Panel3D}`).
- `wakes`:          `NamedTuple` of the shed wake panels per surface.
- `fuselage`:       the `Vector{FuselageLine}` slender body, or `nothing`.
- `doublets`:       `ComponentVector` of surface doublet strengths, keyed by component.
- `wake_strengths`: `ComponentVector` of wake doublet strengths (= trailing-edge jumps).
- `fuse_doublets`:  the fuselage cross-flow doublet strengths `λ`, or `nothing`.
- `influence_matrix`, `boundary_vector`: the assembled linear system.
- `freestream :: Freestream`, `reference :: References`.
"""
struct DoubletSourcePanelSystem{S,W,F,D,WS,FD,M,N,P,R}
    surfaces         :: S
    wakes            :: W
    fuselage         :: F
    doublets         :: D
    wake_strengths   :: WS
    fuse_doublets    :: FD
    influence_matrix :: M
    boundary_vector  :: N
    freestream       :: P
    reference        :: R
end

function Base.show(io :: IO, sys :: DoubletSourcePanelSystem)
    println(io, "DoubletSourcePanelSystem —")
    for (k, p) in pairs(sys.surfaces)
        println(io, "  ", rpad(k, 8), size(p), " ", eltype(p))
    end
    isnothing(sys.fuselage) || println(io, "  ", rpad(:fuse, 8), length(sys.fuselage), " FuselageLine")
    println(io, "  freestream α = ", rad2deg(sys.freestream.alpha), "°, β = ", rad2deg(sys.freestream.beta), "°")
end

_surfaces_only(aircraft :: NamedTuple) = Base.structdiff(aircraft, NamedTuple{(:fuse,)})
_panelview(p) = permutedims(p)[:]

"""
    solve_case(aircraft :: NamedTuple, fs :: Freestream, refs :: References; wake_length = 1e3)

Solve the coupled 3D doublet-source panel system for an `aircraft` — a `NamedTuple` whose
entries are surface panel matrices (`surface_panels(mesh)`), plus an optional `fuse` entry
built with [`make_fuselage_line`](@ref). Returns a [`DoubletSourcePanelSystem`](@ref).

```julia
aircraft = (
    wing  = surface_panels(wing_mesh),
    htail = surface_panels(htail_mesh),
    fuse  = make_fuselage_line(fuselage),
)
system = solve_case(aircraft, Freestream(alpha = 3.0), refs)
```
"""
function solve_case(aircraft :: NamedTuple, fs :: Freestream, refs :: References; wake_length = 1e3)
    V∞       = velocity(fs)
    surfaces = _surfaces_only(aircraft)
    fuse     = haskey(aircraft, :fuse) ? aircraft.fuse : nothing

    # Body panels and one wake panel per spanwise strip, organised per component.
    bodies = ComponentArray(map(_panelview, surfaces))
    wakes  = map(p -> [ wake_panel(p[:, j], wake_length, V∞) for j in axes(p, 2) ], surfaces)
    wakevec = ComponentArray(wakes)

    B  = getdata(bodies);  Nb = length(B)
    W  = getdata(wakevec); Nw = length(W)
    Nf = isnothing(fuse) ? 0 : length(fuse)
    N  = Nb + Nw + Nf

    # Eltype follows the freestream so ForwardDiff Duals (∂/∂α, ∂/∂β) propagate through the
    # whole assembly and linear solve.
    T = promote_type(eltype(V∞), Float64)

    # Doublet influence blocks (field points on rows).
    AIC = zeros(T, N, N)
    AIC[1:Nb, 1:Nb]         .= permutedims(doublet_matrix(B, B))   # body → body (α-independent)
    # Wake → body: the wake geometry depends on α (shed along V∞), so use the ForwardDiff-safe
    # solid-angle doublet kernel here (matches the panel-local kernel to ~1e-15).
    for jw in 1:Nw, ib in 1:Nb
        AIC[ib, Nb+jw] = solidangle_doublet_potential(1.0, W[jw], collocation_point(B[ib]))
    end
    boco = zeros(T, N)
    boco[1:Nb] .= [ dot(V∞, collocation_point(p)) for p in B ]     # Φ∞ (Morino RHS)

    # Morino–Kutta rows:  μ(first chord panel) − μ(last chord panel) + μ_wake = 0
    for s in keys(surfaces)
        nc, ns = size(surfaces[s])
        bidx = ComponentArrays.label2index(bodies, s)
        widx = ComponentArrays.label2index(wakevec, s)
        for j in 1:ns
            r = Nb + widx[j]
            AIC[r, bidx[j]]             += 1.0
            AIC[r, bidx[(nc-1)*ns + j]] -= 1.0
            AIC[r, Nb + widx[j]]        += 1.0
        end
    end

    # Fuselage slender-body coupling.
    if !isnothing(fuse)
        seglen = segment_length.(fuse)
        for i in 1:Nb
            ri = collocation_point(B[i])
            # Prescribed source (thickness) potential adds to the body boundary condition.
            boco[i] += sum(el -> _point_source_potential(el.sigma * norm(V∞), ri - el.rc, el.radius), fuse)
            # Unknown cross-flow doublet potential column (negated to match the panel kernel).
            for f in 1:Nf
                el = fuse[f]
                AIC[i, Nb+Nw+f] = -_point_doublet_potential(seglen[f], el.normal, ri - el.rc, el.radius)
            end
        end
        # Cylinder condition  λ_f = -2π R² W_f : the panel/wake-induced cross-flow moves to
        # the LHS, the freestream + prescribed-source cross-flow forms the RHS, and the
        # fuselage-fuselage block is the identity (slender-body cross-planes are independent).
        for f in 1:Nf
            el = fuse[f]; c = 2π * el.radius^2; n̂ = el.normal; rc = el.rc
            r = Nb + Nw + f
            for j in 1:Nb; AIC[r, j]      = c * dot(panel_ring_velocity(B[j], rc), n̂); end
            for w in 1:Nw; AIC[r, Nb + w] = c * dot(panel_ring_velocity(W[w], rc), n̂); end
            AIC[r, r] = 1.0
            Vsrc = sum(g -> g === el ? zero(SVector{3,Float64}) :
                            point_source_velocity(g.sigma * norm(V∞), rc - g.rc, g.radius), fuse)
            boco[r] = -c * dot(V∞ + Vsrc, n̂)
        end
    end

    x = AIC \ boco

    # Repack the solution by component.
    ks       = keys(surfaces)
    doublets = ComponentArray(NamedTuple{ks}(map(s -> x[ComponentArrays.label2index(bodies, s)], ks)))
    wake_str = ComponentArray(NamedTuple{ks}(map(s -> x[Nb .+ ComponentArrays.label2index(wakevec, s)], ks)))
    λ        = isnothing(fuse) ? nothing : x[Nb+Nw+1:end]

    return DoubletSourcePanelSystem(surfaces, wakes, fuse, doublets, wake_str, λ, AIC, boco, fs, refs)
end

# ---- post-processing --------------------------------------------------------------------
"""
    surface_coefficients(system :: DoubletSourcePanelSystem)

Pressure coefficient ``C_p`` on every panel of every surface, returned as a `NamedTuple` of
matrices matching `system.surfaces`. The surface velocity is the tangential gradient of the
doublet (total-potential) distribution, recovered per panel from a robust 2×2 local-frame
finite difference of the neighbouring doublet strengths.
"""
function surface_coefficients(system :: DoubletSourcePanelSystem)
    Vmag = norm(velocity(system.freestream))
    map(keys(system.surfaces)) do s
        _panel_pressures(system.surfaces[s], system.doublets[s], Vmag)
    end |> NamedTuple{keys(system.surfaces)}
end

# Per-panel surface velocity as GLOBAL 3-vectors. The Morino doublet strength is the total
# surface potential, so the tangential surface velocity is its surface gradient, recovered
# per panel from a robust 2×2 local-frame finite difference of neighbouring doublet strengths
# (robust to sweep) and rotated back to global axes. Normal component is zero by construction
# (flow tangency), so these vectors are the wetted-surface inviscid edge velocities directly
# usable as a boundary-layer EIF.
function _panel_velocities(p, μblock, Vmag)
    nc, ns = size(p)
    φs = permutedims(reshape(collect(μblock), ns, nc))
    tup(a, b) = (a, b)
    clp = collocation_point.(p)
    xpair = midpair_map(tup, clp; dims = 1); ypair = midpair_map(tup, clp; dims = 2)
    φxp = midpair_map(tup, φs; dims = 1);    φyp = midpair_map(tup, φs; dims = 2)
    V = Matrix{SVector{3, promote_type(eltype(μblock), Float64)}}(undef, nc, ns)
    for i in 1:nc, j in 1:ns
        tr = get_transformation(p[i, j])
        nbx1, nbx2 = tr.(xpair[i, j]); nby1, nby2 = tr.(ypair[i, j])
        φx1, φx2 = φxp[i, j]; φy1, φy2 = φyp[i, j]
        # local-frame gradient (∂φ/∂x, ∂φ/∂y) from both neighbour offsets (robust to sweep)
        ax = nbx1[1] - nbx2[1]; ay = nbx1[2] - nbx2[2]
        bx = nby1[1] - nby2[1]; by = nby1[2] - nby2[2]
        det = ax * by - ay * bx
        gx = ((φx1 - φx2) * by - (φy1 - φy2) * ay) / det
        gy = (ax * (φy1 - φy2) - bx * (φx1 - φx2)) / det
        # The doublet gradient is genuinely singular at a sharp trailing edge (the Kutta
        # kink) and the finite-difference stencil degenerates at a few thin/tapered edge
        # panels; cap the reconstructed speed so those isolated panels do not produce
        # absurd pressures/velocities. Forces should be taken from `lift_coefficients`.
        vloc = SVector(-gx, -gy, zero(gx))
        sp   = norm(vloc)
        vloc = sp > 3Vmag ? vloc * (3Vmag / sp) : vloc
        # local (ŝ, l̂, n̂) components → global: local_coordinate_system(p) has those as columns.
        V[i, j] = local_coordinate_system(p[i, j]) * vloc
    end
    V
end

_panel_pressures(p, μblock, Vmag) = map(v -> pressure_coefficient(Vmag, v), _panel_velocities(p, μblock, Vmag))

"""
    surface_velocities(system :: DoubletSourcePanelSystem)

Inviscid tangential edge velocity on every panel of every surface, as GLOBAL 3-vectors,
returned as a `NamedTuple` of matrices matching `system.surfaces`. Each vector is the
tangential gradient of the doublet (total-potential) distribution; its normal component is
zero by flow tangency, so these are directly usable as a wetted-surface boundary-layer EIF.
"""
function surface_velocities(system :: DoubletSourcePanelSystem)
    Vmag = norm(velocity(system.freestream))
    map(keys(system.surfaces)) do s
        _panel_velocities(system.surfaces[s], system.doublets[s], Vmag)
    end |> NamedTuple{keys(system.surfaces)}
end

"""
    lift_coefficients(system :: DoubletSourcePanelSystem)

Lift coefficient contribution of each lifting surface, referenced to `system.reference.area`,
computed by Kutta–Joukowsky from the shed wake circulation ``Γ = μ_{wake}``:
``C_L = -\\tfrac{2}{V_\\infty S} \\sum_j Γ_j Δy_j``. Robust for thin, swept surfaces where the
finite-difference surface speeds are noisy. Returns a `ComponentVector` keyed by component.
"""
function lift_coefficients(system :: DoubletSourcePanelSystem)
    V = norm(velocity(system.freestream)); S = system.reference.area
    cls = map(keys(system.surfaces)) do s
        Γ  = system.wake_strengths[s]
        Δy = [ abs(p4(w)[2] - p1(w)[2]) for w in system.wakes[s] ]
        -2 * sum(Γ .* Δy) / (V * S)
    end
    ComponentArray(NamedTuple{keys(system.surfaces)}(cls))
end

"""
    lift_coefficient(system :: DoubletSourcePanelSystem)

Total lift coefficient of the aircraft (sum of the lifting-surface contributions).
"""
lift_coefficient(system :: DoubletSourcePanelSystem) = sum(lift_coefficients(system))

"""
    field_velocity(system :: DoubletSourcePanelSystem, r)

Total induced velocity at a point `r` — the freestream plus every doublet panel (as a vortex
ring), every shed wake panel, and the fuselage source and cross-flow doublet lines. Use it to
trace streamlines through the solved field.
"""
function field_velocity(system :: DoubletSourcePanelSystem, r)
    V∞ = velocity(system.freestream)
    v  = V∞
    for s in keys(system.surfaces)
        pv = permutedims(system.surfaces[s])[:]
        μ  = system.doublets[s]
        for k in eachindex(pv); v += μ[k] * panel_ring_velocity(pv[k], r); end
        for (w, Γ) in zip(system.wakes[s], system.wake_strengths[s])
            v += Γ * panel_ring_velocity(w, r)
        end
    end
    if !isnothing(system.fuselage)
        for (el, λ) in zip(system.fuselage, system.fuse_doublets)
            v += point_doublet_velocity(λ * segment_length(el), el.normal, r - el.rc, el.radius)
            v += point_source_velocity(el.sigma * norm(V∞), r - el.rc, el.radius)
        end
    end
    v
end

# ---- forces, moments and coefficients ---------------------------------------------------
# Force and moment (per unit dynamic pressure, geometry axes) of one surface from its panel
# pressures. `normal_vector` points inward for `surface_panels` meshes, so the outward
# normal is its negation and the panel force is -cp·A·n̂_out.
function _surface_load(p, cp, r_ref)
    T = promote_type(eltype(cp), Float64)
    F = zero(SVector{3,T}); M = zero(SVector{3,T})
    for i in eachindex(p)
        dF = -cp[i] * panel_area(p[i]) * (-normalize(normal_vector(p[i])))
        F += dF
        M += (collocation_point(p[i]) - r_ref) × dF
    end
    F, M
end

# Slender-body (Munk) load of the fuselage line: sectional force ∝ dλ/dx along ẑ. Zero net
# lift over a closed body, but a nonzero pitching moment.
function _fuselage_load(fuse, λ, r_ref)
    T = promote_type(eltype(λ), Float64)
    F = zero(SVector{3,T}); M = zero(SVector{3,T})
    xs = [ el.rc[1] for el in fuse ]; n = length(fuse)
    for i in 1:n
        dλdx = i == 1 ? (λ[2]-λ[1])/(xs[2]-xs[1]) :
               i == n ? (λ[n]-λ[n-1])/(xs[n]-xs[n-1]) :
                        (λ[i+1]-λ[i-1])/(xs[i+1]-xs[i-1])
        dF = -dλdx * segment_length(fuse[i]) * fuse[i].normal
        F += dF
        M += (fuse[i].rc - r_ref) × dF
    end
    F, M
end

# Per-component force/moment coefficients (NF_COEFFS: CX,CY,CZ,Cl,Cm,Cn) in `axes`.
function nearfield_coefficients(system :: DoubletSourcePanelSystem; axes :: AbstractAxisSystem = Wind())
    refs = system.reference; S = refs.area; b = refs.span; c = refs.chord
    α, β = system.freestream.alpha, system.freestream.beta
    cps  = surface_coefficients(system)

    loads = map(s -> _surface_load(system.surfaces[s], cps[s], refs.location), keys(system.surfaces))
    ks    = collect(keys(system.surfaces))
    if !isnothing(system.fuselage)
        push!(ks, :fuse)
        loads = (loads..., _fuselage_load(system.fuselage, system.fuse_doublets, refs.location))
    end

    coeffs = map(loads) do (F, M)
        CF = _vector_to_axes(F, axes, α, β) / S
        Ma = _moment_to_axes(M, axes, α, β)
        CM = SVector(Ma[1]/(S*b), Ma[2]/(S*c), Ma[3]/(S*b))
        NF_COEFFS(CF..., CM...)
    end
    NamedTuple{Tuple(ks)}(coeffs)
end

"""
    nearfield(system :: DoubletSourcePanelSystem; axes = Wind())

Total nearfield force and moment coefficients `(CX, CY, CZ, Cl, Cm, Cn)` of the aircraft,
integrated from the surface pressures plus the fuselage slender-body (Munk) load. Reported
in wind axes by default.
"""
nearfield(system :: DoubletSourcePanelSystem; axes :: AbstractAxisSystem = Wind()) =
    NF_COEFFS(sum(nearfield_coefficients(system; axes)))

# 2-D vortex induced velocity in the Trefftz plane (streamwise x̂): x̂ × (r_i − r_j) / (2π|r|²).
_trefftz_velocity(r_i, r_j) = (r = r_i - r_j; SVector(1.0, 0, 0) × r / (2π * dot(r, r)))

"""
    farfield(system :: DoubletSourcePanelSystem)

Farfield force coefficients `(CDi, CY, CL)` from a Trefftz-plane integration of the shed wake
circulation, in wind axes. Following the vortex-lattice farfield, each lifting surface's wake
is projected into the plane normal to the freestream; the induced drag is
``C_{D_i} = -\\tfrac{1}{V_\\infty^2 S}\\sum Γ_j Δs_j (∂φ/∂n)_j`` with ``∂φ/∂n`` the downwash
from the trailing vorticity.
"""
function _surface_trefftz(system, s)
    fs = system.freestream; α, β = fs.alpha, fs.beta
    V = norm(velocity(fs)); S = system.reference.area; x̂ = SVector(1.0, 0, 0)
    wk = system.wakes[s]; Γ = system.wake_strengths[s]; ns = length(wk)
    T  = promote_type(eltype(Γ), Float64)
    ns == 0 && return FF_COEFFS(zero(T), zero(T), zero(T))
    nodes = [ [ p1(wk[1]) ]; [ p4(w) for w in wk ] ]       # ns+1 spanwise trailing-edge nodes
    nw    = geometry_to_wind_axes.(nodes, α, β)             # into the Trefftz (wind) frame
    proj  = [ (v = nw[i+1] - nw[i]; v - dot(x̂, v) * x̂) for i in 1:ns ]   # project out streamwise
    Δs    = norm.(proj)
    θ     = [ atan(v[3], v[2]) for v in proj ]             # dihedral angle in the plane
    n̂     = [ normalize(x̂ × v) for v in proj ]
    ctr   = [ (nw[i] + nw[i+1]) / 2 for i in 1:ns ]
    AICf  = [ dot(_trefftz_velocity(ctr[i], nw[k+1]) - _trefftz_velocity(ctr[i], nw[k]), n̂[i])
              for i in 1:ns, k in 1:ns ]
    ∂φ_∂n = AICf * Γ
    # Bound circulation is -μ_wake in this convention, hence the sign on CL, CY (consistent
    # with `lift_coefficients`). CDi is sign-independent (quadratic in the wake strength).
    CDi = -sum(Γ .* Δs .* ∂φ_∂n) / (V^2 * S)
    CY  =  2 * sum(Γ .* Δs .* sin.(θ)) / (V * S)
    CL  = -2 * sum(Γ .* Δs .* cos.(θ)) / (V * S)
    FF_COEFFS(CDi, CY, CL)
end

farfield_coefficients(system :: DoubletSourcePanelSystem) =
    NamedTuple{keys(system.surfaces)}(map(s -> _surface_trefftz(system, s), keys(system.surfaces)))

farfield(system :: DoubletSourcePanelSystem) =
    FF_COEFFS(reduce(+, values(farfield_coefficients(system))))

# ---- pretty-printing (reuses the vortex-lattice tables) ---------------------------------
"""
    print_coefficients(system :: DoubletSourcePanelSystem, name = :aircraft; components = false)

Print a pretty table of the nearfield and farfield coefficients, matching the vortex-lattice
output. Pass `components = true` to also print a table per component.
"""
function print_coefficients(system :: DoubletSourcePanelSystem, name = :aircraft; components = false)
    if components
        nfs = nearfield_coefficients(system; axes = Wind())
        ffs = farfield_coefficients(system)
        for key in keys(nfs)
            print_coefficients(nfs[key], haskey(ffs, key) ? ffs[key] : FF_COEFFS(0.0, 0.0, 0.0), key)
        end
    end
    print_coefficients(nearfield(system), farfield(system), name)
    nothing
end

# ---- stability derivatives --------------------------------------------------------------
# Nine coefficients (CX,CY,CZ,Cl,Cm,Cn,CDi,CYff,CL) of every component, plus the total.
function _stability_coefficients(system, axes)
    nfs = nearfield_coefficients(system; axes)
    ffs = farfield_coefficients(system)          # surfaces only (no wake ⇒ no fuse farfield)
    ks  = keys(nfs)
    nine = map(ks) do k
        nf = nfs[k]
        ff = haskey(ffs, k) ? ffs[k] : FF_COEFFS(0.0, 0.0, 0.0)
        SVector(nf.CX, nf.CY, nf.CZ, nf.Cl, nf.Cm, nf.Cn, ff.CDi, ff.CY, ff.CL)
    end |> NamedTuple{ks}
    nine, reduce(+, values(nine))
end

# Reconstruct the aircraft NamedTuple and the wake length used to build the stored system.
_aircraft(system) = (; system.surfaces...,
                     (isnothing(system.fuselage) ? (;) : (fuse = system.fuselage,))...)
_wake_length(system) = (w = first(first(values(system.wakes))); norm(p2(w) - p1(w)))

"""
    freestream_derivatives(
        system :: DoubletSourcePanelSystem;
        axes = Stability(), name = :aircraft,
        print = false, print_components = false, farfield = false,
    )

Force/moment coefficients and their derivatives with respect to angle of attack ``α`` and
sideslip ``β``, obtained by automatic differentiation (`ForwardDiff`) through the coupled
solve. Reported in stability axes by default and returned as a `NamedTuple` of `Derivs` per
component and for the whole aircraft.

NOTE: this potential (Dirichlet) formulation has no rotational-onset model, so the Mach and
quasi-steady rate derivatives (``∂/∂M``, ``∂/∂p̄``, ``∂/∂q̄``, ``∂/∂r̄``) are reported as zero.
"""
function freestream_derivatives(system :: DoubletSourcePanelSystem;
        axes :: AbstractAxisSystem = Stability(), name = :aircraft,
        print = false, print_components = false, farfield = false)

    ac  = _aircraft(system); refs = system.reference; wl = _wake_length(system)
    ks  = (keys(system.surfaces)..., (isnothing(system.fuselage) ? () : (:fuse,))...)

    # Flat coefficient vector over x = [α, β]: [ comp₁(9); comp₂(9); …; total(9) ].
    function coeff_vector(x)
        sysd = solve_case(ac, Freestream(alpha = rad2deg(x[1]), beta = rad2deg(x[2])), refs; wake_length = wl)
        nt, tot = _stability_coefficients(sysd, axes)
        reduce(vcat, (values(nt)..., tot))
    end

    x0 = [ system.freestream.alpha, system.freestream.beta ]
    y  = coeff_vector(x0)
    J  = ForwardDiff.jacobian(coeff_vector, x0)   # rows: 9·(#components + 1), cols: [∂/∂α ∂/∂β]

    z9 = zeros(eltype(y), 9)
    allkeys = (ks..., name)
    comps = ntuple(length(allkeys)) do b
        rows = 9*(b-1)+1 : 9*(b-1)+9
        Derivs(hcat(y[rows], z9, J[rows, 1], J[rows, 2], z9, z9, z9))
    end |> NamedTuple{allkeys}

    if print_components
        for k in allkeys; print_derivatives(comps[k], k; farfield); end
    elseif print
        print_derivatives(comps[name], name; farfield)
    end

    return comps
end
