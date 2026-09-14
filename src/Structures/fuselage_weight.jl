## Fuselage structural weight — reduced TASOPT-style pressurized shell + bending beam
#==========================================================================================#

"""
    fuselage_structural_mass(
        R, L_shell, S_shell, Δp, material :: Material;
        n_press = 2.0, t_min = 1.2e-3,
        P_tail = 0.0, l_tail = L_shell / 2,
        f_frame = 0.25, w_floor = 16.0, n_bulkhead = 2,
    )

Primary structural mass (kg) of a pressurized fuselage from a reduced TASOPT-style model
(Drela, *Simplified Aircraft Design*): a hoop-stress-sized cylindrical skin with a minimum-gauge
floor, frames/stringers as a fraction of the skin, top/bottom bending caps reacting an empennage
load ``P_{tail}`` at moment arm ``l_{tail}``, floor beams over the pressurized length, and fore/aft
pressure bulkheads. Non-structural items (furnishings, insulation, systems) are excluded.

The skin is sized by the cabin hoop stress ``σ = Δp\\,R / t`` (with a proof factor `n_press`),
floored at `t_min`. The bending caps form a couple at ``±R`` reacting the moment
``M = P_{tail}\\,l_{tail}``; the cap area ``A = M/(σR)`` tapers linearly to the tail, giving a
material volume ``A\\,l_{tail}/2``.

# Arguments
- `R`          : fuselage (effective) radius (m)
- `L_shell`    : pressurized shell length (m)
- `S_shell`    : pressurized shell surface (wetted) area (m²)
- `Δp`         : cabin–ambient pressure differential (Pa)
- `material`   : `Material` supplying `yield_stress` (Pa) and `density` (kg/m³)
- `n_press`    : pressure proof/burst factor of safety
- `t_min`      : minimum skin gauge (m), a damage-tolerance floor
- `P_tail`     : limit empennage (down/side) load bending the aft fuselage (N)
- `l_tail`     : moment arm from the wing box to the empennage load (m)
- `f_frame`    : frame + stringer mass as a fraction of the skin mass
- `w_floor`    : floor structural areal density (kg/m²)
- `n_bulkhead` : number of pressure bulkheads (domes)
"""
function fuselage_structural_mass(
        R, L_shell, S_shell, Δp, material :: Material;
        n_press = 2.0, t_min = 1.2e-3,
        P_tail = 0.0, l_tail = L_shell / 2,
        f_frame = 0.25, w_floor = 16.0, n_bulkhead = 2,
    )
    σ = material.yield_stress
    ρ = material.density

    # 1. Pressure shell — hoop stress σ = Δp R / t, floored at the minimum gauge
    t_skin  = max(n_press * Δp * R / σ, t_min)
    m_skin  = ρ * t_skin * S_shell

    # 2. Frames & stringers
    m_frame = f_frame * m_skin

    # 3. Bending caps reacting the tail load: couple at ±R carries M = P_tail·l_tail,
    #    root cap area A = M/(σR) tapering linearly to the tail ⇒ volume ≈ A·l_tail/2
    m_bend  = ρ * (P_tail * l_tail / (σ * R)) * l_tail / 2

    # 4. Floor beams over the pressurized cabin (floor width ≈ 2R)
    m_floor = w_floor * (2R * L_shell)

    # 5. Pressure bulkheads — hemispherical domes, hoop σ = Δp R / 2t, area 2πR²
    t_bulk  = max(n_press * Δp * R / 2σ, t_min)
    m_bulk  = ρ * n_bulkhead * t_bulk * (2π * R^2)

    return m_skin + m_frame + m_bend + m_floor + m_bulk
end

"""
    fuselage_structural_mass(fuse :: Fuselage, Δp, material :: Material; kwargs...)

Estimate the fuselage structural mass (kg) using the maximum radius, the length, and the full
wetted area of a `Fuselage`. Extra keyword arguments are forwarded (e.g. `P_tail`, `l_tail`).
"""
fuselage_structural_mass(fuse :: Fuselage, Δp, material :: Material; kwargs...) =
    fuselage_structural_mass(maximum(fuse.radii), length(fuse), wetted_area(fuse), Δp, material; kwargs...)

"""
    fuselage_structural_mass(fuse :: HyperEllipseFuselage, Δp, material :: Material; ts = 0:0.1:1, kwargs...)

Estimate the fuselage structural mass (kg) using the radius, the length, and the wetted area of a
`HyperEllipseFuselage`. Extra keyword arguments are forwarded (e.g. `P_tail`, `l_tail`).
"""
fuselage_structural_mass(fuse :: HyperEllipseFuselage, Δp, material :: Material; ts = 0:0.1:1, kwargs...) =
    fuselage_structural_mass(fuse.radius, fuse.length, wetted_area(fuse, ts), Δp, material; kwargs...)
