## Wetted-area method
#==========================================================#

# Schlichting averaged skin-friction coefficients
cf_lam(Re_c, k_lam = 1.) = 1.328 / √(Re_c * k_lam)
cf_turb(Re_c, M) = 0.455 / log10(Re_c)^2.58 / (1 + 0.144M^2)^0.65
cf_schlichting(Re, Re_xtr, M) = max(cf_lam(Re), cf_turb(Re, M) - (Re_xtr / 320 - 39) / Re)

"""
    parasitic_drag_coefficient(
        L, x_tr, 
        ρ, V, M, μ, S_ref, 
        S_wet, Kf, fM
    )

Compute the parasitic drag coefficient ``C_{D₀}`` with the following quantities. Uses a transition-based model from laminar to turbulent flow based on Schlichting's skin-friction coefficient formulas.

# Arguments
- L     = Reference length (m)
- x_tr  = Transition location as ratio of reference length
- ρ     = Density (kg/m³)
- V     = Speed (m/s)
- M     = Mach number
- S_ref = Reference area (m²)
- μ     = Dynamic viscosity (kg/(m ⋅ s))
- S_wet = Wetted area (m²)
- Kf    = Form factor
- fM    = Mach number influence
"""
function parasitic_drag_coefficient(L, x_tr, ρ, V, M, μ, S_ref, S_wet, Kf, fM)
    Re_L = reynolds_number(ρ, V, L, μ) # Reynolds number at reference length
    Re_xtr = reynolds_number(ρ, V, L * x_tr, μ) # Reynolds number at transition locations (ratios of chords)
    cf = cf_schlichting(Re_L, Re_xtr, M) # Skin-friction coefficient
    
    return cf * Kf * fM * S_wet / S_ref
end

## Fuselage
#==========================================================#

# Raymer form factor for fuselage and nacelle, 4th edition (Eq. 12.31)
form_factor_fuselage(f) = 1 + 60 / f^3 + f / 400

"""
    form_factor(fuse :: HyperEllipseFuselage)

Compute the form factor of `HyperEllipseFuselage` for parasitic drag estimation based on the DATCOM formula (Raymer - Aircraft Design: A Conceptual Approach, 4th Ed., Eq. 12.31): ``FF = 1 + 60 / f³ + f / 400``
"""
form_factor(fuse :: HyperEllipseFuselage) = form_factor_fuselage(fuse.length / 2fuse.radius)

function wetted_area_drag_coefficient(fuse :: HyperEllipseFuselage, x_tr, ρ, V, M, μ, S_ref, ts = 0:0.01:1)
    # Fuselage quantities
    L = fuse.length
    S_wet = wetted_area(fuse, ts) 
    Kf = form_factor(fuse)
    fM = (1 - 0.08M^1.45)

    # Parasitic drag coefficient
    CDp = parasitic_drag_coefficient(L, x_tr, ρ, V, M, μ, S_ref, S_wet, Kf, fM)

    return CDp
end

## Wing
#==========================================================#

# Computing maximum thickness-to-chord ratios, their locations, and sweep angles
@views function mean_max_thickness_chord_sweeps(wing :: Wing, num)
    # Compute (t/c)_max locations and values
    xs_max_tbyc = maximum_thickness_to_chord(wing, num) # Compute ((x/c), (t/c)_max)
    max_tbyc = getindex.(xs_max_tbyc, 2) # Get maximum (t/c) for each chord
    xs_temp = getindex.(xs_max_tbyc, 1) # Get (x/c) at max. (t/c) for each chord
    
    # Compute leading-edge sweep angles accounting for dihedral
    le = leading_edge(wing) # Get leading edge
    xs = @. le[:,1] + wing.chords * xs_temp # Determine max. thickness x-coordinates in geometry frame
    ds = xs[2:end] - xs[1:end-1] # Compute differences between max. thickness x-coordinates between sections
    widths = @. wing.spans / cosd(wing.dihedrals) # Project planform span onto geometry based on dihedral angles
    Λs = @. atan(ds, widths) # Now compute sweep angles at max. thickness locations

    # Average chord lengths for sections
    xs = (xs_temp[1:end-1] + xs_temp[2:end]) / 2
    tbycs = (max_tbyc[1:end-1] + max_tbyc[2:end]) / 2

    xs, tbycs, Λs
end

# Raymer form factor for wing, tail, strut, and pylon, 4th edition (Eq. 12.30). "_m" refers to corresonding values at (t/c)_max here.
form_factor_wing(x_c_m, t_c_max, Λ_m) = (1 + 0.6t_c_max / x_c_m + 100t_c_max^4) * cos(Λ_m)^0.28

"""
    form_factor(wing :: Wing, M, num)
    form_factor(wing_mesh :: WingMesh, M, num)

Compute the form factor ``K_f`` of a `Wing` or `WingMesh` at Mach number ``M`` for parasitic drag estimation based on Raymer's formula for wing, tail, strut, and pylon (Raymer - Aircraft Design: A Conceptual Approach, 4th Ed., Eq. 12.30): ``FF = (1 + 0.6(t/c)ₘₐₓ / (x/c) + 100(t/c)ₘₐₓ^4) ⨯ 1.34M^0.18 × cos(Λ_m)^0.28``.

An integer `num` is required for interpolating `Wing`'s airfoils' coordinates to determine their maximum thickness-to-chord ratios.
"""
function form_factor(wing :: Wing, num :: Integer)
    # (x/c)_max, (t/c)_max, (t/c)_max sweep angles
    xcs, tcs, sweeps = mean_max_thickness_chord_sweeps(wing, num)
    Kf = @. form_factor_wing(xcs, tcs, sweeps)

    return Kf
end

# Wetted area calculation based on planform area, form factors Kf, and Mach number correction fM.
@views function wetted_area_drag_coefficient(wing :: Wing, x_tr, ρ, V, M, μ, S_ref, num)
    # Averaging chords over spans
    mean_chords = (wing.chords[1:end-1] + wing.chords[2:end]) / 2

    # Wetted areas for integration over averaged chords
    S_wets = @. mean_chords * wing.spans / cos(wing.dihedrals)

    K_fs = form_factor(wing, num) # Form factors
    fM = 1.34M^0.18  # Mach number correction

    # Calculate wetted area drag coefficient with wetted areas and form factors
    CDp = sum(zip(mean_chords, S_wets, K_fs)) do (c, S_wet, K_f)
        parasitic_drag_coefficient(c, x_tr, ρ, V, M, μ, S_ref, S_wet, K_f, fM)
    end

    # Double if symmetric
    CDp = ifelse(wing.symmetry, 2CDp, CDp)

    return CDp
end

# Sato's local-friction and local-dissipation based on power balance method from Mark Drela, Flight Vehicle Aerodynamics, eq. 4.115. Additional reference: Sato's thesis (https://dspace.mit.edu/handle/1721.1/75837), eq. 4.10-11
function local_dissipation_drag_coefficient(wing :: Wing, S_wets, ρ_es, u_es, x_tr, ρ, V, M, μ, S_ref)
    # Averaging chords over spans
    mean_chords = @views (wing.chords[1:end-1] + wing.chords[2:end]) / 2

    # Compute weighted wetted areas based on edge velocity distribution for dissipation.
    S_wets_diss = sum(@. ρ_es * norm(u_es)^3 * S_wets; dims = 1) ./ (ρ * V^3)

    # Now calculate wetted area drag with the weighted wetted areas
    Kf, fM = 1., 1. # Unit form factors and Mach number correction
    CDp = sum(@. parasitic_drag_coefficient(mean_chords, x_tr, ρ, V, M, μ, S_ref, S_wets_diss, Kf, fM))

    # Double if symmetric
    CDp = ifelse(wing.symmetry, 2CDp, CDp) 

    return CDp
end

## Interfaces for public function with References
#==========================================================#

"""
    parasitic_drag_coefficient(
        fuse :: HyperEllipseFuselage, 
        refs :: References, 
        x_tr :: Real, 
        ts = 0:0.01:1
    )

Estimate the profile drag coefficient of a `HyperEllipseFuselage` using the **wetted-area method** based on Schlichting's skin-friction coefficient formula with given `References`, a specified transition location ``xₜᵣ`` as a ratio of the fuselage length, and optionally the parametric distribution `t ∈ [0,1]` for the discretization of the nose, cabin and rear sections.
"""
parasitic_drag_coefficient(fuse :: HyperEllipseFuselage, refs :: References, x_tr :: Real, ts = 0:0.01:1) = wetted_area_drag_coefficient(fuse, x_tr, refs.density, refs.speed, mach_number(refs), refs.viscosity, refs.area, ts)

"""
    parasitic_drag_coefficient(
        wing :: Wing, 
        refs :: References,
        x_tr
    )

Estimate the profile drag coefficient of a `Wing` using the **wetted-area method** based on Schlichting's skin-friction coefficient formula with given `References`, a specified transition location ``xₜᵣ ∈ [0,1]`` as a ratio of the chord lengths.
"""
parasitic_drag_coefficient(wing :: Wing, refs :: References, x_tr :: Real; num = 60) = wetted_area_drag_coefficient(wing, x_tr, refs.speed, refs.density, mach_number(refs), refs.viscosity, refs.area, num)

# Forwarding common methods for WingMesh
MacroTools.@forward WingMesh.surface parasitic_drag_coefficient, form_factor

"""
    parasitic_drag_coefficient(
        wing :: WingMesh, 
        refs :: References
        x_tr :: Real, 
        u_es, 
    )

Estimate the profile drag coefficient of a `WingMesh` using the **local-friction and local-dissipation method** based on Schlichting's skin-friction coefficient formula with given `References`, a specified transition ratio ``xₜᵣ``, and edge velocities ``\\mathbf u_e``. 

At present, the edge velocities would be computed using the vortex lattice method via `VortexLatticeSystem`. For this case, the panels corresponding to the camber distribution are used in the calculation. 
"""
parasitic_drag_coefficient(wing :: WingMesh, refs :: References, x_tr :: Real, u_es) = local_dissipation_drag_coefficient(wing.surface, map(panel_area, camber_panels(wing)), refs.density, u_es, x_tr, refs.density, refs.speed, mach_number(refs), refs.viscosity, refs.area)
# Should it be doubled for upper and lower surfaces?

# FOR FUTURE USE WITH OTHER INVISCID AERODYNAMICS METHODS
#==========================================#

# """
#     parasitic_drag_coefficient(
#         wing :: Wing, 
#         refs :: References
#         x_tr :: Real,
#         u_es, 
#         panels :: Array{Panel3D}, 
#     )

# Estimate the profile drag coefficient of a `Wing` using the **local-friction and local-dissipation method** based on Schlichting's skin-friction coefficient formula with given `References`, a specified transition location ``xₜᵣ ∈ [0,1]`` as a ratio of the chord lengths, and edge velocities ``\\mathbf u_e`` with corresponding panels.
# """
# parasitic_drag_coefficient(wing :: Wing, refs :: References, x_tr :: Real, u_es, panels) = local_dissipation_drag_coefficient(wing, panel_area.(panels), refs.density, u_es, x_tr, refs.density, refs.speed, mach_number(refs), refs.viscosity, refs.area)

## Wave drag
#==========================================================#

function wave_drag(M, Λ, t_by_c, Cl, κ_A)
    M_drag_divergence = κ_A / cosd(Λ) - t_by_c / cosd(Λ)^2 - Cl / (10 * cosd(Λ)^3) # Drag divergence Mach number
    M_crit  = M_drag_divergence - (0.1 / 80)^(1/3) # Critical Mach number
    CDw = ifelse(M > M_crit, 20 * (Mach - M_crit)^4, 0.)

    return CDw
end

function wave_drag_lift(M, S, L, b, CL)
    @assert M > 1 "This only applies to supersonic flow (M > 1)!"

    K_wL = 2 * (S / (b * L))^2

    CDw = K_wL * CL^2 * (M^2 - 1) / (2π * L^2)

    return CDw
end