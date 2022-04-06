## Wetted-area method
#======================================================#

# Raymer form factor for wing, pylon, nacelle, fuselage? (Eq. 12.30)
form_factor_wing(x_c, t_c, Λ_m, M) = (1 + 0.6t_c / x_c + 100t_c^4) * (1.34M^0.18 * cos(Λ_m)^0.28)

function form_factor(wing :: HalfWing, M)
    # (x/c)_max, (t/c)_max, (t/c)_max sweep angles
    xcs, tcs, sweeps = max_thickness_to_chord_ratio_sweeps(wing, 60)
    Kf  = form_factor_wing.(xcs, tcs, sweeps, M)
end

form_factor(wing :: Wing, M) = (form_factor(wing.left, M) + form_factor(wing.right, M)) / 2

# Schlichting averaged skin-friction coefficients
cf_lam(Re_c, k_lam = 1.) = 1.328 / √(Re_c * k_lam)
cf_turb(Re_c, M) = 0.455 / log10(Re_c)^2.58 / (1 + 0.144M^2)^0.65
cf_schlichting(Re, Re_xtr, M) = max(cf_lam(Re), cf_turb(Re, M) - (Re_xtr / 320 - 39) / Re)

function skin_friction_coefficients(mean_chords, x_tr, V, ρ, M, μ)
    # Reynolds numbers for chord lengths
    Re_c = reynolds_number.(ρ, V, mean_chords, μ)

    # Reynolds numbers at transition locations
    Re_xtr = reynolds_number.(ρ, V, mean_chords .* x_tr, μ)

    # Skin-friction coefficients
    cfs = cf_schlichting.(Re_c, Re_xtr, M)
end

function wetted_area_drag(mean_chords, S_wets, K_fs, x_tr, V, ρ, M, μ)
    cfs = skin_friction_coefficients(mean_chords, x_tr, V, ρ, M, μ)
    
    # Profile drag
    Dp_by_q = sum(@. cfs * S_wets * K_fs)
end

function wetted_area_drag(wing :: HalfWing, x_tr, V, ρ, a_ref = 330., μ = 1.5e-5)
    # Chord processing
    mean_chords = (forward_sum ∘ chords)(wing) / 2

    # Wetted areas
    S_wets  = @. mean_chords * wing.spans / cos(wing.dihedrals)

    # Form factors
    M       = V / a_ref
    K_fs    = form_factor(wing, M)

    wetted_area_drag(mean_chords, S_wets, K_fs, x_tr, V, ρ, M, μ)
end

function wetted_area_drag(wing :: WingMesh, x_tr, V, ρ, a_ref = 330., μ = 1.5e-5)
    # Chord processing
    mean_chords = (forward_sum ∘ chords)(wing.surface) / 2

    # Wetted areas
    surf_pans = camber_panels(wing)
    S_wets    = sum(panel_area, surf_pans, dims = 1)

    # Form factors
    M       = V / a_ref
    K_fs    = form_factor(wing.surface, M)

    wetted_area_drag(mean_chords, S_wets, K_fs, x_tr, V, ρ, M, μ)
end

# Sato's local-friction and local-dissipation based on power balance method from Mark Drela, Flight Vehicle Aerodynamics, eq. 4.115.
function local_dissipation_drag(wing :: Wing, wetted_areas, ρ_es, u_es, x_tr, V, ρ, M, μ)
    # Chord processing
    mean_chords = (forward_sum ∘ chords)(wing) / 2

    # Compute weighted wetted areas based on inviscid edge velocity distribution.
    weighted_S_wets = sum(@. ρ_es * norm(u_es)^3 * wetted_areas; dims = 1) ./ (ρ * V^3)

    wetted_area_drag(mean_chords, weighted_S_wets, 1., x_tr, V, ρ, M, μ)
end

function local_dissipation_drag(wing :: WingMesh, ρ_es, u_es, x_tr, V, ρ, M, μ)
    # Chord processing
    mean_chords = (forward_sum ∘ chords)(wing.surface) / 2

    # Compute weighted wetted areas based on inviscid edge velocity distribution.
    S_wets = panel_area.(camber_panels(wing))
    weighted_S_wets = sum(@. ρ_es * u_es^3 * S_wets; dims = 1) ./ (ρ * V^3)

    wetted_area_drag(mean_chords, weighted_S_wets, 1., x_tr, V, ρ, M, μ)
end

abstract type AbstractProfileDrag end

struct FormFactor <: AbstractProfileDrag end
struct Dissipation <: AbstractProfileDrag end

profile_drag_coefficient(wing :: HalfWing, x_tr, V, rho_ref, a_ref, area_ref, μ) = wetted_area_drag(wing, x_tr, V, rho_ref, a_ref, μ) / area_ref
profile_drag_coefficient(wing :: WingMesh, x_tr, V, rho_ref, a_ref, area_ref, μ) = wetted_area_drag(wing, x_tr, V, rho_ref, a_ref, μ) / area_ref
profile_drag_coefficient(wing :: Wing, x_tr, V, rho_ref, a_ref, area_ref, μ) = profile_drag_coefficient(wing.left, x_tr, V, rho_ref, a_ref, area_ref, μ) + profile_drag_coefficient(wing.right, x_tr, V, rho_ref, a_ref, area_ref, μ)

profile_drag_coefficient(wing :: AbstractWing, x_tr, refs :: References) = profile_drag_coefficient(wing, x_tr, refs.speed, refs.density, refs.sound_speed, refs.area, refs.viscosity)

profile_drag_coefficient(wing :: AbstractWing, x_tr, edge_speeds, panels, refs :: References) = local_dissipation_drag(wing, panel_area.(panels), refs.density, edge_speeds, x_tr, refs.speed, refs.density, mach_number(refs), refs.viscosity) / refs.area

profile_drag_coefficient(wing :: WingMesh, x_tr, edge_speeds, refs :: References) = local_dissipation_drag(wing, refs.density, edge_speeds, x_tr, refs.speed, refs.density, mach_number(refs), refs.viscosity) / refs.area

function wave_drag(M, Λ, t_by_c, Cl, κ_A)
    M_drag_divergence = κ_A / cos(Λ) - t_by_c / cos(Λ)^2 - Cl / (10 * cos(Λ)^3) # Drag divergence Mach number
    M_crit   = M_drag_divergnece - (0.1 / 80)^(1/3) # Critical Mach number
    CD_wave = ifelse(M > M_crit, 20 * (Mach - M_crit)^4, 0.)
end