## Wetted-area method
#======================================================#

# Raymer form factor for wing, pylon, nacelle, fuselage? (Eq. 12.30)
form_factor_wing(x_c, t_c, Λ_m, M) = (1 + 0.6t_c / x_c + 100t_c^4) * (1.34M^0.18 * cos(Λ_m)^0.28)

function form_factor(wing :: HalfWing, M)
    # (x/c)_max, (t/c)_max, (t/c)_max sweep angles
    xcs, tcs, sweeps = max_thickness_to_chord_ratio_sweeps(wing, 60)
    Kf  = form_factor_wing.(xcs, tcs, sweeps, M)
end

# Schlichting averaged skin-friction coefficients
cf_lam(Re_c, k_lam = 1.) = 1.328 / √(Re_c * k_lam)
cf_turb(Re_c, M) = 0.455 / log10(Re_c)^2.58 / (1 + 0.144M^2)^0.65
cf_schlichting(Re, Re_xtr, M) = max(cf_lam(Re), cf_turb(Re, M) - (Re_xtr / 320 - 39) / Re)

function wetted_area_drag(mean_chords, S_wets, K_fs, x_tr, V, ρ, M, μ)
    # Skin-friction coefficients
    Re_c    = reynolds_number.(ρ, V, mean_chords, μ)
    Re_xtr  = reynolds_number.(ρ, V, mean_chords .* x_tr, μ)
    cfs     = cf_schlichting.(Re_c, Re_xtr, M)

    # Profile drag
    Dp_by_q = sum(@. cfs * S_wets * K_fs)
end

function wetted_area_drag(wing :: HalfWing, x_tr, V, ρ, a_ref = 330., μ = 1.5e-5)
    # Chord processing
    mean_chords = (fwdsum ∘ chords)(wing) / 2

    # Wetted areas
    S_wets  = @. mean_chords * wing.spans / cos(wing.dihedrals)

    # Form factors
    M       = V / a_ref
    K_fs    = form_factor(wing, M)

    wetted_area_drag(mean_chords, S_wets, K_fs, x_tr, V, ρ, M, μ)
end

# Sato's local-friction and local-dissipation based on power balance method from Mark Drela, Flight Vehicle Aerodynamics, eq. 4.115.
function local_dissipation_drag(wing :: Wing, wetted_areas, ρ_es, u_es, x_tr, V, ρ, M, μ)
    # Chord processing
    mean_chords = (fwdsum ∘ chords)(wing) / 2

    # Compute weighted wetted areas based on inviscid edge velocity distribution.
    weighted_S_wets = sum(x -> x[1] * norm(x[2])^3 * x[3], zip(ρ_es, u_es, wetted_areas), dims = 1) ./ (ρ * V^3)

    wetted_area_drag(mean_chords, weighted_S_wets, 1., x_tr, V, ρ, M, μ)
end

profile_drag_coefficient(wing :: HalfWing, x_tr, V, rho_ref, a_ref, area_ref, μ) = wetted_area_drag(wing, x_tr, V, rho_ref, a_ref, μ) / area_ref
profile_drag_coefficient(wing :: Wing, x_tr, V, rho_ref, a_ref, area_ref, μ) = profile_drag_coefficient(wing.left, x_tr, V, rho_ref, a_ref, area_ref, μ) + profile_drag_coefficient(wing.right, x_tr, V, rho_ref, a_ref, area_ref, μ)