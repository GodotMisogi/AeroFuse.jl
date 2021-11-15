## Residual BL equations using finite differences
#============================================#

# mutable struct BoundaryLayer{T <: Real}
#   U           :: Real
#   Θ           :: Real
#   δ_star      :: Real
#   H_star      :: Real
#   H_star_star :: Real
#   CF          :: Real
#   CD          :: Real
#   C_τ_EQ      :: Real
#   U_s         :: Real
# end

# BoundaryLayer(U :: T, θ :: T, δ_star :: T, H_star :: T, H_star_star :: T, CF :: T, CD :: T, c_τ_EQ :: T = 0., U_s :: T = 0.) where T <: Real = BoundaryLayer{T}(U, θ, δ_star, H_star, H_star_star, CF, CD, c_τ_EQ
# , U_s)

# Momentum equation
momentum_finite_difference(dlogθ, dlogU, ΔX, θ, H, CF, M_e_sq) = dlogθ + (H + 2 - M_e_sq) * dlogU - CF * ΔX / 2θ # Check Mₑ² term later

# Shape parameter equation
shape_finite_difference(dlogH_star, dlogU, ΔX, θ, H, H_star, H_star_star, CF, CD) = dlogH_star + (2H_star_star / H_star + 1 - H) * dlogU + (CF / 2 - 2CD / H_star) * ΔX / θ

## Turbulence closure equations using finite differences
#============================================#

# Rate equation for maximum shear stress coefficient, c_τ (Green)
stress_transport_finite_difference(c_τ, dlogc_τ, dlogU, Δx, δ, δ_star, CF, c_τ_EQ, H_k) = dlogc_τ / Δx - 2 * (4. / 3δ_star * (CF / 2 - ((H_k - 1) / 6.7H_k)^2)) - dlogU / Δx - 5.6(√c_τ_EQ - √c_τ) / δ

# Amplification ratio equation for Tollmien-Schlichting waves using envelope eⁿ, ñ
amplification_finite_difference(Δn, Δx, H_k, θ) = Δn/Δx - dñdReθ(H_k) * (m(H_k) + 1) / 2 * l(H_k) / θ

## Helper functions
#============================================#

# Laminar quantity correlations
function laminar_body_quantities(H_k, M_e, Re_θ)
    H_star      = lam_H_star(H_k)
    H_star_star = turb_H_star_star(H_k, M_e)
    CF          = 2 * lam_CF(H_k) / Re_θ
    CD          = H_star * lam_CD(H_k)  / 2Re_θ

    H_k, H_star, H_star_star, CF, CD
end

# Turbulent quantity correlations
function turbulent_body_quantities(H_k, H, M_e, Re_θ)
    H_star      = turb_H_star(H_k, Re_θ)
    H_star_star = turb_H_star_star(H_k, M_e)
    CF          = turb_CF(H_k, Re_θ, M_e)
    U_s         = wall_slip(H_star, H_k, H)
    c_τ_EQ      = c_τ_equilibrium(H, H_k, H_star, U_s)
    CD          = turb_CD(CF, U_s, c_τ_EQ)

    H_star, H_star_star, CF, CD, c_τ_EQ, U_s
end

function check_transition(panel1, panel2, θ, p1, p2, δ_star_1, δ_star_2, H_1, H_2, H_k_1, H_k_2, M_e_1, M_e_2, Re_θ_1, Re_θ_2, n_crit)

    p_a   = (p1 + p2) / 2

    if p_a < n_crit
        # println("Laminar")
        f1 = laminar_body_quantities(H_k_1, M_e_1, Re_θ_1)
        f2 = laminar_body_quantities(H_k_2, M_e_2, Re_θ_2)

        dlogf   = @. log(f2 / f1)
        fa      = weighted_vector(f1, f2, 0.5)

        dlogH_star, dlogH_star_star, dlogCF, dlogCD = dlogf
        H_star, H_star_star, CF, CD = fa

        Δn = p2 - p1

        R3 = amplification_finite_difference(Δn, ΔX, H, θ)
    else
        # println("Turbulent")
        H     = (H_1 + H_2) / 2
        H_k   = (H_k_1 + H_k_2) / 2

        f1 = turbulent_body_quantities(H_k, H_1, M_e_1, Re_θ_1)
        f2 = turbulent_body_quantities(H_k, H_2, M_e_2, Re_θ_2)

        dlogf   = @. log(f2 / f1)[1:end-2]
        fa      = weighted_vector(f1, f2, 0.5)[1:end-1]

        dlogH_star, dlogH_star_star, dlogCF, dlogCD = dlogf
        H_star, H_star_star, CF, CD, c_τ_EQ         = fa

        dlogc_τ = log(p2 / p1)
        c_τ     = p_a # ???
        δ_star  = (δ_star_1 + δ_star_2) / 2
        δ       = δ_thick(H, θ, δ_star)

        R3 = stress_transport_finite_difference(c_τ, dlogc_τ, dlogU, ΔX, δ, δ_star, CF, c_τ_EQ, H)
    end

    return R3
end


function evaluate_residuals(panel1, panel2, U1, U2, θ1, θ2, p1, p2, x1, x2, δ_star_1, δ_star_2, H_1, H_2, H_k_1, H_k_2, M_e_1, M_e_2, Re_θ_1, Re_θ_2, n_crit)

    # Differences
    dlogθ = log(θ2 / θ1)
    dlogU = log(U2 / U1)

    ΔX     = x2 - x1
    M_e_sq = M_e_2^2 - M_e_1^2

    # Points
    θ     = (θ1 + θ2) / 2
    H_k   = (H_k_1 + H_k_2) / 2

    R1 = momentum_finite_difference(dlogθ, dlogU, ΔX, θ, H_k, CF, M_e_sq)
    R2 = shape_finite_difference(dlogH_star, dlogU, ΔX, θ, H_k, H_star, H_star_star, CF, CD)
    R3 = check_transition(panel1, panel2, θ, p1, p2, δ_star_1, δ_star_2, H_1, H_2, H_k_1, H_k_2, M_e_1, M_e_2, Re_θ_1, Re_θ_2, n_crit)

    R1, R2, R3
end

## Boundary layer system
#============================================#

function boundary_layer_finite_difference(Us, δ_stars, θs, Ts, ΔXs, panels, n_crit, a, ν)
    # Auxillary variables
    M_es    = Us ./ a
    Hs      = δ_stars ./ θs
    Re_θs   = reynolds_number.(Us, θs, ν)

    H_ks    = H_Whitfield.(Hs, M_es)

    # println(length.([ Us, θs, Ts, ΔXs, δ_stars, Hs, H_ks, M_es, Re_θs]))

    Rs = evaluate_residuals.(panels[1:end-1], panels[2:end],     # Foil panels and wake panels
                             Us[1:end-1], Us[2:end],             # Edge velocities
                             θs[1:end-1], θs[2:end],             # Momentum thicknesses
                             Ts[1:end-1], Ts[2:end],             # Amplification ratio or shear stress coefficeint
                             ΔXs[1:end-1], ΔXs[2:end],           # Panel arc lengths
                             δ_stars[1:end-1], δ_stars[2:end],   # Displacement thicknesses
                             Hs[1:end-1], Hs[2:end],             # Momentum shape parameters
                             H_ks[1:end-1], H_ks[2:end],         # Whitfield shape parameters
                             M_es[1:end-1], M_es[2:end],         # Mach numbers
                             Re_θs[1:end-1], Re_θs[2:end],       # Local Reynolds numbers
                             n_crit)                             # Critical amplification ratio

    R1      = getindex.(Rs, 1)                              # Momentum residual
    R2      = getindex.(Rs, 2)                              # Kinetic energy residual
    R3      = getindex.(Rs, 3)                              # Amplification/shear stress residual

    R1, R2, R3
end