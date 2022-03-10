## INTEGRAL BOUNDARY LAYER EQUATIONS - FINITE DIFFERENCE
#======================================================================#

# Boundary layer equations
#============================================#

# Momentum residual equation
momentum_finite_difference(dlogθ, dlogUe, dlogξ, H, cfxθ, Me) = dlogθ + (H + 2 - Me^2) * dlogUe - dlogξ * cfxθ / 2

# Shape parameter residual equation
shape_finite_difference(dlogHs, dlogUe, dlogξ, H, Hs, Hss, cfxθ, cDixθ) = dlogHs + (2Hss / Hs + 1 - H) * dlogUe + dlogξ * (cfxθ / 2 - 2cDixθ / Hs)

## Laminar/turbulence closure equations
#============================================#

stress_lag_transport_finite_difference(dlogUe, dlog_sqrt_cτ, Δξ, δ, uq, ηD, sqrt_cτ, sqrt_cτ_EQ, a) = 2δ * (dlog_sqrt_cτ + dlogUe - uq * Δξ) - a * (sqrt_cτ_EQ - ηD * sqrt_cτ) * Δξ

# Amplification ratio equation for Tollmien-Schlichting waves using envelope eⁿ, ñ
# amplification_finite_difference(Δn, Δξ, Hk, θ) = Δn - (dñdReθ(Hk) * (m(Hk) + 1) / 2 * l(Hk) / θ) * Δξ

## RESIDUAL EVALUATIONS
#======================================================================#

## Boundary layer residuals
#============================================#

# Transition station
function residual_transition(bl_lam, bl_turb, x1, x2, params)

    # Unpack parameters
    γ, μ0, ρ0, h0, Tsrat, KT_β, KT_λ, Re, M∞, V∞, K_lag, GA, GB, GC, ηD, cτ_C, cτ_E, ñ_crit = params

    dx = x2 - x1

    function transition_residual!(xt)
        xt     = xt[1]

        # Upwinding at transition
        upw   = (xt - x1) / dx
        ñt    = upwind(upw, bl_lam.tu, bl_turb.tu) 
        Hkt   = upwind(upw, bl_lam.Hk, bl_turb.Hk)
        Reθt  = upwind(upw, bl_lam.Reθ, bl_turb.Reθ)
        θt    = upwind(upw, bl_lam.θ, bl_turb.θ)

        # Calculate amplification rates 
        dñ1dξ = amplification_rate(bl_lam.tu, bl_lam.Hk, bl_lam.Reθ, bl_lam.θ, ñ_crit)
        dñtdξ = amplification_rate(ñt, Hkt, Reθt, θt, ñ_crit)

        # Calculate amplification factor residual
        R = ñ_crit - bl_lam.tu - (xt - x1) / 2 * (dñ1dξ + dñtdξ)
    end
    
    xt = x1 + dx / 2 # Initial guess for transition location

    # Evaluate transition location using Newton iteration
    res_trans = nlsolve(transition_residual!, [xt],
                        iterations = 20,
                        method     = :newton,
                        autodiff   = :forward,
                        show_trace = true
                       )

    # Get transition state
    @show xt = res_trans.zero[1]

    # Upwinding to transition location
    upw  = (xt - x1) / dx
    Uet  = upwind(upw, bl_lam.Ue, bl_turb.Ue)
    Hkt  = upwind(upw, bl_lam.Hk, bl_turb.Hk)
    θt   = upwind(upw, bl_lam.θ, bl_turb.θ)
    Reθt = upwind(upw, bl_lam.Reθ, bl_turb.Reθ)
    δst  = upwind(upw, bl_lam.δs, bl_turb.δs)
    cτ_EQt = upwind(upw, bl_lam.cτ_EQ, bl_turb.cτ_EQ)

    cτ_tr = shear_stress_coefficient(Hkt, cτ_C, cτ_E, cτ_EQt)

    bl_lam_tran = BoundaryLayer2D(θt, δst, ñ_crit, Uet, params)
    bl_turb_tran = BoundaryLayer2D(θt, δst, √cτ_tr, Uet, params)

    # Laminar region
    Rt_lam = residual_station(bl_lam, bl_lam_tran, x1, xt, params; turb = false)

    # Turbulent region
    Rt_turb = residual_station(bl_turb_tran, bl_turb, xt, x2, params; turb = true)

    # Combine residuals
    return R = Rt_lam + Rt_turb
end

# Evaluate residuals
function residual_station(bl1 :: BoundaryLayer2D, bl2 :: BoundaryLayer2D, x1, x2, params; turb = false, simi = false, wake = false)

    # Unpack parameters
    γ, μ0, ρ0, h0, Tsrat, KT_β, KT_λ, Re, M∞, V∞, K_lag, GA, GB, GC, ηD, cτ_C, cτ_E, ñ_crit = params

    # Differences
    #=================================#

    # Similarity station ⇒ bl1 = bl2, x1 = x2
    if simi
        dlogθ  = 0.
        dlogHs = 0.
        dlogUe = 1.
        dlogξ  = 1.
        Δξ     = upwind(0.5, x1, x2)
    else
        dlogθ  = log(bl2.θ / bl1.θ)
        dlogHs = log(bl2.Hs / bl1.Hs)
        dlogUe = log(bl2.Ue / bl1.Ue)
        dlogξ  = log(x2 / x1)
        Δξ     = x2 - x1
    end

    # Wake gap shape parameter, Hw
    #=================================#

    # WRONG TREATMENT
    # Hw1 = wgap / bl1.θ
    # Hw2 = wgap / bl2.θ
    # Hw  = upwind(0.5, Hw1, Hw2)

    # Amplification or shear stress residual
    #=================================#

    Hss   = 0.
    R_lag = 0.

    # Upwinding
    upw = upwind_factor(bl1.Hk, bl2.Hk, wake)
    cf  = upwind(upw, bl1.cf, bl2.cf)
    H   = upwind(upw, bl1.H, bl2.H)
    Hk  = upwind(upw, bl1.Hk, bl2.Hk)
    Hs  = upwind(upw, bl1.Hs, bl2.Hs)
    Hss = upwind(upw, bl1.Hss, bl2.Hss)
    Me  = upwind(upw, bl1.Me, bl2.Me)
    ξ   = upwind(upw, x1, x2)
    θ   = upwind(upw, bl1.θ, bl2.θ)
    Reθ = upwind(upw, bl1.Reθ, bl2.Reθ)

    # Special treatment for cf
    cfxθ1 = bl1.cf * x1 / bl1.θ
    cfxθ2 = bl2.cf * x2 / bl2.θ
    blm   = BoundaryLayer2D((bl1.θ  + bl2.θ ) / 2,
                            (bl1.δs + bl2.δs) / 2,
                            (bl1.tu + bl2.tu) / 2,
                            (bl1.Ue + bl2.Ue) / 2,
                            params)

    cfxθm = blm.cf * ξ / blm.θ

    cfxθ  = cfxθ1 / 4 + cfxθm / 2 + cfxθ2 / 4
    
    # Upwinding
    cfxθu = upwind(upw, cfxθ1, cfxθ2)
    cDixθ = upwind(upw, bl1.cDi * x1 / bl1.θ, bl2.cDi * x2 / bl2.θ)

    # Laminar/turbulent residuals
    #=================================#

    if turb # Turbulent
        # Shear stress coefficient
        dlog_sqrt_cτ = log(bl2.tu / bl1.tu)
        sqrt_cτ      = upwind(upw, bl1.tu, bl2.tu)
        sqrt_cτ_EQ   = upwind(upw, √bl1.cτ_EQ, √bl2.cτ_EQ)
    
        # Boundary layer displacement thickness
        δs = upwind(upw, bl1.δs, bl2.δs)
        δ  = displacement_thickness(Hk, θ, δs)

        # Equilibirum 1/Ue * dUe/dx
        Us = upwind(upw, bl1.Us, bl2.Us)
        uq = equilibrium_uq(δs, cf, Hk, Reθ, GA, GB, GC, ηD)
        a  = K_lag / (GB * (1 + Us))

        # Shear stress lag transport equation residual
        R_lag = 2δ * (dlog_sqrt_cτ + dlogUe - uq * Δξ) - a * (sqrt_cτ_EQ - ηD * sqrt_cτ) * Δξ
    else # Laminar
        if simi
            R_lag = bl1.tu + bl2.tu
        else 
            # Evaluate amplification rates and upwind
            dñ1dξ = amplification_rate(bl1.tu, bl1.Hk, bl1.Reθ, bl1.θ, ñ_crit)
            dñ2dξ = amplification_rate(bl2.tu, bl2.Hk, bl2.Reθ, bl2.θ, ñ_crit)
            dñdξ  = upwind(upw, dñ1dξ, dñ2dξ)

            # Amplification rate equation residual
            R_lag = bl2.tu - bl1.tu - dñdξ * Δξ
        end
    end

    # Momentum thickness and shape parameters residuals
    #=================================#

    R_mom = dlogθ + (H + 2 - Me^2) * dlogUe - dlogξ * cfxθ / 2

    R_shape = dlogHs + (2Hss / Hs + 1 - H) * dlogUe + dlogξ * (cfxθ / 2 - 2cDixθ / Hs)

    [R_mom, R_shape, R_lag]
end

# Augmented displacement layer thicknesses
# residual_station(bl1, bl2, x1, x2, Δδs1, Δδs2, params; turb = false, simi = false, wake = false) = residual_station(Ue1, Ue2, θ1, θ2, p1, p2, x1, x2, δs1 - Δδs1, δs2 - Δδs2, params; turb = turb, simi = simi, wake = wake)

## Wake residuals
#============================================#

function wake_system(θ, δs, cτ, bl_l, bl_u, hTE, params)
    cτ_C, cτ_E = params

    cτ_u = shear_stress_coefficient(bl_u.Hk, cτ_C, cτ_E, bl_u.cτ_EQ)
    cτ_l = shear_stress_coefficient(bl_l.Hk, cτ_C, cτ_E, bl_l.cτ_EQ)

    # θ = upwind(0.5, θ_u, θ_l)
    cτ_wake = (cτ_u * bl_u.θ + cτ_l * bl_l.θ) / (bl_u.θ + bl_l.θ)

    # Match momentum thicknesses
    R1 = θ - (bl_u.θ + bl_l.θ)

    # Match displacement thicknesses including trailing edge gap
    R2 = δs - (bl_u.δs + bl_l.δs + hTE)

    # Match shear stress coefficients in the wake
    R3 = cτ - cτ_wake

    [R1, R2, R3]
end

## Initialization
#============================================#

function thwaites_initialization(K, ν) 
    θ = √(0.45ν /6K)
    δ_st = 2.2θ

    θ, δ_st
end

function initialize_boundary_layer()

end