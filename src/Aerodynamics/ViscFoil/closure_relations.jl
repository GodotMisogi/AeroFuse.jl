## CLOSURE RELATIONS
#============================================#

# Momentum shape parameter, H
whitfield_shape_parameter(H, Me) = (H - 0.290Me^2) / (1 + 0.113Me^2)

## Laminar closure for walls - Falkner-Skan profiles
#============================================#

# Kinetic energy shape factor, H*_lam
function laminar_kinetic_energy_shape_factor(Hk) 
    if Hk < 4.35
        Hs = (0.0111(Hk - 4.35)^2 - 0.0278(Hk - 4.35)^3) / (Hk + 1) - 0.0002((Hk - 4.35) * Hk)^2 + 1.528
    else
        Hs = 0.015 * (Hk - 4.35)^2 / Hk + 1.528
    end

    return Hs
end

# Skin-friction coefficient, cf_lam
laminar_skin_friction_coefficient(Hk, Reθ) = (ifelse(Hk < 5.5, 0.0727 * (5.5 - Hk)^3 / (Hk + 1), 0.015 * (1 - 1. / (Hk - 4.5))^2) - 0.07) / Reθ

## Laminar closure for wakes
#============================================#

# Kinetic energy shape factor, H*_lam,wake
laminar_wake_kinetic_energy_shape_factor(Hk) = 1.528 + ifelse(Hk < 3.5, 0.025(3.5 - Hk)^3 + 0.001(3.5 - Hk)^5, 0.015(Hk - 3.5)^2 / Hk)

## Dissipation coefficient relations
#============================================#

# Dissipation coefficient, cDi_lam
function laminar_dissipation_coefficient(Hk, Reθ) 
    if Hk < 4.
        cDi_Reθ = 0.00205(4 - Hk)^5.5 + 0.207
    else 
        cDi_Reθ = -0.0016(Hk - 4.)^2 / (1. + 0.02(Hk - 4.)^2) + 0.207
    end

    return cDi_Reθ / Reθ
end

# Dissipation coefficient due to stress, cDi_lam,stress
laminar_stress_dissipation_coefficient(Hs, Us, Reθ) = 2 * 0.15(0.995 - Us) / (Hs * Reθ)

# Dissipation coefficient, cDi_lam,wake
laminar_wake_dissipation_coefficient(Hk, Hs, Reθ) = 2.2 * (1 - 1 / Hk)^2 * (1 / Hk) / (Hs * Reθ)

# Turbulent outer dissipation coefficient, cDi_outer
turbulent_outer_dissipation_coefficient(Hs, Us, cτ) = cτ * (0.995 - Us) * 2 / Hs

# Turbulent wall dissipation coefficient, cDi_wall
function turbulent_wall_dissipation_coefficient(Hk, Hs, Reθ, cf, Us)
    H_min = 1 + 2.1 / log(Reθ)
    aa    = tanh((Hk - 1) / (H_min - 1))
    fac   = (1 + aa) / 2

    cDi   = cf / 2 * Us * (2 / Hs) * fac
end

function turbulent_dissipation_coefficient(cτ, Hk, Hs, Reθ, cf, Us, wake = false)
    cDi_lam = 0.
    cDi_turbwall = 0.
    if wake
        cDi = laminar_wake_dissipation_coefficient(Hk, Hs, Reθ)
    else
        cDi_turbwall = turbulent_wall_dissipation_coefficient(Hk, Hs, Reθ, cf, Us)

        cDi_outer = turbulent_outer_dissipation_coefficient(Hs, Us, cτ)

        cDi_lam = laminar_dissipation_coefficient(Hk, Reθ)

        cDi = cDi_outer + cDi_lam + cDi_turbwall
    end

    cDi_lamstress = laminar_stress_dissipation_coefficient(Hs, Us, Reθ)

    cDi += cDi_lamstress

    cDi = ifelse(cDi_lam > cDi, cDi_lam, cDi)
    cDi = ifelse(wake, 2cDi, cDi)
end

## Turbulent closure relations
#============================================#

# Skin-friction coefficient, cf_turb
function turbulent_skin_friction_coefficient(Hk, Reθ, Me, γ = 1.4) 
    Fc = √isentropic_temperature_ratio(Me, γ)

    # Smooth limiting
    aa = -1.33Hk
    aa = ifelse(aa < -17, -20 + 3exp((aa + 17) / 3), aa)

    cf0 = 0.3exp(aa) * log10(Reθ / Fc)^(-1.74 - 0.31Hk) / Fc + 1.1e-4(tanh(4 - Hk / 0.875) - 1)
end

# Density shape parameter, H**_turb (negligible for incompressible flow)
turbulent_density_shape_parameter(Hk, Me) = (0.064 / (Hk - 0.8) + 0.251) * Me^2

# Y+
y_plus(u_t, μ_e, ρ_e, η) = ρ_e * u_t * η / μ_e

# Kinetic energy shape parameters, H*_turb
function turbulent_kinetic_energy_shape_parameter(Hk, Reθ)
    Hs_min  = 1.5
    dHs_inf = 0.015

    H0  = ifelse(Reθ < 400, 4, 3 + 400/Reθ)
    Reβ = ifelse(Reθ < 200, 200, Reθ)

    Hs = Hs_min + 4 / Reβ + 
        ifelse(Hk < H0,
            (2 - Hs_min  - 4 / Reβ) * ((H0 - Hk) / H0 - 1)^2 * 1.5 / (Hk + 0.5),
            (Hk - H0)^2 * (0.007 * log(Reβ) / (Hk - H0 + 4 / Reβ)^2 + dHs_inf / Hk)
        )
end


# Normalized wall slip speed, Us
function normalized_wall_slip_speed(Hs, Hk, H, β, wake = false)
    # Limits
    Hk = ifelse(wake && Hk < 1.00005, 1.00005, Hk)
    Hk = ifelse(!wake && Hk < 1.05, 1.05, Hk)

    Us = Hs / 2 * (1 - 1 / β * (Hk - 1) / H)

    # WTF hax0rz
    Us = ifelse(!wake && Us > 0.95, 0.98, Us)
    Us = ifelse(wake && Us > 0.99995, 0.99995, Us)
end

## Boundary layer displacement thickness from Green's correlation, δ
function displacement_thickness(Hk, θ, δ_star) 
    δ = θ * (3.15 + 1.72 / (Hk - 1)) + δ_star
    return ifelse(δ > 12θ, 12θ, δ)
end 

# Equilibirum shear stress coefficient root, √cτ_EQ
function equilibrium_shear_stress_coefficient(H, Hk, Hs, Us, Reθ,  GA, GB, GC, wake = false)
    
    CC = 1 / (2 * GA^2 * GB)

    if wake
        Hkc = ifelse(Hk < 1.00005, 1.00005, Hk) - 1
    else 
        Hkc = ifelse(Hk < 1.05, 1.05, Hk) - 1 - GC / Reθ
    end

    Hkc = ifelse(Hkc < 0.01, 0.01, Hkc)

    cτ_EQ = (CC * Hs * (Hk - 1) * Hkc^2 / ((1 - Us) * H * Hk^2))
end

# Equilibrium uq quantity defined by Fidkowski, 1/Ue * dUe/dx
function equilibrium_uq(δs, cf, Hk, Reθ, GA, GB, GC, ηD, wake = false)
    if wake
        A  = GA * ηD
        B  = GB
        C  = 0.
        Hk = ifelse(Hk < 1.00005, 1.00005, Hk)
    else
        A  = GA
        B  = GB
        C  = GC
        Hk = ifelse(Hk < 1.05, 1.05, Hk)
    end

    Hkc = Hk - 1 - C / Reθ
    Hkc = ifelse(Hkc < 0.01, 0.01, Hkc)
    
    uq = (cf / 2 - (Hkc / (A * Hk))^2) / (B * δs)
end

## Transition
#============================================#

function amplification_rate(ñ, Hk, Reθ, θ, ñ_crit)
    Hk      = ifelse(Hk < 1.05, 1.05, Hk)
    Hmi     = 1 / (Hk - 1)
    lrc     = 2.291 * Hmi^0.43 + 0.7 * (tanh(14Hmi - 9.24) + 10)
    logr    = log10(Reθ)
    dl      = 0.1
    dñdξ    = 0.
    if (logr >= lrc)
        rn   = (logr - (lrc - dl)) / (2 * dl)
        rf   = ifelse(rn >= 1, 1., 3 * rn^2 - 2 * rn^3)
        da   = 0.028(Hk - 1) - 0.0345exp(-(3.87Hmi - 2.52)^2)
        af   = -0.05 + 2.7Hmi - 5.5Hmi^2 + 3Hmi^3 + 0.1exp(-20Hmi)
        dñdξ += rf * af * da * th
    end

    dñdξ += (1 + tanh(5. * (ñ - ñ_crit)) ) * 0.01 / θ

    return dñdξ
end

function shear_stress_coefficient(Hk, cτ_C, cτ_E, cτ_EQ)
    # Limits
    Hk = ifelse(Hk < 1.05, 1.05, 0.)

    # Shear stress coefficient
    a = cτ_C^2 * exp(-2cτ_E / (Hk - 1))

    cτ_tr = a * cτ_EQ
end

## Dispatch methods
#============================================#

# get_shape_parameter(Hk) = laminar_kinetic_energy_shape_factor(Hk)
# get_shape_parameter(Hk, Reθ) = turbulent_kinetic_energy_shape_parameter(Hk, Reθ)

get_shape_parameter(U) = U[2] / U[1]

# get_kinetic_energy_shape_parameter() = kinetic

function get_equilibrium_shear_stress_coefficient(U, Me, param; wake = false)
    GA, GB, GC = param
    H   = get_shape_parameter(U)
    Hk  = whitfield_shape_parameter(H, Me)
    Hs  = get_kinetic_energy_shape_parameter(U, param)
    Reθ = get_reynolds_number(U, param)
    Us  = get_wall_slip_speed(U, param)
    
    return equilibrium_shear_stress_coefficient(H, Hk, Hs, Us, Reθ,  GA, GB, GC, wake)
end

# Laminar quantity correlations
function laminar_body_quantities(Hk, Reθ)
    Hs  = laminar_kinetic_energy_shape_factor(Hk)
    cf  = laminar_skin_friction_coefficient(Hk, Reθ)
    cDi = laminar_dissipation_coefficient(Hk, Reθ)

    Hs, cf, cDi
end

# Turbulent quantity correlations
function turbulent_body_quantities(cτ, Hk, H, Me, Reθ, GA, GB, GC, wake = false)
    Hs    = turbulent_kinetic_energy_shape_parameter(Hk, Reθ)
    cf    = turbulent_skin_friction_coefficient(Hk, Reθ, Me)
    Us    = normalized_wall_slip_speed(Hs, Hk, H, 1. / GB)
    cτ_EQ = equilibrium_shear_stress_coefficient(H, Hk, Hs, Us, Reθ, GA, GB, GC, wake)
    cDi   = turbulent_dissipation_coefficient(cτ, Hk, Hs, Reθ, cf, Us, wake)

    Hs, cf, cDi, cτ_EQ, Us
end

## Generic conversions
#============================================#

# Squire-Young formula for drag coefficient
squire_young_drag_coefficient(θ_wake, u_e_wake, H_wake, V_inf) = 2θ_wake * (u_e_wake / V_inf)^((5 + H_wake) / 2)

speed_of_sound(g, H0, Uk) = √((g - 1) * (H0 - Uk^2 / 2))