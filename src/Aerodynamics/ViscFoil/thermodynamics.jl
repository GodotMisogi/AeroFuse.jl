## THERMODYNAMIC RELATIONS
#============================================#

karman_tsien_factor(M) = (M / (1 + √(1 - M^2)))^2

karman_tsien_correction(V, V∞, λ) = V * (1 - λ) / (1 - λ * (V / V∞)^2)

isentropic_temperature_ratio(M, γ) = 1 / (1 + (γ - 1) / 2 * M^2)

isentropic_density_ratio(M, γ) = isentropic_temperature_ratio(M, γ)^(1 / (γ - 1))

isentropic_pressure_ratio(M, γ) = isentropic_density_ratio(M, γ)^(-γ)

function mach_number(Uk, h0, γ)
    # Local speed of sound
    c = √((γ - 1) * (h0 - Uk^2 / 2))

    M = Uk / c
end

sutherland_ratio(Tr, Ts) = Tr^1.5 * (1 + Ts) / (Tr + Ts)

function reynolds_number(δ, Uk, M, h0, γ, Ts, μ0, ρ0)
    # Edge/stagnation temperature ratio
    Tr = 1 - 0.6 * Uk^2 / h0

    # Sutherland's ratio
    f = sutherland_ratio(Tr, Ts)

    # Local dynamic viscosity
    μ = μ0 * f

    # Isentropic density scaling
    ρ = ρ0 / isentropic_density_ratio(M, γ)

    Reθ = ρ * Uk * δ / μ
end

reynolds_number(ρ, V, c, μ) = ρ * V * c / μ