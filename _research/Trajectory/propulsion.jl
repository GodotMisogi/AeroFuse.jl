abstract type AbstractPropulsion end

# Propulsion
actuator_disk_induced_velocity(T, V_perp, A_disk, ρ, κ) = V_perp + κ / 2 * (-V_perp + √(V_perp^2 + 2T / (ρ * A_disk)))

momentum_power(T, V_perp, A_disk, ρ, κ) = T * actuator_disk_induced_velocity(T, V_perp, A_disk, ρ, κ)

profile_power_coefficient(σ, Cd0p, μ) = σ * Cd0p / 8 * (1 + 4.6μ^2)

propeller_normal_force(σ_e, β, f, q_perp, A_disk, α_in) = (4.25σ_e * sin(β + deg2rad(8)) * f * q_perp * A_disk / (1 + 2σ_e)) * tan(α_in)

solidity(B, c_b, R) = 2B * c_b / (3π * R)
thrust_factor(C_T) = 1 + (√(1 + C_T) - 1) / 2 + C_T / (4 * (2 + C_T))
thrust_coefficient(T, q, A_disk) = force_coefficient(T, q, A_disk)

struct ActuatorDisk{T <: Real}
    thrust :: T
    power  :: T
    speed  :: T
    area   :: T
end

struct ActuatorPropulsion{T <: Real}
    model :: AbstractPropulsion
    throttle :: T
end