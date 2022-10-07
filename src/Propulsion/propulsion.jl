module Propulsion

using Roots
using StaticArrays

## Propulsion models
#==========================================================================================#

abstract type AbstractPropulsion{T <: Real} end

## Actuator disc theory
#=======================#

actuator_disc_induced_velocity(T, V_perp, A_disk, ρ, κ) = V_perp + κ / 2 * (-V_perp + √(V_perp^2 + 2T / (ρ * A_disk)))

momentum_power(T, V_perp, A_disk, ρ, κ) = T * actuator_disc_induced_velocity(T, V_perp, A_disk, ρ, κ)

profile_power_coefficient(σ, Cd0p, μ) = σ * Cd0p / 8 * (1 + 4.6μ^2)

propeller_normal_force(σ_e, β, f, q_perp, A_disk, α_in) = (4.25σ_e * sin(β + deg2rad(8)) * f * q_perp * A_disk / (1 + 2σ_e)) * tan(α_in)

solidity(B, c_b, R) = 2B * c_b / (3π * R)
thrust_factor(C_T) = 1 + (√(1 + C_T) - 1) / 2 + C_T / (4 * (2 + C_T))
thrust_coefficient(T, q, A_disc) = T / (q * A_disc)

struct ActuatorDisc{T <: Real}
    thrust :: T
    power  :: T
    speed  :: T
    area   :: T
    radius :: T
    density :: T
end

induced_velocity(disc :: ActuatorDisc) = actuator_disc_induced_velocity(disc.thrust, disc.speed, disc.area, disc.density, κ)

momentum_power(disc :: ActuatorDisc) = momentum_power(disc.thrust, disc.speed, disc.area, disc.density, κ)

# struct ActuatorPropulsion{T}
#     model :: AbstractPropulsion{T}
#     throttle :: T
# end

## Blade-element momentum theory
#================================#

induced_speed(Vx, Vy, ax, ay) = √((Vx * (1 + ax))^2 + (Vy * (1 - ay))^2)

inflow_angle(φ0, Vx, Vy, Cl, Cd, F, σ) = find_zero(φ0) do φ
    (sφ, cφ) = sincos(φ)
    # Coordinate transformation
    Cx = cφ * Cd - sφ * Cl
    Cy = sφ * Cd + cφ * Cl
    # Residual equation
    4F * sφ * (Vy * (sφ - Cx * σ) - Vx * (cφ + Cy * σ))
end

blade_solidity(B, c_b, R) = B * c_b / (2π * R)

slipstream_contraction(x, a_x, R) = R * √((1 + a_x) / (1 + a_x * (1 + x / √(R^2 + x^2))))

function induced_velocity(x, r, a_x, Vx0)
    R = slipstream_contraction(x, a_x, R)

    t1 = R^2 - x^2 - r^2
    a = √((√(t1^2 + 4R^2 * x^2) + t1) / 2R^2)
    b = √(x^2 + (R + r)^2) + √(x^2 + (R - r)^2)

    t2 = √(R^2 - r^2)
    Vx_r0 = Vx0 * √t2  / R

    V_r = Vx_r0 * (abs(x) * R / (2r * t2) * (1 / a - a) - r / (2 * t2) * asin(2R / b))
    V_x = Vx_r0 * (2 - a * R / t2 + x / t2 * asin(2R / b))

    return SVector(V_r, V_x)
end

end