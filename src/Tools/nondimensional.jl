module NonDimensional

using LinearAlgebra: norm

export dynamic_pressure, force_coefficient, moment_coefficient, rate_coefficient, pressure_coefficient, aerodynamic_coefficients, print_dynamics, reynolds_number

"""
Computes the dynamic pressure given density ρ and speed V.
"""
dynamic_pressure(ρ, V) = 0.5 * ρ * V^2

"""
Computes the non-dimensional force coefficient corresponding to standard aerodynamics.
"""
force_coefficient(force, q, S) = force / (q * S)

"""
Computes the non-dimensional moment coefficient corresponding to standard flight dynamics.
"""
moment_coefficient(moment, q, S, ref_length) = force_coefficient(moment, q, S) / ref_length

"""
Computes the non-dimensional angular velocity coefficient corresponding to standard flight dynamics.
"""
rate_coefficient(angular_speed, V, ref_length) = angular_speed * ref_length / 2V

"""
Computes the pressure coefficient given force, reference density, speed and area.
"""
pressure_coefficient(force, ρ, V, S) = force_coefficient(force, dynamic_pressure(ρ, V), S)

"""
Computes the incompressible pressure coefficient given a magnitude and a velocity vector.
"""
pressure_coefficient(mag, vels) = 1 - norm(vels)^2 / mag^2

"""
Prints the relevant aerodynamic/flight dynamics information.
"""
function aerodynamic_coefficients(force, moment, Ω, V, S, b, c, ρ)
    q = dynamic_pressure(ρ, V)
    CDi = force_coefficient(force[1], q, S)
    CY = force_coefficient(force[2], q, S)
    CL = force_coefficient(force[3], q, S)

    Cl = moment_coefficient(moment[1], q, S, b)
    Cm = moment_coefficient(moment[2], q, S, c)
    Cn = moment_coefficient(moment[3], q, S, b)

    p̄ = rate_coefficient(Ω[1], V, b)
    q̄ = rate_coefficient(Ω[2], V, c)
    r̄ = rate_coefficient(Ω[3], V, b)

    CL, CDi, CY, Cl, Cm, Cn, p̄, q̄, r̄
end

"""
Prints the relevant aerodynamic/flight dynamics information with a separate entry for drag.
"""
function aerodynamic_coefficients(force, moment, drag, Ω, V, S, b, c, ρ)
    q = dynamic_pressure(ρ, V)
    CL, CDi, CY, Cl, Cm, Cn, p̄, q̄, r̄ = aerodynamic_coefficients(force, moment, Ω, V, S, b, c, ρ)
    CDi = force_coefficient(drag, q, S)

    CL, CDi, CY, Cl, Cm, Cn, p̄, q̄, r̄
end

function print_dynamics(CL, CDi, CY, Cl, Cm, Cn, p̄, q̄, r̄)
    # println("Lift Coefficient")
    println("CL: $CL")
    # println("Drag Coefficient")
    println("CDi: $CDi")
    # println("Side Force Coefficient")
    println("CY: $CY")
    # println("Lift-to-Drag Ratio")
    println("L/D: $(CL/CDi)")
    # println("Rolling Moment Coefficient")
    println("Cl: $Cl")
    # println("Pitching Moment Coefficient")
    println("Cm: $Cm")
    # println("Yawing Moment Coefficient")
    println("Cn: $Cn")
    # println("Rolling Rate Coefficient")
    println("p̄: $p̄")
    # println("Pitching Rate Coefficient")
    println("q̄: $q̄")
    # println("Yawing Rate Coefficient")
    println("r̄: $r̄")
end

reynolds_number(ρ, V, c, μ) = ρ * V * c / μ
reynolds_number(V, c, ν) = V * c / ν

end