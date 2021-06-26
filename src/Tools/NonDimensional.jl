module NonDimensional

using LinearAlgebra: norm

# export dynamic_pressure, force_coefficient, moment_coefficient, moment_coefficients, rate_coefficient, pressure_coefficient, aerodynamic_coefficients, print_coefficients, reynolds_number, print_derivatives

"""
    dynamic_pressure(ρ, V)

Compute the dynamic pressure given density ``\\rho`` and speed ``V``.
"""
dynamic_pressure(ρ, V) = 1/2 * ρ * V^2

"""
    force_coefficient(force, q, S)

Compute the non-dimensional force coefficient given force, dynamic pressure ``q``, and area ``S``.
"""
force_coefficient(force, q, S) = force / (q * S)

"""
    moment_coefficient(moment, q, S, L)

Compute the non-dimensional moment coefficient given force, dynamic pressure ``q``, area ``S``, and reference length ``L``.
"""
moment_coefficient(moment, q, S, L) = force_coefficient(moment, q, S) / L

"""
    rate_coefficient(Ω, V, L)

Compute the non-dimensional angular velocity coefficient given angular speed ``\\Omega``, reference speed ``V``, and reference length ``L``.
"""
rate_coefficient(Ω, V, L) = Ω * L / 2V

"""
    pressure_coefficient(force, ρ, V, S)

Compute the pressure coefficient given force, reference density ``\\rho``, speed ``V`` and area ``S``.
"""
pressure_coefficient(force, ρ, V, S) = force_coefficient(force, dynamic_pressure(ρ, V), S)

pressure_coefficient(mag, vels) = 1 - (norm(vels) / mag)^2

moment_coefficient(moment, q, S, b, c) = moment_coefficient.(moment, q, S, [b, c, b])

rate_coefficient(Ω, V, b, c) = rate_coefficient.(Ω, V, [b, c, b])


"""
    aerodynamic_coefficients(force, moment, Ω, V, S, b, c, ρ)

Compute the relevant aerodynamic coefficients given a net force, moment and angular rates ``\\Omega`` with reference speed ``V``, area ``S``, span ``b``, chord ``c``, and density ``\\rho``.
"""
function aerodynamic_coefficients(force, moment, V, S, b, c, ρ)
    q = dynamic_pressure(ρ, V)

    [ force_coefficient(force, q, S) 		 ;
      moment_coefficient(moment, q, S, b, c) ]
end

reynolds_number(ρ, V, c, μ) = ρ * V * c / μ
reynolds_number(V, c, ν) = V * c / ν

force(CF, q, S) = CF * q * S
moment(CM, q, S, c) = CM * q * S * c
moment(CM, q, S, b, c) = moment.(CM, q, S, [b, c, b])

end