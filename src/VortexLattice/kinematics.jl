using Rotations
using LinearAlgebra

linear_momentum_body_axes(F, g, mass, U, U̇, Ω) = F + mass * g - mass * (U̇ + Ω × U)
angular_momentum_body_axes(M, I, Ω, Ω̇, h) = M - I * Ω̇ + Ω × (I * Ω + h)