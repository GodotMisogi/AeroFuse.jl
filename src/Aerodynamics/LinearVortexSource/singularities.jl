## SINGULARITY DISTRIBUTIONS
#============================================#

## Helper functions

function polar_transform(x, z, x1, x2)
    s1, s2 = x - x1, x - x2
    r1, r2 = √(s1^2 + z^2), √(s2^2 + z^2)
    θ1, θ2 = atan(z, s1), atan(z, s2)

    r1, r2, θ1, θ2
end

# "Fundamental" terms
term_1(xp, x, z, d, r1, r2, θ1, θ2) = d - z * (θ2 - θ1) + (x - xp) * log(r2 / r1)
term_2(xp, x, z, r1, r2, θ1, θ2) = (x - xp) * (θ2 - θ1) + z * log(r2 / r1)

boundary_term(xp, x, d) = (1/2π, 1/2 * (x - xp) / d)

# VORTEX SINGULARITIES
#============================================#

## Linear-strength singularities

# Streamfunctions
function linear_vortex_stream_minus(γ, x, z, x1, x2)
    r1, r2, θ1, θ2 = polar_transform(x, z, x1, x2)
    d = x2 - x1

    # Singularity checks
    ϵ = 1e-10
    if r1 < ϵ; logr1 = 0.; else logr1 = log(r1) end
    if r2 < ϵ; logr2 = 0.; else logr2 = log(r2) end

    ψ_m = γ / 2π * (z * (θ2 - θ1) - d + x * logr1 - (x - d) * logr2)
end

function linear_vortex_stream_plus(γ, x, z, x1, x2)
    r1, r2, _, _= polar_transform(x, z, x1, x2)
    d   = x2 - x1
    ψ_m = linear_vortex_stream_minus(1., x, z, x1, x2)

    ϵ = 1e-10
    if r1 < ϵ; logr1 = 0.; else logr1 = log(r1) end
    if r2 < ϵ; logr2 = 0.; else logr2 = log(r2) end
    
    ψ_p = γ / d * (x * ψ_m + 1/4π * (r2^2 * (logr2 - 1/2) - r1^2 * (logr1 - 1/2)))
end

# Velocities
function linear_vortex_velocity_a(γ1, x, z, x1, x2)
    r1, r2, θ1, θ2 = polar_transform(x, z, x1, x2)
    d = x2 - x1

    # ua = term_2(x2, x, z, r1, r2, θ1, θ2)
    # wa = term_1(x2, x, z, d, r1, r2, θ1, θ2)

    dlogr = log(r1 / r2)
    dθ    = (θ2 - θ1)

    u =     dθ + (z * dlogr     - x * dθ) / d
    w = -dlogr + (x * dlogr - d + z * dθ) / d
    
    γ1 / 2π * SVector(u, w)
end

function linear_vortex_velocity_b(γ2, x, z, x1, x2)
    r1, r2, θ1, θ2 = polar_transform(x, z, x1, x2)
    d = x2 - x1

    # ub = term_2(x1, x, z, r1, r2, θ1, θ2)
    # wb = term_1(x1, x, z, d, r1, r2, θ1, θ2)

    dlogr = log(r1 / r2)
    dθ    = (θ2 - θ1)

    u = -(z * dlogr     - x * dθ) / d
    w = -(x * dlogr - d + z * dθ) / d

    γ2 / 2π * SVector(u, w)
end


# SOURCE SINGULARITIES
#============================================#

## Constant-strength singularities

# Streamfunctions
function constant_source_stream(σ, x, z, x1, x2)
    r1, r2, θ1, θ2 = polar_transform(x, z, x1, x2)
    d = x2 - x1

    # Singularity checks
    ϵ = 1e-10;
    if r1 < ϵ; logr1 = 0.; θ1 = π; θ2 = π; else logr1 = log(r1); end
    if r2 < ϵ; logr2 = 0.; θ1 = 0.; θ2 = 0.; else logr2 = log(r2); end

    ψ = σ / 2π * (x * (θ1 - θ2) + d * θ2 + z * (logr1 - logr2)) 

    # Branch cuts
    dψ = x2
    if (θ1 + θ2) > π
        ψ -= 0.25 * dψ
    else
        ψ += 0.75 * dψ
    end
end

# Velocity
function constant_source_velocity(σ, x, z, x1, x2)
    r1, r2, θ1, θ2 = polar_transform(x, z, x1, x2)
    @. σ / 2π * SVector(-log(r2 / r1), θ2 - θ1)
end

## Linear-strength singularities

# Streamfunctions
function linear_source_stream_minus(σ, x, z, x1, x2)
    r1, r2, θ1, θ2 = polar_transform(x, z, x1, x2)
    d = x2 - x1

    # Branch cuts
    if θ1 < 0.; θ1 += 2π; end
    if θ2 < 0.; θ2 += 2π; end

    # Singularity checks
    ϵ = 1e-10;
    if r1 < ϵ; logr1 = 0.; θ1 = π; θ2 = π; else logr1 = log(r1); end
    if r2 < ϵ; logr2 = 0.; θ1 = 0.; θ2 = 0.; else logr2 = log(r2); end

    ψ = σ / 2π * (x * (θ1 - θ2) + d * θ2 + z * (logr1 - logr2))
end

function linear_source_stream_plus(σ, x, z, x1, x2)
    r1, r2, θ1, θ2 = polar_transform(x, z, x1, x2)
    d = x2 - x1
    ψ_m = linear_source_stream_minus(1., x, z, x1, x2)

    # Branch cuts
    if θ1 < 0.; θ1 += 2π; end
    if θ2 < 0.; θ2 += 2π; end

    # Singularity checks
    ϵ = 1e-10;
    if r1 < ϵ; logr1 = 0.; θ1 = π; θ2 = π; else logr1 = log(r1); end
    if r2 < ϵ; logr2 = 0.; θ1 = 0.; θ2 = 0.; else logr2 = log(r2); end

    ψ_p = σ / d * (x * ψ_m + 1/4π * (r2^2 * θ2 - r1^2 * θ1 - z * d))
end

# Velocities
function linear_source_velocity_a(σ1, x, z, x1, x2)
    r1, r2, θ1, θ2 = polar_transform(x, z, x1, x2)
    d = x2 - x1
    # if abs(z) <= 1e-12
        # ua, wa = # boundary_term(x2, x, d)
        # σ1 .* SVector(ua, -wa)
    # else
        # ua = term_1(x2, x, z, d, r1, r2, θ1, θ2)
        # wa = term_2(x2, x, z, r1, r2, θ1, θ2)

        # @. σ1 / (2π * (x2 - x1)) * SVector(ua, -wa)
    # end

    dlogr = log(r1 / r2)
    dθ    = (θ2 - θ1)

    u = dlogr / 2π - (x * dlogr - d + z * dθ) / (2π * d)
    w =    dθ / 2π + (z * dlogr     - x * dθ) / (2π * d)

    σ1 * SVector(u, w)
end

function linear_source_velocity_b(σ2, x, z, x1, x2)
    r1, r2, θ1, θ2 = polar_transform(x, z, x1, x2)
    d = x2 - x1
    # if abs(z) <= 1e-12
    #     ub, wb = # boundary_term(x1, x, d)
    #     σ2 .* SVector(-ub, wb)
    # else
    #     ub = term_1(x1, x, z, d, r1, r2, θ1, θ2)
    #     wb = term_2(x1, x, z, r1, r2, θ1, θ2)

    #     @. σ2 / (2π * (x2 - x1)) * SVector(-ub, wb)
    # end

    dlogr = log(r1 / r2)
    dθ    = (θ2 - θ1)

    u =  (x * dlogr - d + z * dθ) / (2π * d)
    w = -(z * dlogr     - x * dθ) / (2π * d)

    σ2 * SVector(u, w)
end
