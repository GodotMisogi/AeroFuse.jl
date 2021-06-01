## Influence coefficients for velocities
#============================================#

# Helper
function transform(x, z, x1, x2) 
    s1, s2 = x - x1, x - x2	
    r1, r2 = √(s1^2 + z^2), √(s2^2 + z^2)
    θ1, θ2 = atan(z, s1), atan(z, s2)

    r1, r2, θ1, θ2
end

term_1(xp, x, z, x1, x2, r1, r2, θ1, θ2) = (x2 - x1) - z * (θ2 - θ1) + (x - xp) * log(r2 / r1)
term_2(xp, x, z, r1, r2, θ1, θ2) = (x - xp) * (θ2 - θ1) + z * log(r2 / r1)

boundary_term(xp, x, x1, x2) = (1/2π, 1/2 * (x - xp) / (x2 - x1))

# Source velocities
#============================================#

function constant_source_velocity(σ, x, z, x1, x2)
    r1, r2, θ1, θ2 = transform(x, z, x1, x2)
    @. σ / 2π * SVector(-log(r2 / r1), θ2 - θ1)
end

function linear_source_velocity_a(σ1, x, z, x1, x2)
    r1, r2, θ1, θ2 = transform(x, z, x1, x2)
    if abs(z) <= 1e-12
        ua, wa = boundary_term(x2, x, x1, x2)
        σ1 .* SVector(ua, -wa)
    else 
        ua = term_1(x2, x, z, x1, x2, r1, r2, θ1, θ2)
        wa = term_2(x2, x, z, r1, r2, θ1, θ2)
        
        @. σ1 / (2π * (x2 - x1)) * SVector(ua, -wa)
    end
end

function linear_source_velocity_b(σ2, x, z, x1, x2)
    r1, r2, θ1, θ2 = transform(x, z, x1, x2)
    if abs(z) <= 1e-12
        ub, wb = boundary_term(x1, x, x1, x2)
        σ2 .* SVector(-ub, wb)
    else 
        ub = term_1(x1, x, z, x1, x2, r1, r2, θ1, θ2)
        wb = term_2(x1, x, z, r1, r2, θ1, θ2)
        
        @. σ2 / (2π * (x2 - x1)) * SVector(-ub, wb)
    end
end

# Vortex velocities
#============================================#

function linear_vortex_velocity_a(γ1, x, z, x1, x2) 
    r1, r2, θ1, θ2 = transform(x, z, x1, x2)
    # if z == 0.
    # 	ua, wa = (reverse ∘ boundary_term)(x2, x, x1, x2)
    # 	γ1 .* SVector(-ua, -wa)
    # else
        ua = term_2(x2, x, z, r1, r2, θ1, θ2)
        wa = term_1(x2, x, z, x1, x2, r1, r2, θ1, θ2)

        @. γ1 / (2π * (x2 - x1)) * SVector(-ua, -wa)
    # end
end

function linear_vortex_velocity_b(γ2, x, z, x1, x2)
    r1, r2, θ1, θ2 = transform(x, z, x1, x2)
    # if z == 0.
    # 	ub, wb = (reverse ∘ boundary_term)(x1, x, x1, x2)
    # 	γ2 .* SVector(ub, wb)
    # else
        ub = term_2(x1, x, z, r1, r2, θ1, θ2)
        wb = term_1(x1, x, z, x1, x2, r1, r2, θ1, θ2)

        @. γ2 / (2π * (x2 - x1)) * SVector(ub, wb)
    # end
end