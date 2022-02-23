
## Vorticity-streamfunction formulation
#===========================================================================#

# Somewhat unnecessary, as the point is transformed to coordinates with the first point of the panel as the origin
function transform(x, z, x1, x2)
    s1, s2 = x - x1, x - x2

    r1, r2 = √(s1^2 + z^2), √(s2^2 + z^2)
    θ1, θ2 = atan(z, s1), atan(z, s2)

    x, z, r1, r2, θ1, θ2
end

function linear_vortex_stream_minus(γ, x, z, x1, x2)
    a, h, r1, r2, θ1, θ2 = transform(x, z, x1, x2)
    d = x2 - x1

    # Singularity checks
    ϵ = 1e-10
    if r1 < ϵ; logr1 = 0; else logr1 = log(r1) end
    if r2 < ϵ; logr2 = 0; else logr2 = log(r2) end

    ψ_m = γ / 2π * (h * (θ2 - θ1) - d + a * logr1 - (a - d) * logr2)
end

function linear_vortex_stream_plus(γ, x, z, x1, x2)
    a, _, r1, r2, _, _= transform(x, z, x1, x2)
    d   = x2 - x1
    ψ_m = linear_vortex_stream_minus(1., x, z, x1, x2)

    ϵ = 1e-10
    if r1 < ϵ; logr1 = 0; else logr1 = log(r1) end
    if r2 < ϵ; logr2 = 0; else logr2 = log(r2) end
    
    ψ_p = γ / d * (a * ψ_m + 1/4π * (r2^2 * (logr2 - 1/2) - r1^2 * (logr1 - 1/2)))
end


function constant_source_stream(σ, x, n, x1, x2)
    a, h, r1, r2, θ1, θ2 = transform(x, n, x1, x2)
    d = x2 - x1

    # Singularity checks
    ϵ = 1e-9;
    if r1 < ϵ; logr1 = 0; θ1 = π; θ2 = π; else logr1 = log(r1); end
    if r2 < ϵ; logr2 = 0; θ1 = 0; θ2 = 0; else logr2 = log(r2); end

    ψ = σ / 2π * (a * (θ1 - θ2) + d * θ2 + h * (logr1 - logr2)) 

    # Branch cuts
    dψ = x2
    if (θ1 + θ2) > π
      ψ -= 0.25 * dψ
    else
      ψ += 0.75 * dψ
    end
end

# two_piece_source(λ1, λ2, r1, r2, r) = (r - r1) / (r2 - r1) * λ1 + (r2 - r) / (r2 - r1) * λ2

# function source_stream_plus(str, x, n, x_p, x_q)
#   s_p, s_q, r_p, r_q, θ_p, θ_q = transform(x, n, x_p, x_q)
#   str / 4π * (x_p * θ_p - x_q * θ_q + n * log(r_p / r_q))
# end

# function source_stream_minus(str, x, n, x_p, x_q)
#   s_p, s_q, r_p, r_q, θ_p, θ_q = transform(x, n, x_p, x_q)
#   λ_plus = source_stream_1(1., x, n, x_p, x_q)

#   str / (4π * (x_p - x_q)) * ((x_p + x_q) * λ_plus + r_p^2 * θ_p - r_q^2 * θ_q + n * (x_q - x_p))
# end