
## Vorticity-streamfunction formulation
#===========================================================================#

# Somewhat unnecessary, as the point is transformed to coordinates with the first point of the panel as the origin
function transform(x, n, x1, x2)
    s1, s2 = x - x1, x - x2

    # println("ses: $((s1, s2))")

    r1, r2 = √(s1^2 + n^2), √(s2^2 + n^2)
    β1, β2 = atan(s1, n), atan(s2, n)
    # r1, r2, β1, β2 = ifelse(r1 == 0 || r2 == 0, (1., 1., 0., 0.), (r1, r2, β1, β2))

    s1, s2, r1, r2, β1, β2
end

function vortex_stream_plus(str, x, n, x1, x2)
    s1, s2, r1, r2, β1, β2 = transform(x, n, x1, x2)
    str / 4π * (s1 * log(r1) - s2 * log(r2) + s2 - s1 + n * (β1 - β2))
end

function vortex_stream_minus(str, x, n, x1, x2)
    s1, s2, r1, r2, β1, β2 = transform(x, n, x1, x2)
    ψ_plus = vortex_stream_plus(1., x, n, x1, x2)

    str / (4π * (s1 - s2)) * ((s1 + s2) * ψ_plus + r2^2 * log(r2) - r1^2 * log(r1) + 0.5 * (s1^2 - s2^2))
end

function source_stream_trailing_edge(str, x, n, x1, x2)
    s1, s2, r1, r2, β1, β2 = transform(x, n, x1, x2)
    θ1, θ2 = atan(n, s1), atan(n, s2)

    str / 4π * ((s2 * θ2 - s1 * θ1) + n * log(r1 / r2))
end

# two_piece_source(λ1, λ2, r1, r2, r) = (r - r1) / (r2 - r1) * λ1 + (r2 - r) / (r2 - r1) * λ2

# function source_stream_plus(str, x, n, x_p, x_q)
#   s_p, s_q, r_p, r_q, β_p, β_q = transform(x, n, x_p, x_q)
#   str / 4π * (x_p * β_p - x_q * β_q + n * log(r_p / r_q))
# end

# function source_stream_minus(str, x, n, x_p, x_q)
#   s_p, s_q, r_p, r_q, β_p, β_q = transform(x, n, x_p, x_q)
#   λ_plus = source_stream_1(1., x, n, x_p, x_q)

#   str / (4π * (x_p - x_q)) * ((x_p + x_q) * λ_plus + r_p^2 * β_p - r_q^2 * β_q + n * (x_q - x_p))
# end