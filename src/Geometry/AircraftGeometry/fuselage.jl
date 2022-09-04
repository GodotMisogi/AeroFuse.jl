struct Fuselage{T <: Real} <: AbstractAircraft
    length   :: T
    weights  :: Vector{T}
    radii    :: Vector{T}
    position :: SVector{3,T}
end

function Fuselage(L, weights, radii, position = zeros(3))
    T = promote_type(eltype(L), eltype(weights), eltype(radii), eltype(position))
    @assert all(x -> 0 <= x <= 1., weights) "The distribution of weights must be lie in [0,1] ⊂ ℝ."
    perms = sortperm(weights); sort!(weights)
    Fuselage{T}(L, weights, radii[perms], position)
end

projected_area(fuse :: Fuselage) = forward_sum(fuse.radii) / 2 .* fuse.weights
Base.length(fuse :: Fuselage) = fuse.length

truncated_cone_curved_area(r, R, H) = (R + r) * √(H^2 + (R - r)^2) * π

truncated_cone_volume(r, R, H) = H * π/3 * (R^2 + R * r + r^2) 

wetted_area(fuse :: Fuselage) = sum(x -> truncated_cone_curved_area(x...), zip(fuse.radii[1:end-1], fuse.radii[2:end], fuse.weights * fuse.length))

volume(fuse :: Fuselage) = sum(x -> truncated_cone_volume(x...), zip(fuse.radii[1:end-1], fuse.radii[2:end], fuse.weights * fuse.length))

function coordinates(fuse :: Fuselage, n)
    ws_rads = cosine_interpolation(fuse, n)

    [ ws_rads[:,1] .* fuse.length ws_rads[:,2] ]
end

function cosine_interpolation(fuse :: Fuselage, n)
    xs = fuse.weights
    ys = fuse.radii

    x_min, x_max = extrema(xs)
    x_circ = cosine_spacing((x_min + x_max) / 2, x_max - x_min, n)

    y_u = @views LinearInterpolation(xs, ys).(x_circ)

    return [ x_circ y_u ]
end

struct NewFuselage{T <: Real, N <: AbstractAffineMap, F <: Function} <: AbstractAffineMap
    radius :: T
    ell_a  :: T
    ell_b  :: T
    x_nose :: T
    x_end  :: T
    x_blend1 :: T
    x_blend2 :: T
    affine :: N
    droop_LE :: F
    droop_TE :: F
    power_coeffs :: Vector{T} 
end

hyperellipse(ξ, R_fuse, a) = R_fuse * (1 - ξ^a)^(1/a) 

function shape(fuse :: NewFuselage)
    # Linear spacing
    x1s = 0:0.01:1

    # Radii
    R_LE = hyperellipse.(x1s, fuse.radius, ell_a)

    # Longitudinal
    R_yp_LE = R_LE
    R_ym_LE = -R_yp_LE

    R_yp_TE = hyperellipse.(x1s, fuse.radius, ell_a)
    R_ym_TE = -R_yp_TE

    # Centers
    c_LE = fuse.droop_LE.(reverse(x1s))
    c_TE = fuse.droop_TE.(x1s)

    
    

end