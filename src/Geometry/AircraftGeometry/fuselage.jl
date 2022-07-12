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

truncated_cone_curved_area(r, R, H) = (R + r)√(H^2 + r^2)π

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