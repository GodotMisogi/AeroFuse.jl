struct Fuselage{T <: Real} <: AbstractFuselage
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

@views wetted_area(fuse :: Fuselage) = sum(x -> truncated_cone_curved_area(x...), zip(fuse.radii[1:end-1], fuse.radii[2:end], fuse.weights * fuse.length))

@views volume(fuse :: Fuselage) = sum(x -> truncated_cone_volume(x...), zip(fuse.radii[1:end-1], fuse.radii[2:end], fuse.weights * fuse.length))

@views function coordinates(fuse :: Fuselage, n)
    ws_rads = cosine_interpolation(fuse, n)

    [ ws_rads[:,1] .* fuse.length ws_rads[:,2] ]
end

function cosine_interpolation(fuse :: Fuselage, n)
    xs = fuse.weights
    ys = fuse.radii

    x_min, x_max = extrema(xs)
    x_circ = cosine_spacing((x_min + x_max) / 2, x_max - x_min, n)

    y_u = linear_interpolation(xs, ys).(x_circ)

    return [ x_circ y_u ]
end

struct HyperEllipseFuselage{T <: Real, N <: AbstractAffineMap} <: AbstractFuselage
    radius :: T
    length :: T
    x_a    :: T
    x_b    :: T
    xi_a  :: T
    xi_b  :: T
    d_nose :: T
    d_rear :: T
    affine :: N
    function HyperEllipseFuselage(R, L, x1, x2, ell1, ell2, d_nose, d_rear, affine)
        # Type promotion for autodiff support
        T = promote_type(eltype(R), eltype(L), eltype(x1), eltype(x2),  eltype(ell1), eltype(ell2), eltype(affine.linear), eltype(affine.translation), eltype(d_nose), eltype(d_rear))
        N = typeof(affine)

        @assert 0 < x1 < x2 < 1 "Ellipse blending points must lie between 0 and 1!"

        new{T,N}(R, L, x1, x2, ell1, ell2, d_nose, d_rear, affine)
    end
end

    
function HyperEllipseFuselage(;
        radius      = 1., 
        length      = 1., 
        x_a         = 0.2, 
        x_b         = 0.7,
        c_nose      = 1.6,
        c_rear  = 1.2,
        d_nose      = 0.,
        d_rear      = 0.,
        position    = zeros(3),
        angle       = 0.,
        axis        = [0., 1., 0.],
        affine      = AffineMap(AngleAxis(deg2rad(angle), axis...), position),
    )
   
    return HyperEllipseFuselage(radius, length, x_a, x_b, c_nose, c_rear, d_nose, d_rear, affine)
end

hyperellipse(ξ, a) = (1 - ξ^a)^(1/a) 

# Generate circles in the y-z plane with optionally shifted z-coordinates of the centers
function circle3D_yz(x, R, n, zs = zeros(n)) 
    ts = LinRange(0, 2π, n)
    @. [ fill(x, n) R * cos(ts) R * sin(ts) + zs ]
end

function undrooped_curve(fuse :: HyperEllipseFuselage, ts)
    @assert ts[1] >= 0
    @assert ts[end] <= 1

    x_a = fuse.x_a
    x_b = fuse.x_b
    R_f  = fuse.radius
    L_f  = fuse.length

    x_nose = fuse.affine.translation[1]
    x_end  = x_nose + fuse.length

    L_nose  = (x_a - 0) * L_f   # Nose length
    L_cabin = (x_b - x_a) * L_f # Cabin length
    L_rear  = (1 - x_b) * L_f   # Rear length

    # x-coordinates
    x_nose  = @. ts * L_nose
    x_cabin = @. ts * L_cabin + x_nose[end]
    x_rear  = @. ts * L_rear + x_cabin[end]

    y_nose  = @. R_f * hyperellipse(ts, fuse.xi_a)[end:-1:1]
    y_cabin = fill(R_f, length(x_cabin))
    y_rear  = @. R_f * hyperellipse(ts, fuse.xi_b)

    return x_nose, x_cabin, x_rear, y_nose, y_cabin, y_rear
end

function curve(fuse :: HyperEllipseFuselage, ts)
	x_nose, x_cabin, x_rear, z_nose, z_cabin, z_rear = undrooped_curve(fuse, ts)

	# Droop/rise
	z_nose += LinRange(fuse.d_nose, 0, length(x_nose))
	z_rear += LinRange(0, fuse.d_rear, length(x_rear))

	return x_nose, x_cabin, x_rear, z_nose, z_cabin, z_rear
end

function coordinates(fuse :: HyperEllipseFuselage, ts, n_circ = 20)
    x_nose, x_cabin, x_rear, z_nose, z_cabin, z_rear = undrooped_curve(fuse, ts)

	# Droop/rise
	dy_nose = LinRange(fuse.d_nose, 0, length(x_nose))
	dy_rear = LinRange(0, fuse.d_rear, length(x_rear))
	
	coords = reduce(vcat, [ 
		circle3D_yz.(x_nose, z_nose, n_circ, dy_nose); 
		circle3D_yz.(x_cabin, z_cabin, n_circ); 
		circle3D_yz.(x_rear, z_rear, n_circ, dy_rear)
	])

    return coords
end

function volume(fuse :: HyperEllipseFuselage)
    x_nose, x_cabin, x_rear, z_nose, z_cabin, z_rear = curve(fuse, ts)

    # Droop/rise
    z_nose += LinRange(fuse.d_nose, 0, length(x_nose))
    z_rear += LinRange(0, fuse.d_rear, length(x_rear))

    xs = [x_nose; x_cabin; x_rear]
    zs = [z_nose; z_cabin; z_rear]

    @views sum(x -> truncated_cone_volume(x...), zip(zs[1:end-1], zs[2:end], diff(xs)))
end