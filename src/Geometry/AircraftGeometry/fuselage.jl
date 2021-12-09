struct Fuselage{T <: Real} <: AbstractAircraft
    length   :: T
    weights  :: Vector{T}
    radii    :: Vector{T}
    position :: SVector{3,T}
end

function Fuselage(L :: T, weights :: AbstractVector{T}, radii :: AbstractVector{T}, position = zeros(3)) where T <: Real
    @assert all(x -> 0 <= x <= 1., weights) "The distribution of weights must be lie in [0,1] ⊂ ℝ."
    perms = sortperm(weights); sort!(weights)
    Fuselage{T}(L, weights, radii[perms], position)
end

projected_area(fuse :: Fuselage) = fwdsum(fuse.radii) / 2 .* fuse.weights
Base.length(fuse :: Fuselage) = fuse.length

function coordinates(fuse :: Fuselage, n)
    ws_rads = cosine_spacing(fuse, n)

    [ ws_rads[:,1] .* fuse.length ws_rads[:,2] ]
end

function cosine_spacing(fuse :: Fuselage, n)
    xs = fuse.weights
    ys = fuse.radii
    cosine_interp([ xs ys ], n)
end