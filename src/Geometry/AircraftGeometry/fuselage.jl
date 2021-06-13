struct Fuselage{T <: Real} <: Aircraft
    lengths :: Vector{T}
    radii   :: Vector{T}
end

Fuselage(lengths :: AbstractVector{T}, radii :: AbstractVector{T}) where T <: Real = Fuselage{T}(lengths, radii)

projected_area(fuse :: Fuselage) = fwdsum(fuse.radii) / 2 .* fuse.lengths
Base.length(fuse :: Fuselage) = sum(fuse.lengths)

function cosine_spacing(fuse :: Fuselage, n)
    xs = [0.; cumsum(fuse.lengths) ]
    ys = fuse.radii
    cosine_interp([ xs ys ], n)
end