"""
An abstract type to define custom spacing distributions.
"""
abstract type AbstractSpacing end

struct Cosine  <: AbstractSpacing end

struct Sine    <: AbstractSpacing 
    dir :: Bool
end

Sine() = Sine(0)

struct Uniform <: AbstractSpacing end

spacing(n, :: Uniform) = uniform_spacing(0., 1., n + 1)

spacing(n, :: Cosine) = cosine_spacing(0.5, 1., n + 1)

spacing(n, s :: Sine) = sine_spacing(0., 1., (n + 1) * ifelse(s.dir, -1, 1))

spacing(n, distribution :: AbstractSpacing) = spacing(n, distribution)

Base.reverse(q :: AbstractSpacing) = q
Base.Broadcast.broadcastable(q::AbstractSpacing) = Ref(q)

chop_sections(set1, set2, n :: Integer, distribution) = @views [ weighted_vector.(set1, set2, μ) for μ ∈ spacing(n, distribution) ][1:end-1]

chop_coordinates(coords, ns :: Vector{T}, spacings :: Vector{N}) where T <: Integer where N <: AbstractSpacing =
    @views [ mapreduce((c1, c2, n, space) -> chop_sections(c1, c2, n, space), vcat, coords[1:end-1], coords[2:end], ns, spacings); [ coords[end] ] ]

chop_coordinates(coords, n :: Integer, space :: AbstractSpacing) =
    @views [ mapreduce((c1, c2) -> chop_sections(c1, c2, n, space), vcat, coords[1:end-1], coords[2:end]); [ coords[end] ] ]

chop_coordinates(coords, n :: Vector{<:Integer}, space :: AbstractSpacing) = chop_coordinates(coords, n, fill(space, length(n)))

chop_spans(xyzs, span_num, spacing) = combinedimsview(map(xyz -> chop_coordinates(xyz, span_num, spacing), eachrow(xyzs)), (1))

chop_chords(xyzs, div, spacing) = combinedimsview(map(xyz -> chop_coordinates(xyz, div, spacing), eachcol(xyzs)))

function chop_wing(xyz_coords, span_nums, chord_num; span_spacing, chord_spacing)

    # Chordwise chop
    xyz_coords = chop_chords(xyz_coords, chord_num, chord_spacing)
    
    # Spanwise chop
    xyz_coords = chop_spans(xyz_coords, span_nums, span_spacing)

    return xyz_coords
end
# chop_coordinates

# Maybe switch to tensors?
transform_coordinates(xyz, twist, section) = eachrow(xyz * RotY(-twist)') .+ Ref(section)

function chop_spanwise_sections(scaled_foils, twisties, leading_xyz, span_nums, spacings, symmetry = false, flip = false)

    if symmetry
        scaled_foils = [ scaled_foils[end:-1:2]; scaled_foils ]
        twisties     = [ twisties[end:-1:2]; twisties ]
    elseif flip
        scaled_foils = reverse(scaled_foils)
        twisties     = reverse(twisties)
        spacings     = reverse(spacings)
    end

    # Rotate and translate coordinates
    foil_coords = combinedimsview(transform_coordinates.(scaled_foils, twisties, leading_xyz))

    # Chop up spanwise sections
    xyz_span = chop_spans(foil_coords, span_nums, spacings)

    return xyz_span
end