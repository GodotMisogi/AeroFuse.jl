"""
An abstract type to define custom spacing distributions.
"""
abstract type AbstractSpacing end

struct Cosine <: AbstractSpacing end

struct Sine <: AbstractSpacing 
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

# Rotations.jl uses StaticArrays which doesn't appear to be compatible with the pullbacks in ReverseDiff.
y_rotation(θ) = ((sθ, cθ) = sincos(θ);
    [ cθ 0 -sθ;
      0  1   0;
      sθ 0  cθ ]
)

# TODO: Maybe switch to multidimensional arrays
transform_coordinates(xyz, rot, section) = eachrow(xyz * rot) .+ Ref(section)

@views function chop_spanwise_sections(scaled_foils, twisties, leading_xyz, span_nums, spacings, symmetry = false, flip = false)

    if symmetry
        scaled_foils = [ scaled_foils[end:-1:2]; scaled_foils ]
        twisties     = [ twisties[end:-1:2]; twisties ]
    elseif flip
        scaled_foils = reverse(scaled_foils)
        twisties     = reverse(twisties)
        spacings     = reverse(spacings)
    end

    # Rotate and translate coordinates
    foil_coords = @. transform_coordinates(scaled_foils, permutedims(y_rotation(-twisties)), leading_xyz)

    # Chop up spanwise sections
    foil_coords = combinedimsview(foil_coords)
    xyz_span = chop_spans(foil_coords, span_nums, spacings)

    return xyz_span
end

# @views function chop_spanwise_sections(scaled_foils, twisties, leading_xyz, span_nums, spacings, symmetry = false, flip = false)

#     # if symmetry
#     #     scaled_foils = [ scaled_foils[end:-1:2]; scaled_foils ]
#     #     twisties     = [ twisties[end:-1:2]; twisties ]
#     # elseif flip
#     #     scaled_foils = reverse(scaled_foils)
#     #     twisties     = reverse(twisties)
#     #     spacings     = reverse(spacings)
#     # end

#     # Rotate and translate coordinates
#     # foil_coords = @. transform_coordinates(scaled_foils, permutedims(y_rotation(-twisties)), leading_xyz)

#     # Chop up spanwise sections
#     # foil_coords = combinedimsview(foil_coords)
#     # xyz_span = chop_spans(foil_coords, span_nums, spacings)

#     # @show length(scaled_foils), length(twisties), length(leading_xyz), span_nums, spacings

#     xyz_span = reduce(hcat, map(enumerate(zip(scaled_foils[1:end-1], scaled_foils[2:end], twisties[1:end-1], twisties[2:end], leading_xyz[1:end-1], leading_xyz[2:end]))) do (n, (f1, f2, t1, t2, l1, l2))
#         # Rotate and translate airfoils of current spanwise section
#         bf1 = transform_coordinates(f1, permutedims(y_rotation(-t1)), l1)
#         bf2 = transform_coordinates(f2, permutedims(y_rotation(-t2)), l2)

#         # Chop spanwise section
#         chop_spans(combinedimsview([bf1, bf2]), span_nums[n], spacings[n])
#     end)

#     return [ reverse(xyz_span, dims=2) xyz_span ]
# end
