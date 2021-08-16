function chop_sections(set1, set2, n :: Integer, spacing = "cosine"; flip = false)
    if lowercase(spacing) == "uniform"
        space = uniform_spacing(0., 1., n + 1)
    elseif lowercase(spacing) == "sine"
        space = sine_spacing(0., 1., (n + 1) * ifelse(flip, -1, 1))
    else
        space = cosine_spacing(0.5, 1., n + 1)
    end

    @views [ weighted_vector.(set1, set2, μ) for μ ∈ space ][1:end-1]
end

chop_coordinates(coords, n, spacing = "cosine", flip = false) = @views [ reduce(vcat, chop_sections.(coords[1:end-1], coords[2:end], n, spacing; flip = flip)); [ coords[end] ] ]

chop_spans(xyzs, div, spacing = "cosine", flip = false) = permutedims(reduce(hcat, chop_coordinates(xyz, div, spacing, flip) for xyz in eachrow(xyzs)))

chop_chords(xyzs, div, spacing = "cosine") = reduce(hcat, chop_coordinates(xyz, div, spacing) for xyz in eachcol(xyzs))

chop_wing(xyzs, span_num, chord_num; span_spacing = "cosine", chord_spacing = "cosine", flip = false) = chop_chords(chop_spans(xyzs, span_num, span_spacing, flip), chord_num, chord_spacing)

transform_coordinates(xyz, twist, section) = eachrow(xyz * RotY(-twist)') .+ Ref(section)

function chop_spanwise_sections(scaled_foils, twisties, leading_xyz, span_num, spacings, flip = false)
    # Reverse direction if left side
    if flip
        scaled_foils = reverse(scaled_foils)
        twisties     = reverse(twisties)
    end

    # Rotate and translate coordinates
    foil_coords = reduce(hcat, transform_coordinates.(scaled_foils, twisties, leading_xyz))

    # Chop up spanwise sections
    chop_spans(foil_coords, span_num, spacings, flip)
end