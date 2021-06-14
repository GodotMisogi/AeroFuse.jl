
"""
    wing_coords(wing :: HalfWing, n_s :: Integer, n_c :: Integer, flip = false)

Compute the coordinates of a `HalfWing` consisting of `Foil`s and relevant geometric quantities, given numbers of spanwise ``n_s`` and chordwise ``n_c`` panels, with an option to flip the signs of the ``y``-coordinates.
"""
function wing_coords(wing :: HalfWing, span_num :: Vector{<: Integer}, chord_num :: Integer, flip = false)
    leading_xyz  = leading_edge(wing, flip)
    temp_foils 	 = @. wing.chords * (extend_yz ∘ cosine_foil)(wing.foils, chord_num)

    scaled_foils = ifelse(flip, temp_foils[end:-1:1], temp_foils)
    twists  	 = ifelse(flip, wing.twists[end:-1:1], wing.twists)

    foil_coords  = [ coords * RotY(-twist)' .+ section' for (coords, twist, section) ∈ zip(scaled_foils, twists, leading_xyz) ]

    coords_chopper(foil_coords, span_num)
end

"""
    camber_coords(wing :: HalfWing, n_s :: Integer, n_c :: Integer, flip = false)

Compute the coordinates of a `HalfWing` consisting of camber distributions of `Foil`s and relevant geometric quantities, given numbers of spanwise ``n_s`` and chordwise ``n_c`` panels, with an option to flip the signs of the ``y``-coordinates.
"""
function camber_coords(wing :: HalfWing, span_num :: Vector{<: Integer}, chord_num :: Integer; spacings = "cosine", flip = false)
    leading_xyz  = leading_edge(wing, flip)

    scaled_foils = @. wing.chords * (camber_coordinates ∘ camber_thickness)(wing.foils, chord_num)

    scaled_foils = ifelse(flip, scaled_foils[end:-1:1], scaled_foils)
    twists 		 = ifelse(flip, wing.twists[end:-1:1], wing.twists)

    foil_coords  = [ coords * RotY(-twist)' .+ section' for (coords, twist, section) ∈ zip(scaled_foils, twists, leading_xyz) ]

    coords_chopper(foil_coords, span_num, spacings, flip)
end

"""
    mesh_horseshoes(obj :: HalfWing, n_s :: Integer, n_c :: Integer; flip = false)

Mesh a `HalfWing` into panels of ``n_s`` spanwise divisions per section and ``n_c`` chordwise divisions meant for lifting-line/vortex lattice analyses using horseshoe elements.
"""
mesh_horseshoes(obj :: HalfWing, span_num :: Vector{<: Integer}, chord_num :: Integer, spacings = "cosine"; flip = false) = make_panels(wing_chopper(wing_bounds(obj, flip)..., span_num, chord_num, span_spacing = spacings, flip = flip))

"""
    mesh_wing(wing :: HalfWing, n_s :: Integer, n_c :: Integer; flip = false)

Mesh a `HalfWing` into panels of ``n_s`` spanwise divisions per section and ``n_c`` chordwise divisions meant for 3D analyses using doublet-source elements or equivalent formulations. TODO: Tip meshing.
"""
mesh_wing(wing :: HalfWing, span_num :: Vector{<: Integer}, chord_num :: Integer; flip = false) = (make_panels ∘ wing_coords)(wing, span_num, chord_num, flip)

"""
    mesh_cambers(wing :: HalfWing, n_s :: Integer, n_c :: Integer; flip = false)

Mesh the camber distribution of a `HalfWing` into panels of ``n_s`` spanwise divisions per section and ``n_c`` chordwise divisions meant for vortex lattice analyses.
"""
mesh_cambers(wing :: HalfWing, span_num :: Vector{<: Integer}, chord_num :: Integer, spacings = "cosine"; flip = false) = make_panels(camber_coords(wing, span_num, chord_num; spacings = spacings, flip = flip))

leading_chopper(obj :: HalfWing, span_num; flip = false)  = coords_chopper(leading_edge(obj, flip), span_num)
trailing_chopper(obj :: HalfWing, span_num; flip = false) = coords_chopper(leading_edge(obj, flip), span_num)

leading_chopper(obj :: Wing, span_num :: Integer)  = span_chopper(wing_bounds(obj)..., span_num)[1]
trailing_chopper(obj :: Wing, span_num :: Integer) = span_chopper(wing_bounds(obj)..., span_num)[2]

"""
    mesh_horseshoes(wing :: Wing, n_s :: Integer, n_c :: Integer)

Mesh a `Wing` into panels of ``n_s`` spanwise divisions per section and ``n_c`` chordwise divisions meant for lifting-line/vortex lattice analyses using horseshoe elements.
"""
function mesh_horseshoes(wing :: Wing, span_num, chord_num, spacings = "sine")
    if wing.left === wing.right
        # Need to consider sine distribution at symmetry, and cosine thereafter
        left_panels  = mesh_horseshoes(wing.left, span_num[end:-1:1], chord_num, spacings, flip = true)
        right_panels = mesh_horseshoes(wing.right, span_num, chord_num, spacings)
    else
        left_panels  = mesh_horseshoes(wing.left, span_num[end:-1:1], chord_num, flip = true)
        right_panels = mesh_horseshoes(wing.right, span_num, chord_num)
    end

    [ left_panels right_panels ] 
end

"""
    mesh_wing(wing :: Wing, n_s :: Integer, n_c :: Integer)

Mesh a `Wing` into panels of ``n_s`` spanwise divisions per section and ``n_c`` chordwise divisions meant for 3D analyses using doublet-source elements or equivalent formulations. TODO: Tip meshing.
"""
function mesh_wing(wing :: Wing, span_num, chord_num)
    left_panels  = mesh_wing(wing.left, span_num[end:-1:1], chord_num, flip = true)
    right_panels = mesh_wing(wing.right, span_num, chord_num)

    [ left_panels right_panels ] 
end

"""
    mesh_cambers(wing :: Wing, n_s :: Integer, n_c :: Integer)

Mesh the camber distribution of a `Wing` into panels of ``n_s`` spanwise divisions per section and ``n_c`` chordwise divisions meant for vortex lattice analyses.
"""
function mesh_cambers(wing :: Wing, span_num, chord_num, spacings = "sine")
    if wing.left === wing.right
        # Need to consider sine distribution at symmetry, and cosine thereafter
        left_panels  = mesh_cambers(wing.left, span_num[end:-1:1], chord_num, spacings, flip = true)
        right_panels = mesh_cambers(wing.right, span_num, chord_num, spacings)
    else
        left_panels  = mesh_cambers(wing.left, span_num[end:-1:1], chord_num, flip = true)
        right_panels = mesh_cambers(wing.right, span_num, chord_num)
    end

    [ left_panels right_panels ] 
end

vlmesh_wing(wing :: Union{HalfWing, Wing}, span_num, chord_num, spacings = "sine") = mesh_horseshoes(wing, span_num, chord_num, spacings), mesh_cambers(wing, span_num, chord_num, spacings)

function paneller(wing :: Union{Wing, HalfWing}, span_num :: Union{Integer, Vector{<: Integer}}, chord_num :: Integer; rotation = one(RotMatrix{3, Float64}), translation = zeros(3), spacings = "sine")
    horseshoes, cambers = vlmesh_wing(wing, span_num, chord_num, spacings) 
    horseshoe_panels    = [ transform(panel, rotation, translation)               for panel in horseshoes ]
    panel_normals       = [ panel_normal(transform(panel, rotation, translation)) for panel in cambers 	  ]

    horseshoe_panels, panel_normals
end

panel_wing(comp :: Union{Wing, HalfWing}, span_panels :: Union{Integer, Vector{<: Integer}}, chord_panels :: Integer; position = zeros(3), angle = 0., axis = [1., 0., 0.]) = paneller(comp, span_panels, chord_panels, rotation = AngleAxis{Float64}(angle, axis...), translation = position)

number_of_spanwise_panels(wing :: HalfWing, span_num :: Integer) = ceil.(Int, span_num .* wing.spans / span(wing))
number_of_spanwise_panels(wing :: Wing, span_num :: Integer) = number_of_spanwise_panels(wing.right, ceil(Int, span_num / 2))

function number_of_spanwise_panels(wing :: HalfWing, span_num :: Vector{<: Integer}) 
    @assert length(wing.spans) > 1 "Provide an integer number of spanwise panels for 1 wing section."
    span_num
end

function number_of_spanwise_panels(wing :: Wing, span_num :: Vector{<: Integer}) 
    @assert length(wing.right.spans) > 1 "Provide an integer number of spanwise panels for 1 wing section."
    ceil.(Int, span_num ./ 2)
end