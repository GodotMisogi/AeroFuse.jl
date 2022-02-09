## Coordinates
#==========================================================================================#

wing_bounds(lead, trail) = permutedims([ lead trail ])

chop_leading_edge(obj  :: HalfWing, span_num; y_flip = false) = chop_coordinates(leading_edge(obj, y_flip),  span_num)
chop_trailing_edge(obj :: HalfWing, span_num; y_flip = false) = chop_coordinates(trailing_edge(obj, y_flip), span_num)

chop_leading_edge(obj :: Wing, span_num :: Integer)  = chop_coordinates([ leading_edge(left(obj), true)[1:end-1]; leading_edge(right(obj)) ], span_num)
chop_trailing_edge(obj :: Wing, span_num :: Integer) = chop_coordinates([ trailing_edge(left(obj), true)[1:end-1]; trailing_edge(right(obj)) ], span_num)

coordinates(wing :: HalfWing, y_flip = false) = let (lead, trail) = wing_bounds(wing, y_flip); affine_transformation(wing).(wing_bounds(lead, trail)) end
coordinates(wing :: Wing) = let (lead, trail) = wing_bounds(wing); affine_transformation(wing).(wing_bounds(lead, trail)) end

"""
    coordinates(wing :: AbstractWing, n_s :: Integer, n_c :: Integer, flip = false)

Compute the planform coordinates of a `HalfWing` given numbers of spanwise ``n_s`` and chordwise ``n_c`` panels, with an option to flip the signs of the ``y``-coordinates.
"""
coordinates(wing :: AbstractWing, span_num, chord_num; span_spacing = Cosine(), chord_spacing = Cosine()) = chop_wing(coordinates(wing), span_num, chord_num; span_spacing = span_spacing, chord_spacing = chord_spacing)

"""
    chord_coordinates(wing :: AbstractWing, n_s :: Integer, n_c :: Integer, flip = false)

Compute the chord coordinates of a `HalfWing` given numbers of spanwise ``n_s`` and chordwise ``n_c`` panels, with an option to flip the signs of the ``y``-coordinates.
"""
chord_coordinates(wing :: HalfWing, span_num :: Vector{<: Integer}, chord_num :: Integer; spacings = symmetric_spacing(wing), flip = false) = chop_wing(coordinates(wing, flip), span_num, chord_num; span_spacing = spacings, flip = flip)

"""
    surface_coordinates(wing :: HalfWing, n_s :: Integer, n_c :: Integer, flip = false)

Compute the surface coordinates of a `HalfWing` consisting of `Foil`s and relevant geometric quantities, given numbers of spanwise ``n_s`` and chordwise ``n_c`` panels, with an option to flip the signs of the ``y``-coordinates.
"""
function surface_coordinates(wing :: HalfWing, span_num :: Vector{<: Integer}, chord_num :: Integer; spacings = symmetric_spacing(wing), flip = false)
    leading_xyz  = leading_edge(wing, flip)
    scaled_foils = @. wing.chords * (extend_yz ∘ coordinates ∘ cosine_spacing)(wing.foils, chord_num)
    affine_transformation(wing).(chop_spanwise_sections(scaled_foils, twists(wing), leading_xyz, span_num, spacings, flip))
end

"""
    camber_coordinates(wing :: HalfWing, n_s :: Integer, n_c :: Integer, flip = false)

Compute the camber coordinates of a `HalfWing` consisting of camber distributions of `Foil`s and relevant geometric quantities, given numbers of spanwise ``n_s`` and chordwise ``n_c`` panels, with an option to flip the signs of the ``y``-coordinates.
"""
function camber_coordinates(wing :: HalfWing, span_num :: Vector{<: Integer}, chord_num :: Integer; spacings = symmetric_spacing(wing), flip = false)
    leading_xyz  = leading_edge(wing, flip)
    scaled_foils = @. wing.chords * (camber_coordinates ∘ camber_thickness)(wing.foils, chord_num)
    affine_transformation(wing).(chop_spanwise_sections(scaled_foils, twists(wing), leading_xyz, span_num, spacings, flip))
end

number_of_spanwise_panels(wing :: HalfWing, span_num :: Integer) = ceil.(Int, span_num .* spans(wing) / span(wing))

function number_of_spanwise_panels(wing :: HalfWing, span_num :: Vector{<: Integer})
    @assert (length ∘ spans)(wing) > 1 "Provide a positive integer of spanwise panels for 1 wing section."
    span_num
end

# Spacing
symmetric_spacing(wing :: HalfWing) = [ Sine(); fill(Cosine(), (length ∘ spans)(wing) - 1) ]


## Wing variants
#==========================================================================================#

function chord_coordinates(wing :: Wing, span_num :: Vector{<: Integer}, chord_num :: Integer; spacings = symmetric_spacing(wing))
    left_coord  = chord_coordinates(left(wing), reverse(span_num), chord_num; spacings = reverse(spacings), flip = true)
    right_coord = chord_coordinates(right(wing), span_num, chord_num; spacings = spacings)

    [ left_coord[:,1:end-1] right_coord ]
end

function camber_coordinates(wing :: Wing, span_num :: Vector{<: Integer}, chord_num :: Integer; spacings = symmetric_spacing(wing))
    left_coord  = camber_coordinates(left(wing), reverse(span_num), chord_num; spacings = reverse(spacings), flip = true)
    right_coord = camber_coordinates(right(wing), span_num, chord_num; spacings = spacings)

    [ left_coord[:,1:end-1] right_coord ]
end

function surface_coordinates(wing :: Wing, span_num :: Vector{<: Integer}, chord_num :: Integer; spacings = symmetric_spacing(wing))
    left_coord  = surface_coordinates(left(wing), reverse(span_num), chord_num; spacings = reverse(spacings), flip = true)
    right_coord = surface_coordinates(right(wing), span_num, chord_num; spacings = spacings)

    [ left_coord[:,1:end-1] right_coord ]
end

number_of_spanwise_panels(wing :: Wing, span_num :: Integer) = number_of_spanwise_panels(right(wing), span_num ÷ 2)

symmetric_spacing(wing :: Wing)     = [ Sine(); fill(Cosine(), (length ∘ spans ∘ right)(wing) - 1) ]

function number_of_spanwise_panels(wing :: Wing, span_num :: Vector{<: Integer})
    @assert (length ∘ spans ∘ right)(wing) > 1 "Provide a positive integer of spanwise panels for 1 wing section."
    span_num .÷ 2
end

# Coordinates
chord_coordinates(wing :: Wing, span_num :: Integer, chord_num :: Integer; spacings = symmetric_spacing(wing)) = chord_coordinates(wing, number_of_spanwise_panels(wing, span_num), chord_num; spacings = spacings)

camber_coordinates(wing :: Wing, span_num :: Integer, chord_num :: Integer; spacings = symmetric_spacing(wing)) = camber_coordinates(wing, number_of_spanwise_panels(wing, span_num), chord_num; spacings = spacings)

## Panelling
#==========================================================================================#

"""
    mesh_chords(wing :: AbstractWing, n_s :: Vector{Integer}, n_c :: Integer; flip = false)

Mesh the span and chord distributions of an `AbstractWing` with ``n_s`` spanwise divisions per section and ``n_c`` chordwise divisions.
"""
mesh_chords(wing :: AbstractWing, span_num, chord_num; spacings = symmetric_spacing(wing)) = mesh_chords(wing, span_num, chord_num; spacings = spacings)

mesh_chords(wing :: HalfWing, span_num :: Vector{<: Integer}, chord_num :: Integer; spacings = symmetric_spacing(wing), flip = false) = make_panels(chord_coordinates(wing, span_num, chord_num; spacings = spacings, flip = flip))

function mesh_chords(wing :: Wing, span_num, chord_num; spacings = symmetric_spacing(wing))
    left_panels  = mesh_chords(left(wing), reverse(span_num), chord_num; spacings = reverse(spacings), flip = true)
    right_panels = mesh_chords(right(wing), span_num, chord_num; spacings = spacings)

    [ left_panels right_panels ]
end

"""
    mesh_wing(wing :: AbstractWing, n_s :: Vector{Integer}, n_c :: Integer; flip = false)

Mesh the span and airfoil coordinate distributions of an `AbstractWing` with ``n_s`` spanwise divisions per section and ``n_c`` chordwise divisions.
"""
mesh_wing(wing :: AbstractWing, span_num, chord_num; spacings = symmetric_spacing(wing)) = mesh_wing(wing, span_num, chord_num; spacings = spacings)

mesh_wing(wing :: HalfWing, span_num :: Vector{<: Integer}, chord_num :: Integer; spacings = symmetric_spacing(wing), flip = false) = make_panels(surface_coordinates(wing, span_num, chord_num, spacings = spacings, flip = flip))

function mesh_wing(wing :: Wing, span_num, chord_num; spacings = symmetric_spacing(wing))
    left_panels  = mesh_wing(left(wing), reverse(span_num), chord_num; spacings = reverse(spacings), flip = true)
    right_panels = mesh_wing(right(wing), span_num, chord_num; spacings = spacings)

    [ left_panels right_panels ]
end

"""
    mesh_cambers(wing :: AbstractWing, n_s :: Integer, n_c :: Integer; spacings = symmetric_spacing(wing))

Mesh the camber distribution of a `Wing` into panels of ``n_s`` spanwise divisions per section and ``n_c`` chordwise divisions with an `AbstractSpacing`` distribution.
"""
mesh_cambers(wing :: AbstractWing, span_num, chord_num; spacings = symmetric_spacing(wing)) = mesh_cambers(wing, span_num, chord_num; spacings = spacings)

mesh_cambers(wing :: HalfWing, span_num :: Vector{<: Integer}, chord_num :: Integer; spacings = symmetric_spacing(wing), flip = false) = make_panels(camber_coordinates(wing, span_num, chord_num; spacings = spacings, flip = flip))

function mesh_cambers(wing :: Wing, span_num, chord_num; spacings = symmetric_spacing(wing))
    left_panels  = mesh_cambers(left(wing), reverse(span_num), chord_num; spacings = reverse(spacings), flip = true)
    right_panels = mesh_cambers(right(wing), span_num, chord_num; spacings = spacings)

    [ left_panels right_panels ]
end

function make_panels(wing :: AbstractWing, span_num :: Vector{<: Integer}, chord_num :: Integer; spacings = symmetric_spacing(wing))
    horseshoe_panels = mesh_chords(wing, span_num, chord_num; spacings = spacings)
    camber_panels    = mesh_cambers(wing, span_num, chord_num; spacings = spacings)
    horseshoe_panels, panel_normal.(camber_panels)
end

make_panels(wing :: AbstractWing, span_num :: Integer, chord_num :: Integer; spacings = symmetric_spacing(wing)) = make_panels(wing, [span_num], chord_num; spacings = spacings)

panel_wing(comp :: AbstractWing, span_panels :: Union{Integer, Vector{<: Integer}}, chord_panels :: Integer; spacing = symmetric_spacing(comp)) = make_panels(comp, span_panels, chord_panels, spacings = spacing)

## Meshing type for convenience
#==========================================================================================#

mutable struct WingMesh{M <: AbstractWing, N <: Integer, P, Q, T} <: AbstractWing
    surf          :: M
    num_span      :: Vector{N}
    num_chord     :: N
    chord_spacing :: P
    span_spacing  :: Q
    chord_mesh    :: Matrix{T}
    camber_mesh   :: Matrix{T}
end

function WingMesh(surf :: M, n_span :: AbstractVector{N}, n_chord :: N; chord_spacing :: P = Cosine(), span_spacing :: Q = symmetric_spacing(surf)) where {M <: AbstractWing, N <: Integer, P <: AbstractSpacing, Q <: Union{AbstractSpacing, Vector{<:AbstractSpacing}}}
    check_definition(surf, n_span)
    chord_mesh  = chord_coordinates(surf, n_span, n_chord; spacings = span_spacing)
    camber_mesh = camber_coordinates(surf, n_span, n_chord; spacings = span_spacing)
    T = promote_type(eltype(chord_mesh), eltype(camber_mesh))
    WingMesh{M,N,P,Q,T}(surf, n_span, n_chord, chord_spacing, span_spacing, chord_mesh, camber_mesh)
end

check_definition(surf :: HalfWing, n_span) = @assert length(n_span) == length(surf.spans) "The spanwise number vector's length must be the same as the number of sections of the surface."
check_definition(surf :: Wing, n_span) = @assert length(n_span) == length(surf.right.spans) == length(surf.left.spans) "The spanwise number vector's length must be the same as the number of sections of the surface."

##
chord_coordinates(wing :: WingMesh, n_span = wing.num_span, n_chord = wing.num_chord) = chord_coordinates(wing.surf, n_span, n_chord)
camber_coordinates(wing :: WingMesh, n_span = wing.num_span, n_chord = wing.num_chord) = camber_coordinates(wing.surf, n_span, n_chord)
surface_coordinates(wing :: WingMesh, n_span = wing.num_span, n_chord = wing.num_chord) = surface_coordinates(wing.surf, n_span, n_chord)

surface_panels(wing :: WingMesh, n_span = wing.num_span, n_chord = wing.num_chord)  = (make_panels ∘ surface_coordinates)(wing, n_span, n_chord)

chord_panels(wing :: WingMesh) = make_panels(wing.chord_mesh)
camber_panels(wing :: WingMesh) = make_panels(wing.camber_mesh)
normal_vectors(wing :: WingMesh) = panel_normal.(camber_panels(wing))

function Base.show(io :: IO, mesh :: WingMesh)
    n_c, n_s = size(mesh.chord_mesh) .- 1
    println(io, "WingMesh —")
    println(io, "Spanwise panels: ", n_s)
    println(io, "Chordwise panels: ", n_c)
    println(io, "Spanwise spacing: ", mesh.span_spacing)
    println(io, "Chordwise spacing: ", mesh.chord_spacing)

    nothing
end