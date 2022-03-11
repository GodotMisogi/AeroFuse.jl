## Coordinates
#==========================================================================================#

wing_bounds(lead, trail) = permutedims([ lead trail ])

chop_leading_edge(obj  :: HalfWing, span_num; y_flip = false) = chop_coordinates(leading_edge(obj, y_flip),  span_num)
chop_trailing_edge(obj :: HalfWing, span_num; y_flip = false) = chop_coordinates(trailing_edge(obj, y_flip), span_num)

chop_leading_edge(obj :: Wing, span_num :: Integer)  = chop_coordinates([ leading_edge(left(obj), true)[1:end-1]; leading_edge(right(obj)) ], span_num)
chop_trailing_edge(obj :: Wing, span_num :: Integer) = chop_coordinates([ trailing_edge(left(obj), true)[1:end-1]; trailing_edge(right(obj)) ], span_num)

coordinates(wing :: HalfWing, y_flip = false) = let (lead, trail) = wing_bounds(wing, y_flip); affine_transformation(wing).(wing_bounds(lead, trail)) end
coordinates(wing :: Wing) = let (lead, trail) = wing_bounds(wing); affine_transformation(wing).(wing_bounds(lead, trail)) end

chord_coordinates(wing :: HalfWing, span_num :: Vector{<: Integer}, chord_num :: Integer; spacings = symmetric_spacing(wing), flip = false) = chop_wing(coordinates(wing, flip), span_num, chord_num; span_spacing = spacings, flip = flip)

function camber_coordinates(wing :: HalfWing, span_num :: Vector{<: Integer}, chord_num :: Integer; spacings = symmetric_spacing(wing), flip = false)
    leading_xyz  = leading_edge(wing, flip)
    scaled_foils = @. wing.chords * (camber_coordinates ∘ camber_thickness)(wing.foils, chord_num)
    affine_transformation(wing).(chop_spanwise_sections(scaled_foils, twists(wing), leading_xyz, span_num, spacings, flip))
end

function surface_coordinates(wing :: HalfWing, span_num :: Vector{<: Integer}, chord_num :: Integer; spacings = symmetric_spacing(wing), flip = false)
    leading_xyz  = leading_edge(wing, flip)
    scaled_foils = @. wing.chords * (extend_yz ∘ coordinates ∘ cosine_interpolation)(wing.foils, chord_num)
    affine_transformation(wing).(chop_spanwise_sections(scaled_foils, twists(wing), leading_xyz, span_num, spacings, flip))
end

function number_of_spanwise_panels(wing :: HalfWing, span_num :: Integer) 
    # Compute contribution of each section to total span length
    weights = spans(wing) / span(wing)

    weights[findall(<(0.2), weights)] .*= 3

    # Heuristic (aka hAx0rZ) check to ensure small sections also get some panel love
    # weights = ifelse(any(<(0.2), weights), fill(1. / length(spans(wing)), length(spans(wing))), weights)

    # Generate spanwise panel distribution
    ceil.(Int, span_num .* weights)
end

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

mesh_chords(wing :: HalfWing, span_num :: Vector{<: Integer}, chord_num :: Integer; spacings = symmetric_spacing(wing), flip = false) = make_panels(chord_coordinates(wing, span_num, chord_num; spacings = spacings, flip = flip))

function mesh_chords(wing :: Wing, span_num, chord_num; spacings = symmetric_spacing(wing))
    left_panels  = mesh_chords(left(wing), reverse(span_num), chord_num; spacings = reverse(spacings), flip = true)
    right_panels = mesh_chords(right(wing), span_num, chord_num; spacings = spacings)

    [ left_panels right_panels ]
end

mesh_wing(wing :: HalfWing, span_num :: Vector{<: Integer}, chord_num :: Integer; spacings = symmetric_spacing(wing), flip = false) = make_panels(surface_coordinates(wing, span_num, chord_num, spacings = spacings, flip = flip))

function mesh_wing(wing :: Wing, span_num, chord_num; spacings = symmetric_spacing(wing))
    left_panels  = mesh_wing(left(wing), reverse(span_num), chord_num; spacings = reverse(spacings), flip = true)
    right_panels = mesh_wing(right(wing), span_num, chord_num; spacings = spacings)

    [ left_panels right_panels ]
end

mesh_cambers(wing :: HalfWing, span_num :: Vector{<: Integer}, chord_num :: Integer; spacings = symmetric_spacing(wing), flip = false) = make_panels(camber_coordinates(wing, span_num, chord_num; spacings = spacings, flip = flip))

function mesh_cambers(wing :: Wing, span_num, chord_num; spacings = symmetric_spacing(wing))
    left_panels  = mesh_cambers(left(wing), reverse(span_num), chord_num; spacings = reverse(spacings), flip = true)
    right_panels = mesh_cambers(right(wing), span_num, chord_num; spacings = spacings)

    [ left_panels right_panels ]
end

function make_panels(wing :: AbstractWing, span_num :: Vector{<: Integer}, chord_num :: Integer; spacings = symmetric_spacing(wing))
    horseshoe_panels = mesh_chords(wing, span_num, chord_num; spacings = spacings)
    camber_panels    = mesh_cambers(wing, span_num, chord_num; spacings = spacings)
    horseshoe_panels, normal_vector.(camber_panels)
end

make_panels(wing :: AbstractWing, span_num :: Integer, chord_num :: Integer; spacings = symmetric_spacing(wing)) = make_panels(wing, [span_num], chord_num; spacings = spacings)

panel_wing(comp :: AbstractWing, span_panels :: Union{Integer, Vector{<: Integer}}, chord_panels :: Integer; spacing = symmetric_spacing(comp)) = make_panels(comp, span_panels, chord_panels, spacings = spacing)

## Meshing type for convenience
#==========================================================================================#

struct WingMesh{M <: AbstractWing, N <: Integer, P, Q, T} <: AbstractWing
    surface       :: M
    num_span      :: Vector{N}
    num_chord     :: N
    chord_spacing :: P
    span_spacing  :: Q
    chord_mesh    :: Matrix{T}
    camber_mesh   :: Matrix{T}
end

"""
    WingMesh(
             surface :: AbstractWing, 
             n_span :: Vector{Integer}, n_chord :: Integer;
             span_spacing :: AbstractSpacing = symmetric_spacing(surface)
            )

Define a container to generate meshes and panels for a given `AbstractWing` with a specified distribution of number of spanwise panels, and a number of chordwise panels.

Optionally a combination of `AbstractSpacing` types (`Sine(), Cosine(), Uniform()`) can be provided to the **named argument** `span_spacing`, either as a singleton or as a vector with length equal to the number of spanwise sections. By default, the combination is `[Sine(), Cosine(), ..., Cosine()]`.
"""
function WingMesh(surface :: M, n_span :: AbstractVector{N}, n_chord :: N; chord_spacing :: P = Cosine(), span_spacing :: Q = symmetric_spacing(surface)) where {M <: AbstractWing, N <: Integer, P <: AbstractSpacing, Q <: Union{AbstractSpacing, Vector{<:AbstractSpacing}}}
    check_definition(surface, n_span)
    chord_mesh  = chord_coordinates(surface, n_span, n_chord; spacings = span_spacing)
    camber_mesh = camber_coordinates(surface, n_span, n_chord; spacings = span_spacing)
    T = promote_type(eltype(chord_mesh), eltype(camber_mesh))
    WingMesh{M,N,P,Q,T}(surface, n_span, n_chord, chord_spacing, span_spacing, chord_mesh, camber_mesh)
end

WingMesh(surface, n_span :: Integer, n_chord :: Integer; chord_spacing = Cosine(), span_spacing = symmetric_spacing(surface)) = WingMesh(surface, number_of_spanwise_panels(surface, n_span), n_chord; chord_spacing = chord_spacing, span_spacing = span_spacing)

check_definition(surf :: HalfWing, n_span) = @assert length(n_span) == length(surf.spans) "The spanwise number vector's length must be the same as the number of sections of the surface."
check_definition(surf :: Wing, n_span) = @assert length(n_span) == length(surf.right.spans) == length(surf.left.spans) "The spanwise number vector's length must be the same as the number of sections of the surface."

##
chord_coordinates(wing :: WingMesh, n_span = wing.num_span, n_chord = wing.num_chord) = chord_coordinates(wing.surface, n_span, n_chord)
camber_coordinates(wing :: WingMesh, n_span = wing.num_span, n_chord = wing.num_chord) = camber_coordinates(wing.surface, n_span, n_chord)
surface_coordinates(wing :: WingMesh, n_span = wing.num_span, n_chord = wing.num_chord) = surface_coordinates(wing.surface, n_span, n_chord)

surface_panels(wing :: WingMesh, n_span = wing.num_span, n_chord = wing.num_chord)  = (make_panels ∘ surface_coordinates)(wing, n_span, n_chord)

"""
    chord_panels(wing_mesh :: WingMesh)

Generate the chord panel distribution from a `WingMesh`.
"""
chord_panels(wing :: WingMesh) = make_panels(wing.chord_mesh)

"""
    camber_panels(wing_mesh :: WingMesh)

Generate the camber panel distribution from a `WingMesh`.
"""
camber_panels(wing :: WingMesh) = make_panels(wing.camber_mesh)

function Base.show(io :: IO, mesh :: WingMesh)
    n_c, n_s = size(mesh.chord_mesh) .- 1
    println(io, "WingMesh —")
    println(io, "Spanwise panels: ", n_s)
    println(io, "Chordwise panels: ", n_c)
    println(io, "Spanwise spacing: ", mesh.span_spacing)
    println(io, "Chordwise spacing: ", mesh.chord_spacing)

    nothing
end