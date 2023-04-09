## Coordinates
#==========================================================================================#

chop_leading_edge(obj :: Wing, span_num) = @views chop_coordinates(coordinates(obj)[1,:], span_num, Uniform())
chop_trailing_edge(obj :: Wing, span_num) = @views chop_coordinates(coordinates(obj)[end,:], span_num, Uniform())

function coordinates(wing :: Wing, affine = true)
    bounds = permutedims(wing_bounds(wing), (3,2,1))

    # Symmetry
    if wing.symmetry
        sym_bounds = bounds[:,:,end:-1:2]
        sym_bounds[:,2,:] .*= -1

        bounds = cat(sym_bounds, bounds, dims = Val(3))
    # Reflection
    elseif wing.flip
        sym_bounds = bounds[:,:,end:-1:1]
        sym_bounds[:,2,:] .*= -1

        bounds = sym_bounds
    end

    # Reshape into matrix of vectors
    vec_bounds = splitdimsview(bounds,(1,3))

    # TODO: Weird hack
    if affine
        aff = wing.affine
    else
        aff = AffineMap(AngleAxis(0., 1, 0, 0.), zeros(3))
    end

    # @tullio t1[i,k,l] := aff.linear[i,j] * bounds[k,j,l]
    # @tullio verbose=true t2[i,k,l] := bounds[k,i,l] + aff.translation[i]
    # @tullio t3[i,k,l] := (aff.linear[i,j] * coo2[j,k,l]) + aff.translation[i]

    return map(aff, vec_bounds)
end

# coo2: [[x,y,z],[le,te],secs]

function chord_coordinates(wing :: Wing, span_num :: Vector{<: Integer}, chord_num :: Integer; span_spacing = symmetric_spacing(wing), chord_spacing = Cosine())
    coords = chop_wing(coordinates(wing, false), span_num, chord_num; span_spacing = span_spacing, chord_spacing = chord_spacing)

    return map(wing.affine, coords)
end

function camber_coordinates(wing :: Wing, span_num :: Vector{<: Integer}, chord_num :: Integer; span_spacing = symmetric_spacing(wing))
    # Get leading edge coordinates
    leading_xyz = @views coordinates(wing, false)[1,:]

    # Scale camber distribution of airfoil
    # scaled_foils = @. wing.chords * camber_coordinates(camber_thickness(wing.foils, chord_num))
    scaled_foils =  map((c, foil) -> c * camber_coordinates(camber_thickness(foil, chord_num)), wing.chords, wing.foils) 

    # # Discretize spanwise sections
    coords = chop_spanwise_sections(scaled_foils, -deg2rad.(wing.twists), leading_xyz, span_num, span_spacing, wing.symmetry, wing.flip)

    # # Transform
    return map(wing.affine, coords)
end

function surface_coordinates(wing :: Wing, span_num :: Vector{<: Integer}, chord_num :: Integer; span_spacing = symmetric_spacing(wing))
    # Get leading edge coordinates
    leading_xyz  = @views coordinates(wing, false)[1,:]

    # Scale airfoil coordinates
    scaled_foils = @. wing.chords * (extend_yz ∘ coordinates ∘ cosine_interpolation)(reflect.(wing.foils), chord_num)

    # Discretize spanwise sections
    coords = chop_spanwise_sections(scaled_foils, -deg2rad.(wing.twists), leading_xyz, span_num, span_spacing, wing.symmetry, wing.flip)

    # Transform
    return map(wing.affine, coords)
end

function number_of_spanwise_panels(wing :: Wing, span_num :: Integer) 
    # Compute contribution of each section to total span length
    weights = spans(wing) / span(wing)

    weights[findall(<(0.2), weights)] .*= 3

    # Heuristic (aka hAx0rZ) check to ensure small sections also get some panel love
    # weights = ifelse(any(<(0.2), weights), fill(1. / length(spans(wing)), length(spans(wing))), weights)

    # Generate spanwise panel distribution
    return ceil.(Int, span_num .* weights)
end

function number_of_spanwise_panels(wing :: Wing, span_num :: Vector{<: Integer})
    @assert (length ∘ spans)(wing) > 1 "Provide a positive integer of spanwise panels for 1 wing section."
    return span_num
end

# Spacing
function symmetric_spacing(wing :: Wing)
    sins = AbstractSpacing[Sine(1); Sine(0)]
    if length(wing.spans) == 1 
        if wing.symmetry
            return sins
        else 
            return Cosine()
        end
    else
        cosins = fill(Cosine(), (length(spans(wing)) - 1))
        return [ cosins; sins; cosins ]
    end
end

## Panelling
#==========================================================================================#

mesh_chords(wing :: Wing, span_num :: Vector{<: Integer}, chord_num :: Integer; spacings = symmetric_spacing(wing)) = make_panels(chord_coordinates(wing, span_num, chord_num; span_spacing = spacings))

mesh_wing(wing :: Wing, span_num :: Vector{<: Integer}, chord_num :: Integer; spacings = symmetric_spacing(wing)) = make_panels(surface_coordinates(wing, span_num, chord_num, span_spacing = spacings))

mesh_cambers(wing :: Wing, span_num :: Vector{<: Integer}, chord_num :: Integer; spacings = symmetric_spacing(wing)) = make_panels(camber_coordinates(wing, span_num, chord_num; span_spacing = spacings))

function make_panels(wing :: Wing, span_num :: Vector{<: Integer}, chord_num :: Integer; spacings = symmetric_spacing(wing))
    horseshoe_panels = mesh_chords(wing, span_num, chord_num; spacings = spacings)
    camber_panels = mesh_cambers(wing, span_num, chord_num; spacings = spacings)
    horseshoe_panels, normal_vector.(camber_panels)
end

make_panels(wing :: Wing, span_num :: Integer, chord_num :: Integer; spacings = symmetric_spacing(wing)) = make_panels(wing, [span_num], chord_num; spacings = spacings)

panel_wing(comp :: AbstractWing, span_panels :: Union{Integer, Vector{<: Integer}}, chord_panels :: Integer; spacing = symmetric_spacing(comp)) = make_panels(comp, span_panels, chord_panels, spacings = spacing)

## Meshing type for convenience
#==========================================================================================#

"""
    WingMesh(
        wing :: AbstractWing, 
        n_span :: Vector{Integer}, n_chord :: Integer;
        span_spacing :: AbstractSpacing = symmetric_spacing(wing)
    )

Define a container to generate meshes and panels for a given `Wing` with a specified distribution of number of spanwise panels, and a number of chordwise panels.

Optionally a combination of `AbstractSpacing` types (`Sine(), Cosine(), Uniform()`) can be provided to the **named argument** `span_spacing`, either as a singleton or as a vector with length equal to the number of spanwise sections. By default, the combination is `[Sine(), Cosine(), ..., Cosine()]`.

For surface coordinates, the wing mesh will have (n_chord - 1) * 2 chordwise panels from TE-LE-TE and (n_span * 2) spanwise panels.
"""
struct WingMesh{M <: AbstractWing, N <: Integer, P, Q} <: AbstractWing
    surface       :: M
    num_span      :: Vector{N}
    num_chord     :: N
    chord_spacing :: P
    span_spacing  :: Q

    # Main constructor
    WingMesh(surface :: M, n_span :: AbstractVector{N}, n_chord :: N, chord_spacing :: P, span_spacing :: Q) where {M <: AbstractWing, N <: Integer, P <: AbstractSpacing, Q <: Union{AbstractSpacing, Vector{<:AbstractSpacing}}} = new{M,N,P,Q}(surface, n_span, n_chord, chord_spacing, span_spacing)
end

function check_definition(surf :: Wing, n_span) 
    @assert length(n_span) == length(surf.spans) "$(length(n_span)) ≂̸ $(length(surf.spans)), the spanwise number vector's length must be the same as the number of sections of the surface."
end

function WingMesh(surface, n_span, n_chord :: Integer; chord_spacing = Cosine(), span_spacing = symmetric_spacing(surface)) 
    check_definition(surface, n_span)

    if surface.symmetry
        n_span = [ reverse(n_span); n_span ] .÷ 2
    elseif surface.flip
        n_span = reverse(n_span)
    end
    
    WingMesh(surface, n_span, n_chord, chord_spacing, span_spacing)
end

WingMesh(surface, n_span :: Integer, n_chord :: Integer; chord_spacing = Cosine(), span_spacing = symmetric_spacing(surface)) = WingMesh(surface, number_of_spanwise_panels(surface, n_span), n_chord, chord_spacing, span_spacing)

# Forwarding functions for Wing type
MacroTools.@forward WingMesh.surface foils, chords, spans, twists, sweeps, dihedrals, mean_aerodynamic_chord, camber_thickness, leading_edge, trailing_edge, wing_bounds, mean_aerodynamic_center, projected_area, span, position, orientation, affine_transformation, maximum_thickness_to_chord, chop_leading_edge

"""
    chord_coordinates(wing :: WingMesh, n_span = wing.num_span, n_chord = wing.num_chord)

Generate the chord coordinates of a `WingMesh` with default spanwise ``n_s`` and chordwise ``n_c`` panel distributions from the mesh.
"""
chord_coordinates(wing :: WingMesh, n_span = wing.num_span, n_chord = wing.num_chord) = chord_coordinates(wing.surface, n_span, n_chord; span_spacing = wing.span_spacing)

"""
    camber_coordinates(wing :: WingMesh, n_span = wing.num_span, n_chord = wing.num_chord)

Generate the camber coordinates of a `WingMesh` with default spanwise ``n_s`` and chordwise ``n_c`` panel distributions from the mesh.
"""
camber_coordinates(wing :: WingMesh, n_span = wing.num_span, n_chord = wing.num_chord) = camber_coordinates(wing.surface, n_span, n_chord; span_spacing = wing.span_spacing)

"""
    surface_coordinates(wing :: WingMesh, n_span = wing.num_span, n_chord = wing.num_chord)

Generate the surface coordinates of a `WingMesh` with default spanwise ``n_s`` and chordwise ``n_c`` panel distributions from the mesh.
"""
surface_coordinates(wing :: WingMesh, n_span = wing.num_span, n_chord = wing.num_chord) = surface_coordinates(wing.surface, n_span, n_chord, span_spacing = wing.span_spacing)

"""
    surface_panels(
        wing_mesh :: WingMesh, 
        n_s = wing_mesh.num_span, 
        n_c = length(first(foils(wing_mesh.surface))).x
    )

Generate the surface panel distribution from a `WingMesh` with the default spanwise ``n_s`` panel distribution from the mesh and the chordwise panel ``n_c`` distribution from the number of points defining the root airfoil.

In case of strange results, provide a higher number of chordwise panels to represent the airfoils more accurately
""" 
surface_panels(wing :: WingMesh, n_span = wing.num_span, n_chord = length(first(foils(wing.surface)).x))  = (make_panels ∘ surface_coordinates)(wing, n_span, n_chord)

"""
    chord_panels(wing_mesh :: WingMesh)

Generate the chord mesh as a matrix of `Panel3D` from a `WingMesh`.
"""
chord_panels(wing :: WingMesh) = make_panels(chord_coordinates(wing))

"""
    camber_panels(wing_mesh :: WingMesh)

Generate the camber mesh as a matrix of `Panel3D` from a `WingMesh`.
"""
camber_panels(wing :: WingMesh) = make_panels(camber_coordinates(wing))

"""
    wetted_area(
        wing_mesh :: WingMesh, 
        n_s = wing_mesh.num_span, 
        n_c = length(first(foils(wing_mesh.surface))).x
    )

Determine the wetted area ``S_{wet}`` of a `WingMesh` by calculating the total area of the surface panels.
"""
wetted_area(wing :: WingMesh, n_span = wing.num_span, n_chord = length(first(foils(wing.surface)).x)) = wetted_area(surface_panels(wing, n_span, n_chord))

"""
    wetted_area_ratio(
        wing_mesh :: WingMesh, 
        n_s = wing_mesh.num_span, 
        n_c = length(first(foils(wing_mesh.surface))).x
    )

Determine the wetted area ratio ``S_{wet}/S`` of a `WingMesh` by calculating the ratio of the total area of the surface panels to the projected area of the `Wing`.

The wetted area ratio should be slightly above 2 for thin airfoils.
"""
wetted_area_ratio(wing :: WingMesh, n_span = wing.num_span, n_chord = length(first(foils(wing.surface)).x)) = wetted_area(wing, n_span, n_chord) / projected_area(wing.surface)

function Base.show(io :: IO, mesh :: WingMesh)
    n_c, n_s = mesh.num_chord, mesh.num_span # size(mesh.chord_mesh) .- 1
    println(io, "WingMesh —")
    show(io, mesh.surface)
    println(io, "Spanwise panels: ", n_s)
    println(io, "Chordwise panels: ", n_c)
    println(io, "Spanwise spacing: ", mesh.span_spacing)
    println(io, "Chordwise spacing: ", mesh.chord_spacing)

    nothing
end