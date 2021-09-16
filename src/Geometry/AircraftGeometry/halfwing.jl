"""
    HalfWing(foils      :: Vector{Foil},
             chords     :: Vector{Real},
             twists     :: Vector{Real},
             spans      :: Vector{Real},
             dihedrals  :: Vector{Real},
             sweeps     :: Vector{Real})

Definition for a `HalfWing` consisting of ``N+1`` airfoils and their associated chord lengths ``c``, twist angles ``\\iota``, for ``N`` sections with span lengths ``b``, dihedrals ``\\delta`` and sweep angles ``\\Lambda``, with all angles in degrees.
"""
struct HalfWing{T <: Real} <: Aircraft
    foils      :: Vector{Foil{T}}
    chords     :: Vector{T}
    twists     :: Vector{T}
    spans      :: Vector{T}
    dihedrals  :: Vector{T}
    sweeps     :: Vector{T}
    position   :: SVector{3,T}
    orientation :: SMatrix{3,3,T}
    function HalfWing(foils :: AbstractVector{Foil{T}}, chords :: AbstractVector{T}, twists :: AbstractVector{T}, spans :: AbstractVector{T}, dihedrals :: AbstractVector{T}, sweeps :: AbstractVector{T}, position = zeros(3), angle = 0., axis = [0.,1.,0.]) where T <: Real
        # Error handling
        check_wing(foils, chords, twists, spans, dihedrals, sweeps)
        # Convert angles to radians, with adjusting twists to leading edge, and generate HalfWing
        new{T}(foils, chords, -deg2rad.(twists), spans, deg2rad.(dihedrals), deg2rad.(sweeps), position, AngleAxis{T}(deg2rad(angle), axis...))
    end
end

# HalfWing(foils :: AbstractVector{<: Foil}, chords :: AbstractVector{<: Real}, twists :: AbstractVector{<: Real}, spans :: AbstractVector{<: Real}, dihedrals :: AbstractVector{<: Real}, sweeps :: AbstractVector{<: Real}) = HalfWing(foils, chords, twists, spans, dihedrals, sweeps)

function check_wing(foils, chords, twists, spans, dihedrals, sweeps)
    # Check if number of sections match up with number of edges
    @assert (length ∘ zip)(foils, chords, twists) == (length ∘ zip)(spans, dihedrals, sweeps) + 1 "N+1 foil sections, chords and twists are required for N section(s)."
    # Check if lengths are positive
    @assert any(x -> x >= 0., chords) || any(x -> x >= 0., spans) "Chord and span lengths must be positive."
    # Check if dihedrals and sweeps are within bounds
    @assert any(x -> x >= -90. || x <= 90., dihedrals) || any(x -> x >= -90. || x <= 90., sweeps) "Dihedrals and sweep angles must not exceed ±90ᵒ."
end

# Named arguments version for ease, with default NACA-4 0012 airfoil shape
HalfWing(; chords, twists, spans, dihedrals, sweep_LEs, foils = fill(Foil(naca4(0,0,1,2), "NACA 0012"), length(chords)), position = zeros(3), angle = 0., axis = [1.,0.,0.]) = HalfWing(foils, chords, twists, spans, dihedrals, sweep_LEs, position, angle, axis)

HalfWingSection(; span = 1., dihedral = 0., sweep_LE = 0., taper = 1., root_chord = 1., root_twist = 0., tip_twist = 0., root_foil = naca4((0,0,1,2)), tip_foil = naca4((0,0,1,2)), position = zeros(3), angle = 0., axis = [1.,0.,0.]) = HalfWing([ Foil(root_foil, "Root"), Foil(tip_foil, "Tip") ], [root_chord, taper * root_chord], [root_twist, tip_twist], [span], [dihedral], [sweep_LE], position, angle, axis)

# Getters
foils(wing     :: HalfWing) = wing.foils
chords(wing    :: HalfWing) = wing.chords
twists(wing    :: HalfWing) = wing.twists
spans(wing     :: HalfWing) = wing.spans
dihedrals(wing :: HalfWing) = wing.dihedrals
sweeps(wing    :: HalfWing) = wing.sweeps

# Affine transformation
Base.position(wing :: HalfWing) = wing.position
orientation(wing :: HalfWing) = wing.orientation
affine_transformation(wing :: HalfWing{T}) where T <: Real = Translation(position(wing)) ∘ LinearMap(orientation(wing))

"""
    span(half_wing :: HalfWing)

Compute the planform span of a `HalfWing`.
"""
span(wing :: HalfWing) = sum(wing.spans)

taper_ratio(wing :: HalfWing) = last(wing.chords) / first(wing.chords)

"""
    projected_area(half_wing :: HalfWing)

Compute the projected area of a `HalfWing` onto the ``x``-``y`` plane.
"""
function projected_area(wing :: HalfWing)
    mean_chords = fwdsum(wing.chords) / 2 # Mean chord lengths of sections.
    mean_twists = fwdsum(wing.twists) / 2 # Mean twist angles of sections
    sum(@. wing.spans * mean_chords * cos(mean_twists))
end

section_macs(wing :: HalfWing) = @views mean_aerodynamic_chord.(wing.chords[1:end-1], fwddiv(wing.chords))

section_projected_areas(wing :: HalfWing) = wing.spans .* fwdsum(wing.chords) / 2

"""
    mean_aerodynamic_chord(half_wing :: HalfWing)

Compute the mean aerodynamic chord of a `HalfWing`.
"""
function mean_aerodynamic_chord(wing :: HalfWing)
    areas = section_projected_areas(wing)
    macs  = section_macs(wing)
    sum(macs .* areas) / sum(areas)
end

function mean_aerodynamic_center(wing :: HalfWing, factor = 0.25)
    # Compute mean aerodynamic chords and projected areas
    macs 		= section_macs(wing)
    areas 		= section_projected_areas(wing)

    # Computing leading edge coordinates
    wing_LE 	= leading_edge(wing)
    x_LEs 		= getindex.(wing_LE, 1)
    y_LEs		= getindex.(wing_LE, 2)

    # Compute x-y locations of MACs
    x_mac_LEs	= @views @. y_mac.(x_LEs[1:end-1], 2 * x_LEs[2:end], wing.chords[2:end] / wing.chords[1:end-1])
    y_macs		= @views @. y_mac.(y_LEs[1:end-1], wing.spans, wing.chords[2:end] / wing.chords[1:end-1])

    mac_coords 	= @. SVector(x_mac_LEs + factor * macs, y_macs, 0.)

    affine_transformation(wing)(sum(mac_coords .* areas) / sum(areas))
end

camber_thickness(wing :: HalfWing, num) = camber_thickness.(wing.foils, num)

max_thickness_to_chord_ratio_location(wing :: HalfWing, num) = max_thickness_to_chord_ratio_location.(camber_thickness(wing, num))

function max_tbyc_sweeps(wing :: HalfWing, num)
    xs_max_tbyc = max_thickness_to_chord_ratio_location(wing, num)
    max_tbyc 	= last.(xs_max_tbyc)
    xs_temp 	= first.(xs_max_tbyc)
    xs 			= first.(leading_edge(wing)) .+ wing.chords .* first.(xs_temp)
    ds 			= fwddiff(xs)
    widths 		= @. wing.spans / cos(wing.dihedrals)

    sweeps 		= @. atan(ds, widths)
    xs 	  		= fwdsum(xs_temp) / 2
    tbycs 		= fwdsum(max_tbyc) / 2

    xs, tbycs, sweeps
end

"""
    leading_edge(wing :: HalfWing, flip = false)

Compute the leading edge coordinates of a `HalfWing`, with an option to flip the signs of the ``y``-coordinates.
"""
function leading_edge(wing :: HalfWing, flip = false)
    spans, dihedrals, sweeps = wing.spans, wing.dihedrals, wing.sweeps

    sweeped_spans 	 = [ 0; cumsum(@. spans * tan(sweeps)) ]
    dihedraled_spans = [ 0; cumsum(@. spans * tan(dihedrals)) ]
    cum_spans 		 = [ 0; cumsum(spans) ]

    leading 		 = SVector.(sweeped_spans, ifelse(flip, -cum_spans, cum_spans), dihedraled_spans)

    ifelse(flip, leading[end:-1:1], leading)
end

"""
    wing_bounds(wing :: HalfWing, flip = false)

Compute the leading and trailing edge coordinates of a `HalfWing`, with an option to flip the signs of the ``y``-coordinates.
"""
function wing_bounds(wing :: HalfWing, flip = false)
    chords 	 		= wing.chords
    twisted_chords 	= @. chords * sin(wing.twists)

    leading  		= leading_edge(wing, flip)
    trailing 		= SVector.(chords, (zeros ∘ length)(chords), twisted_chords)

    shifted_trailing = ifelse(flip, trailing[end:-1:1], trailing) .+ leading

    leading, shifted_trailing
end

trailing_edge(wing :: HalfWing, flip = false) = wing_bounds(wing, flip)[2]

