## Wing definition
#=====================================================================#

# Helper functions
aspect_ratio(span, area) = span^2 / area
mean_geometric_chord(span, area) = area / span
taper_ratio(root_chord, tip_chord) = tip_chord / root_chord
area(span, chord) = span * chord

# Mean aerodynamic chord for a single section
mean_aerodynamic_chord(c_r, λ) = (2/3) * c_r * (1 + λ + λ^2)/(1 + λ)

# Spanwise location of mean aerodynamic chord for a single section
y_mac(y, b, λ) = y + b * (1 + 2λ) / (3(1 + λ))

# Projected area of a single trapezoidal section
section_projected_area(b, c1, c2, t1, t2) = b * (c1 + c2) / 2 * cosd((t1 + t2) / 2)

# Aspect ratio for a single trapezoidal section
aspect_ratio(b, c1, c2) = 2b / (c1 + c2)

# Conversion of sweep angle from leading-edge to any ratio along the chord
sweep_angle(λ, AR, Λ_LE, w) = Λ_LE - atand(2w * (1 - λ), AR * (1 + λ))

# Wing type definition
struct Wing{V,N <: AbstractAffineMap} <: AbstractWing
    foils      :: Vector{<:AbstractFoil}
    chords     :: V
    twists     :: V
    spans      :: V
    dihedrals  :: V
    sweeps     :: V
    affine     :: N
    symmetry   :: Bool
    flip       :: Bool
    # controls   :: Vector{<:AbstractControlSurface}
end

# Default constructor
function Wing(foils, chords, twists, spans, dihedrals, sweeps, affine, symmetry, flip, controls)
    # Error handling
    check_wing(foils, chords, twists, spans, dihedrals, sweeps, controls)

    # TODO: Perform automatic cosine interpolation of foils with minimum number of points for surface construction?
    # foils = cosine_interpolation.(foils, 60)

    T = promote_type(typeof(chords), typeof(twists), typeof(spans), typeof(dihedrals), typeof(sweeps))
    N = typeof(affine)

    # Convert angles to radians, adjust twists to leading edge, and generate Wing
    Wing{T,N}(foils, chords, -twists, spans, dihedrals, sweeps, affine, symmetry, flip)#, controls)
end

# Named arguments version for ease, with default NACA-4 0012 airfoil shape
"""
    Wing(
        chords, 
        foils :: Vector{Foil}, 
        twists, 
        spans, 
        dihedrals, 
        sweeps,
        controls = [],
        symmetry = false,
        flip     = false,
        position = zeros(3),
        angle    = 0.,
        axis     = [0.,1.,0.],
    )

Definition for a `Wing` consisting of ``N+1`` `Foil`s, their associated chord lengths ``c`` and twist angles ``ι``, for ``N`` sections with span lengths ``b``, dihedrals ``δ`` and leading-edge sweep angles ``Λ_{LE}``, with all angles in degrees. Optionally, specify translation and a rotation in angle-axis representation for defining coordinates in a global axis system. Additionally, specify Boolean arguments for symmetry or reflecting in the ``x``-``z`` plane.

# Arguments
- `chords :: Vector{Real}`: Chord lengths (m)
- `foils :: Vector{Foil} = fill(naca4(0,0,1,2), length(chords))`: `Foil` shapes, default is NACA 0012.
- `spans :: Vector{Real} = ones(length(chords) - 1) / (length(chords) - 1)`: Span lengths (m), default yields total span length 1.
- `dihedrals :: Vector{Real} = zero(spans)`: Dihedral angles (deg), default is zero.
- `sweeps :: Vector{Real} = zero(spans)`: Sweep angles (deg), default is zero.
- `chord_ratio :: Real = 0.`: Chord ratio for sweep angle 
                          e.g., 0    = Leading-edge sweep, 
                                1    = Trailing-edge sweep,
                                0.25 = Quarter-chord sweep
- `symmetry :: Bool = false`: Symmetric in the ``x``-``z`` plane
- `flip :: Bool = false`: Flip coordinates in the ``x``-``z`` plane
- `position :: Vector{Real} = zeros(3)`: Position (m)
- `angle :: Real = 0.`: Angle of rotation (degrees)
- `axis :: Vector{Real} = [0.,1.,0.]`: Axis of rotation
- `affine :: AffineMap = AffineMap(AngleAxis(deg2rad(angle), axis...), position)`: Affine mapping for the position and orientation via `CoordinateTransformations.jl` (overrides `angle` and `axis` if specified)
"""
@views function Wing(;
        chords, 
        foils     = fill(naca4(0,0,1,2), length(chords)),
        twists    = zeros(eltype(chords), length(chords)),
        spans     = ones(eltype(chords), length(chords) - 1) / (length(chords) - 1),
        dihedrals = zeros(eltype(chords), length(spans)),
        sweeps    = zeros(eltype(chords), length(spans)),
        # controls  = fill(Flap(0), length(spans)),
        chord_ratio = 0.0,
        position  = zeros(eltype(chords), 3),
        angle     = 0.,
        axis      = eltype(chords)[0., 1., 0.],
        affine    = AffineMap(AngleAxis(deg2rad(angle), axis...), SVector{3}(position)),
        symmetry  = false,
        flip      = false
    )

    # Convert sweep angles to leading-edge
    sweeps = @. sweep_angle(
        chords[2:end] / chords[1:end-1], # Section taper ratios
        2spans / (chords[2:end] + chords[1:end-1]), # Section aspect ratios
        sweeps, # Sweep angles at desired normalized location
        -chord_ratio # Normalized sweep angle location ∈ [0,1]
    )

    Wing(foils, chords, twists, spans, dihedrals, sweeps, affine, symmetry, flip)#, controls)
end

function check_wing(foils, chords, twists, spans, dihedrals, sweeps, controls)
    # Check if number of sections match up with number of edges
    nf, nc, nt = length(foils), length(chords), length(twists)
    nb, nd, ns = length(spans), length(dihedrals), length(spans)
    ncon = length(controls)
    @assert nf == nc == nt "Number of foils, chords and twists specified must be the same!"
    @assert nb == nd == ns "Number of spans, dihedrals, and sweeps specified must be the same!"
    @assert nf == ns + 1 "$(ns+1) foils, chords and twists are required for $ns spanwise section(s)."
    @assert ncon == ns "Number of control surfaces must be the same as number of spanwise section(s)."

    # Check if lengths are positive
    @assert any(x -> x >= zero(eltype(x)), chords) | any(x -> x >= zero(eltype(x)), spans) "Chord and span lengths must be positive."
    # Check if dihedrals and sweeps are within bounds
    @assert all(x -> x > -convert(typeof(x), 90.) && x < convert(typeof(x), 90.), dihedrals) && 
        all(x -> x > -convert(typeof(x), 90.) && x < convert(typeof(x), 90.), sweeps) 
        "Dihedrals and sweep angles must not exceed ±90ᵒ."
end


# Getters
foils(wing :: Wing) = wing.foils
chords(wing :: Wing) = wing.chords
twists(wing :: Wing) = wing.twists
spans(wing :: Wing) = wing.spans
dihedrals(wing :: Wing) = wing.dihedrals

"""
    sweeps(wing :: AbstractWing, w = 0.)

Obtain the sweep angles (in radians) at the corresponding normalized chord length ratio ``w ∈ [0,1]``.
"""
sweeps(wing :: Wing, w = 0.) = @views @. sweep_angle(
    wing.chords[2:end] / wing.chords[1:end-1], # Section tapers
    aspect_ratio(wing.spans, wing.chords[2:end], wing.chords[1:end-1]), # Section aspect ratios
    wing.sweeps, # Section leading-edge sweep angles
    w # Normalized sweep angle location ∈ [0,1]
)

# Affine transformations
Base.position(wing :: Wing) = wing.affine.translation
orientation(wing :: Wing) = wing.affine.linear
affine_transformation(wing :: Wing) = wing.affine

(f :: AffineMap)(wing :: Wing) = @set wing.affine = f ∘ wing.affine

"""
    span(wing :: Wing)

Compute the planform span of a `Wing`.
"""
function span(wing :: Wing)
    b = sum(wing.spans)
    ifelse(wing.symmetry, 2b, b)
end

"""
    taper_ratio(wing :: Wing)

Compute the taper ratio of a `Wing`, defined as the tip chord length divided by the root chord length, independent of the number of sections.
"""
taper_ratio(wing :: Wing) = last(wing.chords) / first(wing.chords)

"""
    aspect_ratio(wing :: AbstractWing)

Compute the aspect ratio of an `AbstractWing`.
"""
aspect_ratio(wing) = aspect_ratio(span(wing), projected_area(wing))

"""
    properties(wing :: AbstractWing)

Compute the generic properties of interest (span, area, etc.) of an `AbstractWing`.
"""
properties(wing :: AbstractWing) = [ aspect_ratio(wing), span(wing), projected_area(wing), mean_aerodynamic_chord(wing), mean_aerodynamic_center(wing) ]

function section_projected_areas(wing :: Wing)
    spans = ifelse(wing.symmetry, 2 * wing.spans, wing.spans)
    @views section_projected_area.(spans, wing.chords[1:end-1], wing.chords[2:end], wing.twists[1:end-1], wing.twists[2:end])
end

"""
    projected_area(wing :: Wing)

Compute the projected area (onto the spanwise plane) of a `Wing`.
"""
projected_area(wing :: Wing) = sum(section_projected_areas(wing))

section_macs(wing :: Wing) = @views @. mean_aerodynamic_chord(wing.chords[1:end-1], wing.chords[2:end] / wing.chords[1:end-1])

"""
    mean_aerodynamic_chord(wing :: Wing)

Compute the mean aerodynamic chord of a `Wing`.
"""
function mean_aerodynamic_chord(wing :: Wing)
    areas = section_projected_areas(wing)
    macs = section_macs(wing)
    sum(macs .* areas) / sum(areas)
end

"""
    mean_aerodynamic_center(wing :: Wing, 
        factor = 0.25; 
        symmetry = wing.symmetry, 
        flip = wing.flip
    )

Compute the mean aerodynamic center of a `Wing`. By default, the factor is assumed to be at 25% from the leading edge, which can be adjusted. Similarly, options are provided to account for symmetry or to flip the location in the ``x``-``z`` plane.
"""
@views function mean_aerodynamic_center(wing :: Wing, factor = 0.25; symmetry = wing.symmetry, flip = wing.flip)
    # Compute mean aerodynamic chords and projected areas
    macs = section_macs(wing)
    areas = section_projected_areas(wing)

    # Get leading edge coordinates
    wing_LE = leading_edge(wing)
    x_LEs = wing_LE[:,1]
    y_LEs = wing_LE[:,2]

    # Compute x-y locations of section MACs
    x_mac_LEs = @views @. y_mac(x_LEs[1:end-1], x_LEs[2:end], wing.chords[2:end] / wing.chords[1:end-1])
    y_macs = @views @. y_mac(y_LEs[1:end-1], wing.spans, wing.chords[2:end] / wing.chords[1:end-1])

    # Calculate section MAC coords
    mac_coords = @. SVector(x_mac_LEs + factor * macs, y_macs, 0.)

    # Calculate weighted MAC by section areas
    mac = sum(mac_coords .* areas) / sum(areas)

    # Symmetry/flipping adjustments
    if symmetry
        mac = SVector(mac[1], 0., mac[3])
    elseif flip
        mac = SVector(mac[1], -mac[2], mac[3])
    end

    # Do affine transformation
    return affine_transformation(wing)(mac)
end

"""
    camber_thickness(wing :: Wing, num :: Integer)

Compute the camber-thickness distribution at each spanwise intersection of a `Wing`. A `num` must be specified to interpolate the internal `Foil` coordinates, which affects the accuracy of ``(t/c)ₘₐₓ`` accordingly.
"""
camber_thickness(wing :: Wing, num :: Integer) = camber_thickness.(wing.foils, num)

"""
    maximum_thickness_to_chord(wing :: Wing, num :: Integer)

Compute the maximum thickness-to-chord ratios ``(t/c)ₘₐₓ`` and their locations ``(x/c)`` at each spanwise intersection of a `Wing`. 

Returns an array of pairs ``[(x/c),(t/c)ₘₐₓ]``, in which the first entry of each pair is the location (normalized to the local chord length at the spanwise intersection) and the corresponding maximum thickness-to-chord ratio at the intersection.

A `num` must be specified to interpolate the internal `Foil` coordinates, which affects the accuracy of ``(t/c)ₘₐₓ`` accordingly.
"""
maximum_thickness_to_chord(wing :: Wing, num :: Integer) = map(maximum_thickness_to_chord, camber_thickness(wing, num))

function wing_bounds(wing :: Wing)
    # Compute y-z points
    spans, dihedrals, sweeps = wing.spans, wing.dihedrals, wing.sweeps
    cum_spans = [ 0; cumsum(spans) ]
    sweeped_spans = [ 0; cumsum(@. spans * tand(sweeps)) ]
    dihedraled_spans = [ 0; cumsum(@. spans * tand(dihedrals)) ]

    # Compute x points
    chords = wing.chords
    twisted_chords  = @. chords * sind(wing.twists)

    # Leading edge
    le = [ sweeped_spans cum_spans dihedraled_spans ]

    # Trailing edge perturbation
    te = [ chords (zeros ∘ length)(chords) twisted_chords ]

    # Get actual bounds by sum-cumming, hehe I wanna die
    # Indexing: [section,xyz,le/te]
    bounds = cumsum(cat(le, te, dims = Val(3)), dims = 3)

    return bounds
end

"""
    leading_edge(wing :: Wing)

Compute the leading edge coordinates of a `Wing`.
"""
leading_edge(wing :: Wing) = @views wing_bounds(wing)[:,:,1]

"""
    trailing_edge(wing :: Wing)

Compute the trailing edge coordinates of a `Wing`.
"""
trailing_edge(wing :: Wing) = @views wing_bounds(wing)[:,:,2]


function Base.show(io :: IO, wing :: AbstractWing)
    sym = ifelse(wing.symmetry, "Symmetric ", "")
    println(io, sym, supertype(typeof(wing)),  " with ", length(spans(wing)), " spanwise section(s).")
    println(io, "Aspect Ratio: ", aspect_ratio(wing))
    println(io, "Span (m): ", span(wing))
    println(io, "Projected Area (m²): ", projected_area(wing))
    println(io, "Mean Aerodynamic Chord (m): ", mean_aerodynamic_chord(wing))
    println(io, "Mean Aerodynamic Center (m): ", mean_aerodynamic_center(wing))
    println(io, "Position (m): ", wing.affine.translation)

    nothing
end