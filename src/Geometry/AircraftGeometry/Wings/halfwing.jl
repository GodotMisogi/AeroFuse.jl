"""
    Wing(foils :: Vector{Foil}, 
             chords, 
             twists, 
             spans, 
             dihedrals, 
             sweeps,
             position = zeros(3),
             angle    = 0.
             axis     = [0.,1.,0.])

Definition for a `Wing` consisting of ``N+1`` `Foil`s, their associated chord lengths ``c`` and twist angles ``ι``, for ``N`` sections with span lengths ``b``, dihedrals ``δ`` and leading-edge sweep angles ``Λ_{LE}``, with all angles in degrees.
"""
struct Wing{T <: Number, V <: AbstractVector{T}, N <: AbstractAffineMap} <: AbstractWing
    foils      :: Vector{<: AbstractFoil}
    chords     :: V
    twists     :: V
    spans      :: V
    dihedrals  :: V
    sweeps     :: V
    affine     :: N
    symmetry   :: Bool
    flip       :: Bool

    # Default constructor
    function Wing(foils, chords, twists, spans, dihedrals, sweeps, affine, w_sweep = 0., symmetry = false, flip = false)
        # Error handling
        # check_wing(foils, chords, twists, spans, dihedrals, sweeps)

        # Convert sweep angles to leading-edge
        sweeps = @. sweep_angle(
            chords[2:end] / chords[1:end-1], # Section tapers
            aspect_ratio(spans, chords[2:end], chords[1:end-1]), # Section aspect ratios
            sweeps, # Sweep angles at desired normalized location
            -w_sweep # Normalized sweep angle location ∈ [0,1]
        )

        # TODO: Perform automatic cosine interpolation of foils with minimum number of points for surface construction?
        # foils = cosine_interpolation.(foils, 60)

        T = promote_type(eltype(chords), eltype(twists), eltype(spans), eltype(dihedrals), eltype(sweeps), typeof(w_sweep))
        N = typeof(affine)
        V = promote_type(typeof(chords),  typeof(twists), typeof(spans), typeof(dihedrals), typeof(sweeps))
        # S = promote_type(typeof(symmetry), typeof(flip))

        # Convert angles to radians, adjust twists to leading edge, and generate Wing
        new{T,V,N}(foils, chords, -twists, spans, dihedrals, sweeps, affine, symmetry, flip)
    end
end

# Named arguments version for ease, with default NACA-4 0012 airfoil shape
function Wing(;
    chords, 
    twists    = zero(chords),
    spans     = ones(length(chords) - 1) / (length(chords) - 1),
    dihedrals = zero(spans),
    sweeps    = zero(spans),
    foils     = fill(naca4(0,0,1,2), length(chords)),
    w_sweep   = 0.0,
    position  = zeros(3),
    angle     = 0.,
    axis      = [0., 1., 0.],
    affine    = AffineMap(QuatRotation(AngleAxis(deg2rad(angle), axis...)), position),
    symmetry  = false,
    flip      = false
)

    Wing(foils, chords, twists, spans, dihedrals, sweeps, affine, w_sweep, symmetry, flip)
end

function check_wing(foils, chords, twists, spans, dihedrals, sweeps)
    # Check if number of sections match up with number of edges (NEEDS WORK)
    @assert (length ∘ zip)(foils, chords, twists) == (length ∘ zip)(spans, dihedrals, sweeps) + 1 "N+1 foils, chords and twists are required for N section(s)."
    # Check if lengths are positive
    @assert any(x -> x >= 0., chords) || any(x -> x >= 0., spans) "Chord and span lengths must be positive."
    # Check if dihedrals and sweeps are within bounds
    @assert any(x -> x >= -90. || x <= 90., dihedrals) || any(x -> x >= -90. || x <= 90., sweeps) "Dihedrals and sweep angles must not exceed ±90ᵒ."
end


# Getters
foils(wing     :: Wing) = wing.foils
chords(wing    :: Wing) = wing.chords
twists(wing    :: Wing) = wing.twists
spans(wing     :: Wing) = wing.spans
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

aspect_ratio(b, c1, c2) = 2b / (c1 + c2)
sweep_angle(λ, AR, Λ_LE, w) = Λ_LE - atand(2w * (1 - λ), AR * (1 + λ))

# Affine transformations
Base.position(wing :: Wing) = wing.affine.translation
orientation(wing :: Wing) = wing.affine.linear
affine_transformation(wing :: Wing) = wing.affine

(f :: AffineMap)(wing :: Wing) = @set wing.affine = f ∘ wing.affine

function span(wing :: Wing)
    b = sum(wing.spans)
    ifelse(wing.symmetry, 2b, b)
end

taper_ratio(wing :: Wing) = last(wing.chords) / first(wing.chords)

section_projected_area(b, c1, c2, t1, t2) = b * (c1 + c2) / 2 * cosd((t1 + t2) / 2)

function section_projected_areas(wing :: Wing)
    spans = ifelse(wing.symmetry, 2 * wing.spans, wing.spans)
    @views section_projected_area.(spans, wing.chords[1:end-1], wing.chords[2:end], wing.twists[1:end-1], wing.twists[2:end])
end

projected_area(wing :: Wing) = sum(section_projected_areas(wing))

section_macs(wing :: Wing) = @views @. mean_aerodynamic_chord(wing.chords[1:end-1], wing.chords[2:end] / wing.chords[1:end-1])

function mean_aerodynamic_chord(wing :: Wing)
    areas = section_projected_areas(wing)
    macs = section_macs(wing)
    sum(macs .* areas) / sum(areas)
end

function mean_aerodynamic_center(wing :: Wing, factor = 0.25; symmetry = wing.symmetry, flip = wing.flip)
    # Compute mean aerodynamic chords and projected areas
    macs = section_macs(wing)
    areas = section_projected_areas(wing)

    # Get leading edge coordinates
    wing_LE = combinedimsview(leading_edge(wing), (1))
    x_LEs   = @views wing_LE[:,1]
    y_LEs   = @views wing_LE[:,2]

    # Compute x-y locations of section MACs
    x_mac_LEs = @views @. y_mac(x_LEs[1:end-1], 2 * x_LEs[2:end], wing.chords[2:end] / wing.chords[1:end-1])
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

camber_thickness(wing :: Wing, num :: Integer) = camber_thickness.(wing.foils, num)

maximum_thickness_to_chord(wing :: Wing, num :: Integer) = maximum_thickness_to_chord.(camber_thickness(wing, num))

function max_thickness_to_chord_ratio_sweeps(wing :: Wing, num)
    # Compute (t/c)_max locations and values
    xs_max_tbyc = maximum_thickness_to_chord(wing, num)
    max_tbyc = getindex.(xs_max_tbyc, 2)
    xs_temp = getindex.(xs_max_tbyc, 1)
    
    # Get leading edge
    le = leading_edge(wing)

    # Determine x-coordinates in geometry frame
    xs = @. getindex(le, 1) + wing.chords * getindex(xs_temp, 1)

    # Compute sectional x-coordinates
    ds = xs[2:end] - xs[1:end-1]
    
    # Compute leading-edge sweep angles accounting for dihedral
    widths = @. wing.spans / cosd(wing.dihedrals)
    sweeps = @. atand(ds, widths)

    # Averaging for sections
    xs = (xs[1:end-1] + xs[2:end]) / 2
    tbycs = (max_tbyc[1:end-1] + max_tbyc[2:end]) / 2

    xs, tbycs, sweeps
end

"""
    wing_bounds(wing :: Wing, flip = false)

Compute the leading and trailing edge coordinates of a `Wing`, with an option to flip the signs of the ``y``-coordinates.
"""
function wing_bounds(wing :: Wing)
    # Compute y-z points
    spans, dihedrals, sweeps = wing.spans, wing.dihedrals, wing.sweeps
    cum_spans = [ 0; cumsum(spans) ]
    sweeped_spans = [ 0; cumsum(@. spans * tand(sweeps)) ]
    dihedraled_spans = [ 0; cumsum(@. spans * tand(dihedrals)) ]

    # Compute x points
    chords          = wing.chords
    twisted_chords  = @. chords * sind(wing.twists)

    # Leading edge
    le = [ sweeped_spans cum_spans dihedraled_spans ]

    # Trailing edge perturbation
    te = [ chords (zeros ∘ length)(chords) twisted_chords ]

    # Get actual bounds by sum-cumming, hehe I wanna die
    # Indexing: [section,xyz,le/te]
    bounds = cumsum([ le ;;; te ], dims = 3)

    return splitdimsview(bounds, (3,1))
end

"""
    trailing_edge(wing :: Wing, flip = false)

Compute the trailing edge coordinates of a `Wing`, with an option to flip the signs of the ``y``-coordinates.
"""
leading_edge(wing :: Wing) = @views wing_bounds(wing)[1,:]
trailing_edge(wing :: Wing) = @views wing_bounds(wing)[2,:]


function Base.show(io :: IO, wing :: AbstractWing)
    sym = ifelse(wing.symmetry, "Symmetric ", "")
    println(io, sym, supertype(typeof(wing)),  " with ", length(spans(wing)), " spanwise section(s).")
    println(io, "Aspect Ratio: ", aspect_ratio(wing))
    println(io, "Span (m): ", span(wing))
    println(io, "Projected Area (m): ", projected_area(wing))
    println(io, "Mean Aerodynamic Chord (m): ", mean_aerodynamic_chord(wing))
    println(io, "Mean Aerodynamic Center (m): ", mean_aerodynamic_center(wing))

    nothing
end