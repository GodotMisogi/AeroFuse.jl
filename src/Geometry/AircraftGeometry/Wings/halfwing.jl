"""
    HalfWing(foils :: Vector{Foil}, 
             chords, 
             twists, 
             spans, 
             dihedrals, 
             sweeps,
             position = zeros(3),
             angle    = 0.
             axis     = [0.,1.,0.])

Definition for a `HalfWing` consisting of ``N+1`` `Foil`s, their associated chord lengths ``c`` and twist angles ``ι``, for ``N`` sections with span lengths ``b``, dihedrals ``δ`` and leading-edge sweep angles ``Λ_{LE}``, with all angles in degrees.
"""
struct HalfWing{T <: Real, N <: AbstractAffineMap} <: AbstractWing
    foils      :: Vector{<: AbstractFoil}
    chords     :: Vector{T}
    twists     :: Vector{T}
    spans      :: Vector{T}
    dihedrals  :: Vector{T}
    sweeps     :: Vector{T}
    affine     :: N
end

function HalfWing(foils, chords, twists, spans, dihedrals, sweeps, w_sweep = 0., position = zeros(3), angle = 0., axis = [0.,1.,0.], affine = AffineMap(AngleAxis(deg2rad(angle), axis...), position))
    # Error handling
    check_wing(foils, chords, twists, spans, dihedrals, sweeps)

    # Convert sweep angles to leading-edge
    sweeps = @. sweep_angle(
                    chords[2:end] / chords[1:end-1], # Section tapers
                    aspect_ratio(spans, chords[2:end], chords[1:end-1]), # Section aspect ratios
                    deg2rad(sweeps), # Sweep angles at desired normalized location
                    -w_sweep # Normalized sweep angle location ∈ [0,1]
                )

    # TODO: Perform automatic cosine interpolation of foils with minimum number of points for surface construction?
    # foils = cosine_interpolation.(foils, 60)

    T = promote_type(eltype(chords), eltype(twists), eltype(spans), eltype(dihedrals), eltype(sweeps), typeof(w_sweep))
    N = typeof(affine)

    # Convert angles to radians, adjust twists to leading edge, and generate HalfWing
    HalfWing{T,N}(foils, chords, -deg2rad.(twists), spans, deg2rad.(dihedrals), sweeps, affine)
end

function check_wing(foils, chords, twists, spans, dihedrals, sweeps)
    # Check if number of sections match up with number of edges (NEEDS WORK)
    @assert (length ∘ zip)(foils, chords, twists) == (length ∘ zip)(spans, dihedrals, sweeps) + 1 "N+1 foils, chords and twists are required for N section(s)."
    # Check if lengths are positive
    @assert any(x -> x >= 0., chords) || any(x -> x >= 0., spans) "Chord and span lengths must be positive."
    # Check if dihedrals and sweeps are within bounds
    @assert any(x -> x >= -90. || x <= 90., dihedrals) || any(x -> x >= -90. || x <= 90., sweeps) "Dihedrals and sweep angles must not exceed ±90ᵒ."
end

# Named arguments version for ease, with default NACA-4 0012 airfoil shape
function HalfWing(;
        chords, 
        twists    = zero(chords),
        spans     = ones(length(chords) - 1) / (length(chords) - 1),
        dihedrals = zero(spans),
        sweeps    = zero(spans),
        foils     = fill(naca4(0,0,1,2), length(chords)),
        w_sweep   = 0.0,
        position  = zeros(3),
        angle     = 0.,
        axis      = [0., 1., 0.]
    )

    HalfWing(foils, chords, twists, spans, dihedrals, sweeps, w_sweep, position, angle, axis)
end

# Getters
foils(wing     :: HalfWing) = wing.foils
chords(wing    :: HalfWing) = wing.chords
twists(wing    :: HalfWing) = wing.twists
spans(wing     :: HalfWing) = wing.spans
dihedrals(wing :: HalfWing) = wing.dihedrals

"""
    sweeps(wing :: AbstractWing, w = 0.)

Obtain the sweep angles (in radians) at the corresponding normalized chord length ratio ``w ∈ [0,1]``.
"""
sweeps(wing :: HalfWing, w = 0.) =
    @views @. sweep_angle(
                    wing.chords[2:end] / wing.chords[1:end-1], # Section tapers
                    aspect_ratio(wing.spans, wing.chords[2:end], wing.chords[1:end-1]), # Section aspect ratios
                    wing.sweeps, # Section leading-edge sweep angles
                    w # Normalized sweep angle location ∈ [0,1]
                )

aspect_ratio(b, c1, c2) = 2b / (c1 + c2)
sweep_angle(λ, AR, Λ_LE, w) = Λ_LE - atan(2w * (1 - λ), AR * (1 + λ))

# Affine transformation
Base.position(wing :: HalfWing) = wing.affine.transformation
orientation(wing :: HalfWing) = wing.affine.linear
affine_transformation(wing :: HalfWing) = wing.affine


span(wing :: HalfWing) = sum(wing.spans)

taper_ratio(wing :: HalfWing) = last(wing.chords) / first(wing.chords)

function projected_area(wing :: HalfWing)
    mean_chords = forward_sum(wing.chords) / 2 # Mean chord lengths of sections.
    mean_twists = forward_sum(wing.twists) / 2 # Mean twist angles of sections
    sum(@. wing.spans * mean_chords * cos(mean_twists))
end

section_macs(wing :: HalfWing) = @views mean_aerodynamic_chord.(wing.chords[1:end-1], forward_division(wing.chords))

section_projected_areas(wing :: HalfWing) = wing.spans .* forward_sum(wing.chords) / 2

function mean_aerodynamic_chord(wing :: HalfWing)
    areas = section_projected_areas(wing)
    macs  = section_macs(wing)
    sum(macs .* areas) / sum(areas)
end

function mean_aerodynamic_center(wing :: HalfWing, factor = 0.25)
    # Compute mean aerodynamic chords and projected areas
    macs        = section_macs(wing)
    areas       = section_projected_areas(wing)

    # Computing leading edge coordinates
    wing_LE     = leading_edge(wing)
    x_LEs       = getindex.(wing_LE, 1)
    y_LEs       = getindex.(wing_LE, 2)

    # Compute x-y locations of MACs
    x_mac_LEs   = @views @. y_mac(x_LEs[1:end-1], 2 * x_LEs[2:end], wing.chords[2:end] / wing.chords[1:end-1])
    y_macs      = @views @. y_mac(y_LEs[1:end-1], wing.spans, wing.chords[2:end] / wing.chords[1:end-1])

    mac_coords  = @. SVector(x_mac_LEs + factor * macs, y_macs, 0.)

    affine_transformation(wing)(sum(mac_coords .* areas) / sum(areas))
end

camber_thickness(wing :: HalfWing, num :: Integer) = camber_thickness.(wing.foils, num)

maximum_thickness_to_chord(wing :: HalfWing, num :: Integer) = maximum_thickness_to_chord.(camber_thickness(wing, num))

function max_thickness_to_chord_ratio_sweeps(wing :: HalfWing, num)
    # Compute (t/c)_max locations and values
    xs_max_tbyc = maximum_thickness_to_chord(wing, num)
    max_tbyc    = getindex.(xs_max_tbyc, 2)
    xs_temp     = getindex.(xs_max_tbyc, 1)
    
    # Get leading edge
    le          = leading_edge(wing)

    # Determine x-coordinates in Cartesian frame
    xs          = @. getindex(le, 1) + wing.chords * getindex(xs_temp, 1)

    # Compute sectional x-coordinates
    ds          = forward_difference(xs)
    
    # Compute leading-edge sweep angles accounting for dihedral
    widths      = @. wing.spans / cos(wing.dihedrals)
    sweeps      = @. atan(ds, widths)

    # Averaging for sections
    xs          = forward_sum(xs_temp) / 2
    tbycs       = forward_sum(max_tbyc) / 2

    xs, tbycs, sweeps
end

"""
    leading_edge(wing :: HalfWing, flip = false)

Compute the leading edge coordinates of a `HalfWing`, with an option to flip the signs of the ``y``-coordinates.
"""
function leading_edge(wing :: HalfWing, flip = false)
    spans, dihedrals, sweeps = wing.spans, wing.dihedrals, wing.sweeps

    sweeped_spans    = [ 0; cumsum(@. spans * tan(sweeps)) ]
    dihedraled_spans = [ 0; cumsum(@. spans * tan(dihedrals)) ]
    cum_spans        = [ 0; cumsum(spans) ]

    leading          = SVector.(sweeped_spans, ifelse(flip, -cum_spans, cum_spans), dihedraled_spans)

    ifelse(flip, leading[end:-1:1], leading)
end

"""
    wing_bounds(wing :: HalfWing, flip = false)

Compute the leading and trailing edge coordinates of a `HalfWing`, with an option to flip the signs of the ``y``-coordinates.
"""
function wing_bounds(wing :: HalfWing, flip = false)
    chords          = wing.chords
    twisted_chords  = @. chords * sin(wing.twists)

    leading         = leading_edge(wing, flip)
    trailing        = SVector.(chords, (zeros ∘ length)(chords), twisted_chords)

    shifted_trailing = ifelse(flip, trailing[end:-1:1], trailing) .+ leading

    leading, shifted_trailing
end

"""
    trailing_edge(wing :: HalfWing, flip = false)

Compute the trailing edge coordinates of a `HalfWing`, with an option to flip the signs of the ``y``-coordinates.
"""
trailing_edge(wing :: HalfWing, flip = false) = wing_bounds(wing, flip)[2]