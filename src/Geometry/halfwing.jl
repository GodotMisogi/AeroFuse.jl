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
    foils
    chords   
    twists   
    spans    
    dihedrals
    sweeps   
    function HalfWing(foils :: Vector{Foil{T}}, chords, twists, spans, dihedrals, sweeps) where T <: Real
        # Checking if sections match up with edges
        num_sections = maximum(length.([spans, dihedrals, sweeps]))
        @assert all(length.([foils, chords, twists]) .== num_sections + 1) "$(num_sections + 1) foils, chords and twists are required for $(num_sections) section(s)."
        new{T}(foils, chords, -deg2rad.(twists), spans, deg2rad.(dihedrals), deg2rad.(sweeps))
    end
end

# Default NACA0012 foil shape, compatible with ReverseDiff
HalfWing(chords :: AbstractVector{<: Real}, twists :: AbstractVector{<: Real}, spans :: AbstractVector{<: Real}, dihedrals :: AbstractVector{<: Real}, sweeps :: AbstractVector{<: Real}) = HalfWing(fill(Foil(naca4((0,0,1,2))), length(chords)), chords, twists, spans, dihedrals, sweeps)

"""
    span(half_wing :: HalfWing)

Compute the planform span of a `HalfWing`.
"""
span(wing :: HalfWing) = sum(wing.spans)

"""
    projected_area(half_wing :: HalfWing)
    
Compute the projected area of a `HalfWing` onto the ``x``-``y`` plane.
"""
function projected_area(wing :: HalfWing)
    mean_chords = fwdsum(wing.chords) / 2 # Mean chord lengths of sections.
    mean_twists = fwdsum(wing.twists) / 2 # Mean twist angles of sections
    sum(@. wing.spans * mean_chords * cos(mean_twists))
end

section_macs(wing :: HalfWing) = mean_aerodynamic_chord.(wing.chords[1:end-1], fwddiv(wing.chords))

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
    x_mac_LEs	= y_mac.(x_LEs[1:end-1], 2 .* fwddiff(x_LEs), fwddiv(wing.chords))
    y_macs		= y_mac.(y_LEs[1:end-1], wing.spans, fwddiv(wing.chords))

    mac_coords 	= SVector.(x_mac_LEs .+ factor .* macs, y_macs, (zeros ∘ length)(y_macs))

    sum(mac_coords .* areas) / sum(areas)
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

