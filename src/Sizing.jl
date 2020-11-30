module Sizing

using StaticArrays

export dynamic_pressure, span_efficiency_factor, wing_loading_stall_speed,
Wing, HalfWing, aspect_ratio, mean_geometric_chord, mean_aerodynamic_chord, taper_ratio, area, span, quarter_chord, projected_area, wing_bounds

# Generic
dynamic_pressure(ρ, V) = 0.5 * ρ * V^2
stall_speed(wing_loading, CL_max, ρ = 1.225) = 2 * wing_loading/(ρ * CL_max)

# Drag polar
span_efficiency_factor(e, AR) = 1 / (π * e * AR)
drag_polar(CD0, k, CL) = CD0 + k * CL^2

# Wing loading
wing_loading_stall_speed(V_stall, CL_max, ρ = 1.225) = 0.5 * ρ * V_stall^2 * CL_max

abstract type Aircraft end

"""
Definition for a half-wing consisting of airfoils, span lengths, dihedrals, and sweep angles.
"""
# chords :: Array{Float64} # Chord lengths (m)
struct HalfWing <: Aircraft
    chords :: Array{<: Real}
    spans :: Array{<: Real}  # Leading-edge to leading-edge distance between foils (m)
    dihedrals :: Array{<: Real} # Dihedral angles (deg)
    sweeps :: Array{<: Real} # Leading-edge sweep angles (deg)
    twists :: Array{<: Real} # Twist angles (deg)
    HalfWing(chords, spans, dihedrals, sweeps, twists) = new(chords, spans, deg2rad.(dihedrals), deg2rad.(sweeps), -deg2rad.(twists)) # Convert to radians
end

"""
Computes the aspect ratio given a span and area.
"""
aspect_ratio(span, area) = span^2 / area

"""
Computes the mean geometric chord given a span and area.
"""
mean_geometric_chord(span, area) = area / span

"""
Computes the taper ratio given root and tip chords.
"""
taper_ratio(root_chord, tip_chord) = tip_chord / root_chord

"""
Computes the area given a span and chord.
"""
area(span, chord) = span * chord

"""
Computes the mean aerodynamic chord, essentially a root-mean-square chord, given a root chord and taper ratio.
"""
mean_aerodynamic_chord(root_chord, taper_ratio) = (2/3) * root_chord * (1 + taper_ratio + taper_ratio^2)/(1 + taper_ratio)

"""
Computes the quarter-chord length given a chord.
"""
quarter_chord(chord) = 0.25 * chord

"""
Computes the planform span of a half-wing.
"""
span(wing :: HalfWing) = sum(wing.spans .* cos.(wing.dihedrals) .* cos.(fwdsum(wing.twists) / 2))

"""
Computes the projected area of a half-wing.
"""
function projected_area(wing :: HalfWing)
    mean_chords = fwdsum(wing.chords) / 2 # Mean chord lengths of sections.
    mean_twists = fwdsum(wing.twists) / 2 # Mean twist angles of sections
    sum(@. wing.spans * cos(wing.dihedrals) * cos(wing.sweeps) * mean_chords * cos(mean_twists))
end

"""
Computes the mean aerodynamic chord of a half-wing.
"""
function mean_aerodynamic_chord(wing :: HalfWing)
    mean_chords = fwdsum(wing.chords) / 2
    taper_ratios = fwddiv(wing.chords)
    areas = mean_chords .* wing.spans
    macs = mean_aerodynamic_chord.(wing.chords[1:end-1], taper_ratios)
    sum(macs .* areas) / sum(areas)
end

"""
Computes the leading edge coordinates of a HalfWing, with an option to flip the signs of the y-coordinates.
"""
function lead_wing(wing :: HalfWing, flip :: Bool = false)
    spans, dihedrals, sweeps = wing.spans, wing.dihedrals, wing.sweeps

    sweeped_spans = [ 0; cumsum(spans .* tan.(sweeps)) ]
    dihedraled_spans = [ 0; cumsum(spans .* tan.(dihedrals)) ]
    cum_spans = [ 0; cumsum(spans) ]
    
    leading = SVector.(sweeped_spans, flip ? -cum_spans : cum_spans, dihedraled_spans)

    flip ? leading[end:-1:1] : leading
end

"""
Computes the leading and trailing edge coordinates of a HalfWing, with an option to flip the signs of the y-coordinates.
"""
function wing_bounds(wing :: HalfWing, flip :: Bool = false)
    chords = wing.chords
    twisted_chords = chords .* sin.(wing.twists)
    
    leading = lead_wing(wing, flip)
    trailing = SVector.(chords, (zeros ∘ length)(chords), twisted_chords) 

    leading, (flip ? trailing[end:-1:1] : trailing) .+ leading
end

"""
A composite type consisting of two HalfWings. `left' is meant to have its y-coordinates flipped ∈ Cartesian representation.
"""
struct Wing <: Aircraft
    left :: HalfWing
    right :: HalfWing
end

span(wing :: Wing) = span(wing.left) + span(wing.right)
projected_area(wing :: Wing) = projected_area(wing.left) + projected_area(wing.right)
mean_aerodynamic_chord(wing :: Wing) = (mean_aerodynamic_chord(wing.left) + mean_aerodynamic_chord(wing.right)) / 2
aspect_ratio(wing :: Union{Wing, HalfWing}) = aspect_ratio(span(wing), projected_area(wing))

"""
Computes the wing bounds for a Wing.
"""
function wing_bounds(wing :: Wing)
    left_lead, left_trail = wing_bounds(wing.left, true)
    right_lead, right_trail = wing_bounds(wing.right)

    leading = [ left_lead; right_lead ]
    trailing = [ left_trail; right_trail ]

    leading, trailing
end

end


# # Power loading
# power_loading(tbW, V, η_p) = tbW * V / η_p

# # Thrust-to-weight ratios
# thrust_to_weight_fw_cruise(q, CD0, wing_loading) = q * CD0 * 1 / (wing_loading) + k / q * wing_loading

# best_climb_rate(wing_loading, k, CD0, ρ = 1.225) = sqrt((2 / ρ) * wing_loading * sqrt(k / (3 * CD0))) 
# thrust_to_weight_fw_climb(rate_of_climb, best_climb_rate, q, CD0) = rate_of_climb / best_climb_rate + q / wing_loading * CD0 + k / q * wing_loading

# thrust_to_weight_vtol_climb(wing_loading, V_vtol, S_total, S_wing, CD, ρ = 1.225) = 1.2(1 + ρ * V_vtol^2 * S_total / S_wing / wing_loading)
# # W + 0.5 * ρ * V_vtol^2 * S_proj * CD

# # Mass estimation
# mass_takeoff(m_vtol_prop, m_fixed_prop, m_payload, ff_batt, ff_struct, ff_subsys, ff_avionics) = (m_vtol_prop + m_fixed_prop + m_payload) / (1 - (ff_batt + ff_struct + ff_subsys + ff_avionics))

# flat_plate_cd(α) = 2(sin(α))^3

# figure_of_merit(thrust) = 0.4742 * thrust^0.0793

# hover_speed(T, S_p, ρ) = sqrt(0.5 * T / (ρ * S_p))

# induced_velocity(climb_rate_vtol, v_h) = let x = climb_rate_vtol / (2v_h); -x + sqrt(x^2 + 1) end

# disk_loading(M_TO) = 3.2261M_TO + 74.991

# prop_area(W_TO, DL, η_prop) = W_TO / (DL * η_prop)


# area(a, b, c, d) = a * b * c / d
# prodsratio(a, b, c, d) = a * b / (c * d) 