"""
	Wing(left :: HalfWing, right :: HalfWing)

A composite type consisting of two `HalfWing`s with fields `left` and `right` for constructing a wing.
"""
struct Wing{T <: Real} <: Aircraft
	left :: HalfWing{T}
	right :: HalfWing{T}
end

# Symmetric wing
Wing(; chords :: Vector{T}, twists :: Vector{T}, spans :: Vector{T}, dihedrals :: Vector{T}, sweep_LEs :: Vector{T}, foils :: Vector{Foil{T}} = fill(Foil(naca4((0,0,1,2))), length(chords))) where T <: Real = let w = HalfWing(foils = foils, chords = chords, twists = twists, spans = spans, dihedrals = dihedrals, sweep_LEs = sweep_LEs); Wing(w, w) end

# Single section for convenience
WingSection(; span = 1., dihedral = 0., sweep_LE = 0., taper = 1., root_chord = 1., root_twist = 0., tip_twist = 0., root_foil = naca4((0,0,1,2)), tip_foil = naca4((0,0,1,2))) = let w = HalfWingSection(span = span, dihedral = dihedral, sweep_LE = sweep_LE, taper = taper, root_chord = root_chord, root_twist = root_twist, tip_twist = tip_twist, root_foil = root_foil, tip_foil = tip_foil); Wing(w, w) end

"""
	span(wing :: Wing)
	
Compute the span of a `Wing`.
"""
span(wing :: Wing) = span(wing.left) + span(wing.right)

taper_ratio(wing :: Wing) = taper_ratio(wing.left), taper_ratio(wing.right)

"""
	projected_area(wing :: Wing)
	
Compute the projected_area of a `Wing`.
"""
projected_area(wing :: Wing) = projected_area(wing.left) + projected_area(wing.right)

"""
	mean_aerodynamic_chord(wing :: Wing)
	
Compute the mean aerodynamic chord of a `Wing`.
"""
mean_aerodynamic_chord(wing :: Wing) = (mean_aerodynamic_chord(wing.left) + mean_aerodynamic_chord(wing.right)) / 2

"""
	wing_bounds(wing :: Wing)

Return the leading and trailing edge coordinates of a `Wing`.
"""
function wing_bounds(wing :: Wing)
	left_lead, left_trail 	= wing_bounds(wing.left, true)
	right_lead, right_trail = wing_bounds(wing.right)

	leading  = [ left_lead ; right_lead  ]
	trailing = [ left_trail; right_trail ]

	leading, trailing
end

mean_aerodynamic_center(wing :: Wing, factor = 0.25) = (mean_aerodynamic_center(wing.right, factor) .+ mean_aerodynamic_center(wing.left, factor)) ./ 2