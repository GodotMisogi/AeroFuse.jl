"""
	Wing(left :: HalfWing, right :: HalfWing)

A composite type consisting of two `HalfWing`s with fields `left` and `right` for constructing a wing.
"""
struct Wing{T <: Real} <: Aircraft
	left :: HalfWing{T}
	right :: HalfWing{T}
end

"""
	span(wing :: Wing)
	
Compute the span of a `Wing`.
"""
span(wing :: Wing) = span(wing.left) + span(wing.right)

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

	leading 	= [ left_lead; right_lead ]
	trailing 	= [ left_trail; right_trail ]

	leading, trailing
end

mean_aerodynamic_center(wing :: Wing, factor = 0.25) = (mean_aerodynamic_center(wing.right, factor) .+ mean_aerodynamic_center(wing.left, factor)) ./ 2
