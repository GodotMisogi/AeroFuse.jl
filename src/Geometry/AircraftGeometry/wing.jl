"""
    Wing(left :: HalfWing, right :: HalfWing)

A composite type consisting of two `HalfWing`s with fields `left` and `right` for constructing a wing.
"""
struct Wing{T <: Real} <: Aircraft
    left :: HalfWing{T}
    right :: HalfWing{T}
end

# function Wing(left :: HalfWing{T}, right :: HalfWing{T}) where T <: Real
#     @assert first(left.chords) == first(right.chords) "Properties of first section must match."
#     # Add more conditions...
#     Wing{T}(left, right)
# end

left(wing :: Wing)      = wing.left
right(wing :: Wing)     = wing.right
foils(wing :: Wing)     = @views [ (foils ∘ left)(wing)[end:-1:2]     ; (foils ∘ right)(wing)     ]
chords(wing :: Wing)    = @views [ (chords ∘ left)(wing)[end:-1:2]    ; (chords ∘ right)(wing)    ]
twists(wing :: Wing)    = @views [ (twists ∘ left)(wing)[end:-1:2]    ; (twists ∘ right)(wing)    ]
spans(wing :: Wing)     = @views [ (reverse ∘ spans ∘ left)(wing)     ; (spans ∘ right)(wing)     ]
dihedrals(wing :: Wing) = @views [ (reverse ∘ dihedrals ∘ left)(wing) ; (dihedrals ∘ right)(wing) ]
sweeps(wing :: Wing)    = @views [ (reverse ∘ sweeps ∘ left)(wing)    ; (sweeps ∘ right)(wing)    ]

Base.position(wing :: Wing)     = (position ∘ right)(wing)
orientation(wing :: Wing)  = (orientation ∘ right)(wing)

affine_transformation(wing :: Wing) = affine_transformation(right(wing))

# Symmetric wing
Wing(; chords, twists, spans, dihedrals, LE_sweeps, foils = fill(Foil(naca4((0,0,1,2))), length(chords)), position = zeros(3), angle = 0., axis = [1.,0.,0.]) = let w = HalfWing(foils = foils, chords = chords, twists = twists, spans = spans, dihedrals = dihedrals, LE_sweeps = LE_sweeps, position = position, angle = angle, axis = axis); Wing(w, w) end

# Single section for convenience
WingSection(; span = 1., dihedral = 0., sweep_LE = 0., taper = 1., root_chord = 1., root_twist = 0., tip_twist = 0., root_foil = naca4((0,0,1,2)), tip_foil = naca4((0,0,1,2)), position = zeros(3), angle = 0., axis = [1.,0.,0.]) = let w = HalfWingSection(span = span / 2, dihedral = dihedral, sweep_LE = sweep_LE, taper = taper, root_chord = root_chord, root_twist = root_twist, tip_twist = tip_twist, root_foil = root_foil, tip_foil = tip_foil, position = position, angle = angle, axis = axis); Wing(w, w) end

"""
    span(wing :: Wing)

Compute the span of a `Wing`.
"""
span(wing :: Wing) = (span ∘ left)(wing) + (span ∘ right)(wing)

taper_ratio(wing :: Wing) = (taper_ratio ∘ left)(wing), (taper_ratio ∘ right)(wing)

"""
    projected_area(wing :: Wing)

Compute the projected_area of a `Wing`.
"""
projected_area(wing :: Wing) = (projected_area ∘ left)(wing) + (projected_area ∘ right)(wing)

"""
    mean_aerodynamic_chord(wing :: Wing)

Compute the mean aerodynamic chord of a `Wing`.
"""
mean_aerodynamic_chord(wing :: Wing) = ((mean_aerodynamic_chord ∘ left)(wing) + (mean_aerodynamic_chord ∘ right)(wing)) / 2

"""
    wing_bounds(wing :: Wing)

Return the leading and trailing edge coordinates of a `Wing`.
"""
function wing_bounds(wing :: Wing)
    left_lead, left_trail   = wing_bounds(left(wing), true)
    right_lead, right_trail = wing_bounds(right(wing))

    leading  = @views [ left_lead[1:end-1,:] ; right_lead  ]
    trailing = @views [ left_trail[1:end-1,:]; right_trail ]

    leading, trailing
end

leading_edge(wing :: Wing)  = @views [ leading_edge(left(wing), true)[1:end-1,:] ; leading_edge(right(wing))  ]
trailing_edge(wing :: Wing) = @views [ trailing_edge(left(wing), true)[1:end-1,:]; trailing_edge(right(wing)) ]

reflect_xz(vector) = SVector(vector[1], -vector[2], vector[3])
mean_aerodynamic_center(wing :: Wing, factor = 0.25) = (mean_aerodynamic_center(right(wing), factor) .+ reflect_xz(mean_aerodynamic_center(left(wing), factor))) ./ 2