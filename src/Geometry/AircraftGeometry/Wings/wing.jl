# """
#     Wing(; foils, chords, twists, 
#            spans, dihedrals, LE_sweeps,
#            position, angle, axis)

#     Wing(left :: HalfWing, right :: HalfWing)
#     Wing(half :: HalfWing)


# A composite type consisting of fields `left` and `right`, each a `HalfWing` for constructing a wing.

# A single argument generates a symmetric `Wing`, viz. `left === right`.

# # Arguments
# - `foils :: Vector{Foil}`
# - `chords`
# - `twists`
# - `spans`
# - `dihedrals`
# - `LE_sweeps `
# - `position = zeros(3) `
# - `angle    = 0 `
# - `axis     = [0.,1.,0.] `

# """
# struct Wing{M <: AbstractWing} <: AbstractWing
#     left  :: M
#     right :: M
#     # function Wing(left :: T, right :: T) where T <: AbstractWing
#     #     @assert first(left.chords) == first(right.chords) && first(left.twists) == first(right.twists) "Properties of first section must match."
#     #     # Add more conditions...
#     #     Wing{T}(left, right)
#     # end
# end

# # Getters
# left(wing :: Wing)      = wing.left
# right(wing :: Wing)     = wing.right
# foils(wing :: Wing)     = @views [ (foils ∘ left)(wing)[end:-1:2]     ; (foils ∘ right)(wing)     ]
# chords(wing :: Wing)    = @views [ (chords ∘ left)(wing)[end:-1:2]    ; (chords ∘ right)(wing)    ]
# twists(wing :: Wing)    = @views [ (twists ∘ left)(wing)[end:-1:2]    ; (twists ∘ right)(wing)    ]
# spans(wing :: Wing)     = @views [ (reverse ∘ spans ∘ left)(wing)     ; (spans ∘ right)(wing)     ]
# dihedrals(wing :: Wing) = @views [ (reverse ∘ dihedrals ∘ left)(wing) ; (dihedrals ∘ right)(wing) ]

# sweeps(wing :: Wing, w = 0.) = [ reverse(sweeps(left(wing), w)); sweeps(right(wing), w) ]

# Base.position(wing :: Wing) = right(wing).affine.translation
# orientation(wing :: Wing)   = right(wing).affine.linear

# affine_transformation(wing :: Wing) = affine_transformation(right(wing))

# function (f :: AffineMap)(wing :: Wing) 
#     wing_r = @set wing.right.affine = f ∘ wing.right.affine
#     wing_l = @set wing_r.left.affine = f ∘ wing_r.right.affine

#     return wing_l
# end


# # Symmetric wings
# Wing(bing :: HalfWing) = Wing(bing, bing)

# Wing(;
#     chords,
#     twists    = zero(chords), 
#     spans     = ones(length(chords) - 1) / (length(chords) - 1), 
#     dihedrals = zero(spans), 
#     sweeps    = zero(spans),
#     w_sweep   = 0.0,
#     foils     = fill(naca4(0,0,1,2), length(chords)),
#     position  = zeros(3),
#     angle     = 0.,
#     axis      = SVector(0., 1., 0.),
#     affine    = AffineMap(AngleAxis(angle, axis...), position),
# ) = Wing(HalfWing(foils, chords, twists, spans / 2, dihedrals, sweeps, w_sweep, affine))

# span(wing :: Wing) = (span ∘ left)(wing) + (span ∘ right)(wing)

# taper_ratio(wing :: Wing) = (taper_ratio ∘ left)(wing), (taper_ratio ∘ right)(wing)

# projected_area(wing :: Wing) = (projected_area ∘ left)(wing) + (projected_area ∘ right)(wing)

# mean_aerodynamic_chord(wing :: Wing) = ((mean_aerodynamic_chord ∘ left)(wing) + (mean_aerodynamic_chord ∘ right)(wing)) / 2

# """
#     wing_bounds(wing :: Wing)

# Return the leading and trailing edge coordinates of a `Wing`.
# """
# function wing_bounds(wing :: Wing)
#     left_lead, left_trail   = wing_bounds(left(wing), true)
#     right_lead, right_trail = wing_bounds(right(wing))

#     leading  = @views [ left_lead[1:end-1,:] ; right_lead  ]
#     trailing = @views [ left_trail[1:end-1,:]; right_trail ]

#     leading, trailing
# end

# leading_edge(wing :: Wing)  = @views [ leading_edge(left(wing), true)[1:end-1,:] ; leading_edge(right(wing))  ]
# trailing_edge(wing :: Wing) = @views [ trailing_edge(left(wing), true)[1:end-1,:]; trailing_edge(right(wing)) ]

# reflect_xz(vector) = SVector(vector[1], -vector[2], vector[3])
# mean_aerodynamic_center(wing :: Wing, factor = 0.25) = (mean_aerodynamic_center(right(wing), factor) .+ reflect_xz(mean_aerodynamic_center(left(wing), factor))) ./ 2

# maximum_thickness_to_chord(wing :: Wing, num :: Integer) = [ maximum_thickness_to_chord(wing.left, num)[1:end-1]; maximum_thickness_to_chord(wing.right, num)[1:end] ]