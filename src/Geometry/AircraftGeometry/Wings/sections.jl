## Section definitions
#==========================================================================================#

# mutable struct WingSection{T <: Real, F <: AbstractFoil, A <: AbstractAffineMap, W <: AbstractWing}
#     S :: T
#     AR :: T
#     b  :: T
#     dihedral :: T
#     sweep :: T
#     w_sweep :: T
#     c_r:: T
#     c_t :: T
#     tau_r :: T
#     tau_t :: T
#     foil_r :: F
#     foil_t :: F
#     affine :: A
#     wing :: W
# end

# Span length from aspect ratio and area
span_length(b, AR) = √(b * AR)

# Root chord
root_chord(S, b, λ) = (2 * S) / (b * (1 + λ))

function Wing(span, dihedral, sweep, w_sweep, taper, root_chord, root_twist, tip_twist, root_control, tip_control, root_foil, tip_foil, affine, symmetry, flip)

    # Control surface deflections at root and tip
    root_foil = control_surface(
        root_foil; 
        angle = root_control[1],
        hinge = root_control[2]
    )

    tip_foil  = control_surface(
        tip_foil; 
        angle = tip_control[1],
        hinge = tip_control[2]
    )

    section = Wing(
        foils     = [root_foil, tip_foil],
        chords    = [root_chord, taper * root_chord],
        twists    = [root_twist, tip_twist],
        spans     = [span],
        dihedrals = [dihedral],
        sweeps    = [sweep],
        w_sweep   = w_sweep,
        affine    = affine,
        symmetry  = symmetry,
        flip      = flip
    )

    return section
end

# function WingSection(S, taper, sweep, w_sweep, dihedral, root_twist, tip_twist, root_foil, tip_foil, affine, symmetry)
#     span = S^2 / AR
#     root_chord = 2 * S / (b * (1 + λ))

#     Wing(span, dihedral, sweep, w_sweep, taper, root_chord, root_twist, tip_twist, root_control, tip_control, root_foil, tip_foil, affine, symmetry, flip)
# end

"""
    WingSection(; span, dihedral, sweep, taper, root_chord,
                  root_twist, tip_twist, root_foil, tip_foil,
                  position, angle, axis)

Define a `Wing` in the ``x``-``z`` plane, with optional Boolean arguments for symmetry and flipping in the plane.

# Arguments
- `span :: Real = 1.`: Span length 
- `dihedral :: Real = 1.`: Dihedral angle (degrees)
- `LE_sweep :: Real = 0.`: Leading-edge sweep angle (degrees)
- `taper :: Real = 1.`: Taper ratio of tip to root chord
- `root_chord :: Real = 1.`: Root chord length
- `root_twist :: Real = 0.`: Twist angle at root (degrees)
- `tip_twist :: Real = 0.`: Twist angle at tip (degrees)
- `root_foil :: Foil = naca4((0,0,1,2))`: Foil coordinates at root
- `tip_foil :: Foil = root_foil`: Foil coordinates at tip
- `root_control :: NTuple{2} = (0., 0.75)`: (Angle, chord-length ratio) for adding a control surface at the root.
- `tip_control :: NTuple{2} = root_control`: (Angle, chord-length ratio) for adding a control surface at the tip. Defaults to root control's settings.
- `position :: Vector{Real} = zeros(3)`: Position
- `angle :: Real = 0.`: Angle of rotation (degrees)
- `axis :: Vector{Real} = [0.,1.,0.]`: Axis of rotation
"""
function WingSection(;
        area         = 1.,
        aspect       = 6.,
        dihedral     = 0.,
        sweep        = 0.,
        w_sweep      = 0.,
        taper        = 1.,
        root_twist   = 0.,
        tip_twist    = 0.,
        root_foil    = naca4((0,0,1,2)),
        tip_foil     = root_foil,
        root_control = (0., 0.75),
        tip_control  = root_control,
        position     = zeros(3),
        angle        = 0.,
        axis         = SVector(0., 1., 0.),
        affine       = AffineMap(AngleAxis(deg2rad(angle), axis...), position),
        symmetry     = false,
        flip         = false
    )

    b_w = span_length(area, aspect)
    c_root_w = root_chord(area, b_w, taper)

    if symmetry
        b_w /= 2
    end
    
    return Wing(b_w, dihedral, sweep, w_sweep, taper, c_root_w, root_twist, tip_twist, root_control, tip_control, root_foil, tip_foil, affine, symmetry, flip)
end

## Standard stabilizers (non-exhaustive)
#==========================================================================================#

# HorizontalTail(; root_chord = 1.0, taper = 1.0, span = 1.0, sweep = 0., w_sweep = 0., root_twist = 0., tip_twist = 0., root_foil = naca4(0,0,1,2), tip_foil = root_foil, angle = 0.) = WingSection(root_chord = root_chord, taper = taper, span = span, sweep = sweep, w_sweep = w_sweep, root_twist = root_twist, tip_twist = tip_twist, root_foil = naca4(0,0,1,2), tip_foil = tip_foil, angle = angle, axis = [0., 1., 0.])

# VerticalTail(; root_chord = 1.0, taper = 1.0, span = 1.0, sweep = 0., w_sweep = 0., root_foil = naca4(0,0,1,2), tip_foil = root_foil) = HalfWingSection(root_chord = root_chord, taper = taper, span = span, sweep = sweep, w_sweep = w_sweep, tip_foil = tip_foil, angle = 90., axis = [1., 0., 0.])