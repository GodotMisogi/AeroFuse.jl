## Section definitions
#==========================================================================================#

# Span length from aspect ratio and area
span_length(b, AR) = √(b * AR)

# Root chord
root_chord(S, b, λ) = (2 * S) / (b * (1 + λ))

# Single-section constructor
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

"""
    WingSection(; 
        area, aspect, taper
        dihedral, sweep, w_sweep,
        root_twist, tip_twist,
        position, angle, axis,
        symmetry, flip,
        root_foil, tip_foil,
        root_control, tip_control,
    )

Define a `Wing` in the ``x``-``z`` plane, with optional Boolean arguments for symmetry and flipping in the plane.

# Arguments
- `area :: Real = 1.`: Area (m²)
- `aspect :: Real = 6.`: Aspect ratio
- `taper :: Real = 1.`: Taper ratio of tip to root chord
- `dihedral :: Real = 1.`: Dihedral angle (degrees)
- `sweep :: Real = 0.`: Sweep angle (degrees)
- `w_sweep :: Real = 0.`: Chord ratio for sweep angle 
                          e.g., 0    = Leading-edge sweep, 
                                1    = Trailing-edge sweep,
                                0.25 = Quarter-chord sweep
- `root_twist :: Real = 0.`: Twist angle at root (degrees)
- `tip_twist :: Real = 0.`: Twist angle at tip (degrees)
- `root_foil :: Foil = naca4((0,0,1,2))`: `Foil` at root
- `tip_foil :: Foil = root_foil`: `Foil` at tip. Defaults to root foil.
- `root_control :: NTuple{2} = (0., 0.75)`: (Angle, hinge ratio) for adding a control surface at the root.
- `tip_control :: NTuple{2} = root_control`: (Angle, hinge ratio) for adding a control surface at the tip. Defaults to root control's settings.
- `symmetry :: Bool = false`: Symmetric in the ``x``-``z`` plane
- `flip :: Bool = false`: Flip coordinates in the ``x``-``z`` plane
- `position :: Vector{Real} = zeros(3)`: Position (m)
- `angle :: Real = 0.`: Angle of rotation (degrees)
- `axis :: Vector{Real} = [0.,1.,0.]`: Axis of rotation
- `affine :: AffineMap = AffineMap(AngleAxis(deg2rad(angle), axis...), position)`: Affine mapping for the position and orientation via `CoordinateTransformations.jl` (overrides `angle` and `axis` if specified)
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
        axis         = [0., 1., 0.],
        affine       = AffineMap(AngleAxis(deg2rad(angle), axis...), SVector(position...)),
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