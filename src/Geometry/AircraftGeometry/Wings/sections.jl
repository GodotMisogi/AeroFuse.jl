## Section definitions
#==========================================================================================#

function WingSection(span, dihedral, sweep, w_sweep, taper, root_chord, root_twist, tip_twist, root_control, tip_control, root_foil, tip_foil, affine, symmetry, flip)

    # Control surface deflections at root and tip
    root_foil = control_surface(root_foil; 
                                angle = root_control[1],
                                hinge = root_control[2])

    tip_foil  = control_surface(tip_foil; 
                                angle = tip_control[1],
                                hinge = tip_control[2])

    section = Wing(
        foils     = [root_foil, tip_foil],
        chords    = [root_chord, taper * root_chord],
        twists    = [root_twist, tip_twist],
        spans     = [span],
        dihedrals = [dihedral],
        sweeps    = [sweep],
        w_sweep   = w_sweep,
        affine    = affine,
        symmetry = symmetry,
        flip      = flip
    )

    return section
end

function WingSection(S, taper, sweep, w_sweep, dihedral, root_twist, tip_twist, root_foil, tip_foil, affine, symmetry)
    span = S^2 / AR
    root_chord = 2 * S / (b * (1 + λ))

    WingSection(span, dihedral, sweep, w_sweep, taper, root_chord, root_twist, tip_twist , root_control, tip_control, root_foil, tip_foil, affine, symmetry, flip)
end

"""
    WingSection(; span, dihedral, sweep, taper, root_chord,
                  root_twist, tip_twist, root_foil, tip_foil,
                  position, angle, axis)

Define a `Wing` consisting of two trapezoidal sections with reflection symmetry in the ``x``-``z`` plane.

# Arguments
- `span       :: Real         = 1.`: Span length 
- `dihedral   :: Real         = 1.`: Dihedral angle (degrees)
- `LE_sweep   :: Real         = 0.`: Leading-edge sweep angle (degrees)
- `taper      :: Real         = 1.`: Taper ratio of tip to root chord
- `root_chord :: Real         = 1.`: Root chord length
- `root_twist :: Real         = 0.`: Twist angle at root (degrees)
- `tip_twist  :: Real         = 0.`: Twist angle at tip (degrees)
- `root_foil  :: Foil         = naca4((0,0,1,2))`: Foil coordinates at root
- `tip_foil   :: Foil         = root_foil`: Foil coordinates at tip
- `position   :: Vector{Real} = zeros(3)`: Position
- `angle      :: Real         = 0.`: Angle of rotation (degrees)
- `axis       :: Vector{Real} = [0.,1.,0.]`: Axis of rotation
"""
function WingSection(;
        span         = 1.,
        dihedral     = 0.,
        sweep        = 0.,
        w_sweep      = 0.,
        taper        = 1.,
        root_chord   = 1.,
        root_twist   = 0.,
        tip_twist    = 0.,
        root_control = (0., 0.75),
        tip_control  = root_control,
        root_foil    = naca4((0,0,1,2)),
        tip_foil     = root_foil,
        position     = zeros(3),
        angle        = 0.,
        axis         = SVector(0., 1., 0.),
        affine       = AffineMap(AngleAxis(angle, axis...), position),
        symmetry     = false,
        flip         = false
    )
    
    return WingSection(span / 2, dihedral, sweep, w_sweep, taper, root_chord, root_twist, tip_twist, root_control, tip_control, root_foil, tip_foil, affine, symmetry, flip)
end

## Standard stabilizers (non-exhaustive)
#==========================================================================================#

# HorizontalTail(; root_chord = 1.0, taper = 1.0, span = 1.0, sweep = 0., w_sweep = 0., root_twist = 0., tip_twist = 0., root_foil = naca4(0,0,1,2), tip_foil = root_foil, angle = 0.) = WingSection(root_chord = root_chord, taper = taper, span = span, sweep = sweep, w_sweep = w_sweep, root_twist = root_twist, tip_twist = tip_twist, root_foil = naca4(0,0,1,2), tip_foil = tip_foil, angle = angle, axis = [0., 1., 0.])

# VerticalTail(; root_chord = 1.0, taper = 1.0, span = 1.0, sweep = 0., w_sweep = 0., root_foil = naca4(0,0,1,2), tip_foil = root_foil) = HalfWingSection(root_chord = root_chord, taper = taper, span = span, sweep = sweep, w_sweep = w_sweep, tip_foil = tip_foil, angle = 90., axis = [1., 0., 0.])