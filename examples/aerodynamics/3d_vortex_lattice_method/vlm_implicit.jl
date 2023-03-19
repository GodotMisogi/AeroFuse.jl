## VLM implicit analysis and derivatives
using AeroFuse
using LinearAlgebra
using ImplicitAD
using NLsolve
using StructArrays

## Surfaces
#================================================#

function make_surfaces(x)
    ## Wing
    wing = Wing(
        foils     = fill(naca4(2,4,1,2), 2),
        chords    = [1.0, 0.6],
        twists    = [2.0, 0.0],
        spans     = [4.0],
        dihedrals = [5.],
        sweeps    = [5.],
        symmetry  = true
    )

    ## Horizontal tail
    htail = Wing(
        foils     = fill(naca4(0,0,1,2), 2),
        chords    = [0.7, 0.42],
        twists    = [0.0, 0.0],
        spans     = [1.25],
        dihedrals = [0.],
        sweeps    = [6.39],
        position  = [4., 0, -0.1],
        angle     = -2.,
        axis      = [0., 1., 0.],
        symmetry  = true
    )

    ## Vertical tail
    vtail = Wing(
        foils     = fill(naca4(0,0,0,9), 2),
        chords    = [0.7, 0.42],
        twists    = [0.0, 0.0],
        spans     = [1.0],
        dihedrals = [0.],
        sweeps    = [7.97],
        position  = [4., 0, 0],
        angle     = 90.,
        axis      = [1., 0., 0.]
    )

    ## Meshing
    #================================================#

    wing_mesh  = WingMesh(wing, [24], 12, span_spacing = Cosine())
    htail_mesh = WingMesh(htail, [12], 6, span_spacing = Cosine())
    vtail_mesh = WingMesh(vtail, [12], 6, span_spacing = Cosine())

    # Assemble
    aircraft = ComponentArray(
        wing  = make_horseshoes(wing_mesh),
        htail = make_horseshoes(htail_mesh),
        vtail = make_horseshoes(vtail_mesh)
    )

    return aircraft
end

ac = make_surfaces([1.0, 0.6])

## VLM analysis
#================================================#

fs  = Freestream(
    alpha = 0.0, 
    beta  = 0.0, 
    omega = [0., 0., 0.]
)

ref = References(
    speed     = 150.0, 
    density   = 1.225,
    viscosity = 1.5e-5,
)

function solve(system)
    rwrap(R, p) = AeroFuse.VortexLattice.residual!(R, p, system)
    res = nlsolve(
            rwrap,
            zero(system.circulations),
            method = :newton,
            autodiff = :forward,
            show_trace = true
        )
    return res.zero
end

##
function program(x)
    ac = make_surfaces(x)
    fs  = Freestream(
        alpha = convert(eltype(x), 5.0), 
        beta  = 0.0, 
        omega = [0., 0., 0.]
    )
    
    ref = References(
        speed     = 150.0, 
        density   = 1.225,
        viscosity = 1.5e-5,
    )

    sys = VortexLatticeSystem(ac, fs, ref)

    solve(sys)
    # [ velocity(fs) sum(ac.rc) ]
end

##
using ForwardDiff

ForwardDiff.jacobian(program, [1.1, 0.6])
