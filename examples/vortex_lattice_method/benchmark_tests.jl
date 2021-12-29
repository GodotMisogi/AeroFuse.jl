## Aircraft analysis benchmarking
using BenchmarkTools
using TimerOutputs

# All subsequent analyses use no symmetry tricks for performance as AeroMDAO hasn't implemented them and apples must be compared to apples.

## AeroMDAO tests: https://github.com/GodotMisogi/AeroMDAO.jl
#=======================================================#
using AeroMDAO

println("AeroMDAO Aircraft Functional -")
function aeromdao_steady_vlm()
    # Wing
    wing = Wing(foils     = Foil.(fill(naca4((0,0,1,2)), 2)),
                chords    = [1.0, 0.6],
                twists    = [2.0, 2.0],
                spans     = [5.0],
                dihedrals = [11.31],
                LE_sweeps = [2.29]);

    # Horizontal tail
    htail = Wing(foils     = Foil.(fill(naca4((0,0,1,2)), 2)),
                 chords    = [0.7, 0.42],
                 twists    = [0.0, 0.0],
                 spans     = [1.25],
                 dihedrals = [0.],
                 LE_sweeps = [6.39],
                 position  = [4., 0, 0],
                 angle     = 0.,
                 axis      = [0., 1., 0.])

    # Vertical tail
    vtail = HalfWing(foils     = Foil.(fill(naca4((0,0,0,9)), 2)),
                     chords    = [0.7, 0.42],
                     twists    = [0.0, 0.0],
                     spans     = [1.0],
                     dihedrals = [0.],
                     LE_sweeps = [7.97],
                     position  = [4., 0, 0],
                     angle     = 90.,
                     axis      = [1., 0., 0.])

    wing_panels, wing_normals   = panel_wing(wing, 20, 10, 
                                             spacing = Cosine()
                                            )

    htail_panels, htail_normals = panel_wing(htail, 12, 6;
                                             spacing  = Cosine()
                                            )

    vtail_panels, vtail_normals = panel_wing(vtail, 10, 5;
                                             spacing  = Cosine()
                                            )

    # Aircraft assembly
    aircraft = ComponentArray(
                              wing  = Horseshoe.(wing_panels,  wing_normals),
                              htail = Horseshoe.(htail_panels, htail_normals),
                              vtail = Horseshoe.(vtail_panels, vtail_normals),
                             );


    # display(size.([ wing_panels[1], htail_panels[1], vtail_panels[1] ])) # Checking sizes

    ρ       = 1.225
    x_ref   = [0.5, 0., 0.]
    V, α, β = 1.0, 5.0, 0.0
    Ω       = [0.0, 0.0, 0.0]
    fs      = AeroMDAO.Freestream(V, α, β, Ω)
    S, b, c = 9.0, 10.0, 0.9
    refs    = References(S, b, c, ρ, x_ref)

    # reset_timer!()
    # @timeit "Solving Case" 
    data = solve_case(aircraft, fs, refs;)

    # @timeit "Computing Forces" Fs = surface_forces(data)
    CFs, CMs = surface_coefficients(data)
    FFs      = farfield_coefficients(data)

    nearfield(data), farfield(data)
    # print_timer()
end

t1 = @benchmark aeromdao_steady_vlm()

## BYU FLOW Lab tests: https://github.com/byuflowlab/VortexLattice.jl
#=======================================================#

# Example taken from: https://flow.byu.edu/VortexLattice.jl/dev/examples/

# Pulled these out to avoid compilation for proper benchmarking
fc   = fill((xc) -> 0, 2) # camberline function for each section
fc_h = fill((xc) -> 0, 2) # camberline function for each section
fc_v = fill((xc) -> 0, 2) # camberline function for each section

##
using VortexLattice

println("BYU FLOW Lab VortexLattice.jl - ")
function byu_steady_vlm()
    # wing
    xle = [0.0, 0.2]
    yle = [0.0, 5.0]
    zle = [0.0, 1.0]
    chord = [1.0, 0.6]
    theta = [2.0*pi/180, 2.0*pi/180]
    phi = [0.0, 0.0]
    ns = 20
    nc = 10
    spacing_s = VortexLattice.Cosine()
    spacing_c = VortexLattice.Cosine()
    mirror = true

    # horizontal stabilizer
    xle_h = [0.0, 0.14]
    yle_h = [0.0, 1.25]
    zle_h = [0.0, 0.0]
    chord_h = [0.7, 0.42]
    theta_h = [0.0, 0.0]
    phi_h = [0.0, 0.0]
    ns_h = 12
    nc_h = 6
    spacing_s_h = VortexLattice.Cosine()
    spacing_c_h = VortexLattice.Cosine()
    mirror_h = true

    # vertical stabilizer
    xle_v = [0.0, 0.14]
    yle_v = [0.0, 0.0]
    zle_v = [0.0, 1.0]
    chord_v = [0.7, 0.42]
    theta_v = [0.0, 0.0]
    phi_v = [0.0, 0.0]
    ns_v = 10
    nc_v = 5
    spacing_s_v = VortexLattice.Cosine()
    spacing_c_v = VortexLattice.Cosine()
    mirror_v = false

    Sref = 9.0
    cref = 0.9
    bref = 10.0
    rref = [0.5, 0.0, 0.0]
    Vinf = 1.0
    ref = Reference(Sref, cref, bref, rref, Vinf)

    alph = 5.0*pi/180
    beta = 0.0
    Omega = [0.0; 0.0; 0.0]
    fs = VortexLattice.Freestream(Vinf, alph, beta, Omega)

    symmetric = [false, false, false]

    # generate surface panels for wing
    wgrid, wing = wing_to_surface_panels(xle, yle, zle, chord, theta, phi, ns, nc;
        mirror=mirror, fc = fc, spacing_s=spacing_s, spacing_c=spacing_c)

    # generate surface panels for horizontal tail
    hgrid, htail = wing_to_surface_panels(xle_h, yle_h, zle_h, chord_h, theta_h, phi_h, ns_h, nc_h;
        mirror=mirror_h, fc=fc_h, spacing_s=spacing_s_h, spacing_c=spacing_c_h)
    VortexLattice.translate!(hgrid, [4.0, 0.0, 0.0]) # Plots also exports translate! -_-
    VortexLattice.translate!(htail, [4.0, 0.0, 0.0]) # Plots also exports translate! -_-

    # generate surface panels for vertical tail
    vgrid, vtail = wing_to_surface_panels(xle_v, yle_v, zle_v, chord_v, theta_v, phi_v, ns_v, nc_v;
        mirror=mirror_v, fc=fc_v, spacing_s=spacing_s_v, spacing_c=spacing_c_v)
    VortexLattice.translate!(vgrid, [4.0, 0.0, 0.0]) # Plots also exports translate! -_-
    VortexLattice.translate!(vtail, [4.0, 0.0, 0.0]) # Plots also exports translate! -_-

    # display(size.([ wing, htail, vtail ])) # Checking sizes

    grids = [wgrid, hgrid, vgrid]
    surfs = [wing, htail, vtail]
    surface_id = [1, 2, 3]

    system = steady_analysis(surfs, ref, fs; symmetric=symmetric, surface_id=surface_id)

    CF, CM = body_forces(system; frame=VortexLattice.Wind())

    CDiff = far_field_drag(system)

    # dCF_b, dCM_b = body_derivatives(system)
    # dCF_s, dCM_s = stability_derivatives(system)

    CF, CM, CDiff #, dCF_b, dCM_b, dCF_s, dCM_s
end

t2 = @benchmark byu_steady_vlm()