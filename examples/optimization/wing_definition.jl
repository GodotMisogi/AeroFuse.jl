# Geometry
function make_wing(xc, w = 0.25)
    # Planform geometry
    n_vars = length(xc)
    bs = fill(1., n_vars - 1)
    wing = Wing(
        foils     = fill(naca4(4,4,1,2), n_vars),
        chords    = xc, # Design variables
        spans     = bs / sum(bs), # Normalizing halfspan to 1
        w_sweep   = w, # Quarter-chord
        symmetry  = true
    )

    # Meshing
    wing_mesh = WingMesh(
        wing, fill(2, n_vars - 1), 6, 
        span_spacing = Uniform(),
    );

    return wing_mesh
end

# VLM analysis
function make_case(α, wing_mesh, refs)
    aircraft = ComponentVector(wing = make_horseshoes(wing_mesh))

    # Freestream conditions
    fs = Freestream(alpha = α)  # Design variable: Angle of attack

    # Solve system
    return VortexLatticeSystem(aircraft, fs, refs)
end

# Aerodynamic forces
function get_forces(system, wing_mesh)
    # Evaluate aerodynamic coefficients
    CDi, CY, CL, Cl, Cm, Cn = nearfield(system)
    # CDi, _, _ = farfield(system)

    # Calculate equivalent flat-plate skin-friction drag
    # CDv = profile_drag_coefficient(wing_mesh, 1.0, system.reference)

    # Calculate local-dissipation/local-friction drag
    CVs = norm.(surface_velocities(system)).wing
    CDv = profile_drag_coefficient(wing_mesh, 0.98, CVs, system.reference)

    return (CDi = CDi, CDv = CDv, CD = CDi + CDv, CL = CL)
end