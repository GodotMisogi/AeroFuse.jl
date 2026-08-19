## Aircraft analysis with the 3D doublet-source (Morino) panel method, coupled to a
## slender-body fuselage line — assembled generically from a `NamedTuple` of components,
## mirroring the vortex-lattice `solve_case` workflow (see `vlm_aircraft.jl`).
using AeroFuse

## Surfaces
#==========================================================================================#

# Wing
wing = Wing(
    foils     = fill(naca4(0, 0, 1, 2), 2),
    chords    = [1.0, 0.6],
    twists    = [2.0, 0.0],
    spans     = [4.0],
    dihedrals = [5.0],
    sweeps    = [5.0],
    symmetry  = true,
)

# Horizontal tail
htail = Wing(
    foils    = fill(naca4(0, 0, 1, 2), 2),
    chords   = [0.7, 0.42],
    spans    = [1.25],
    sweeps   = [6.39],
    position = [4.0, 0, -0.1],
    angle    = -2.0,
    axis     = [0.0, 1, 0],
    symmetry = true,
)

# Vertical tail
vtail = Wing(
    foils    = fill(naca4(0, 0, 0, 9), 2),
    chords   = [0.7, 0.42],
    spans    = [1.0],
    sweeps   = [7.97],
    position = [4.0, 0, 0],
    angle    = 90.0,
    axis     = [1.0, 0, 0],
)

# Fuselage — a slender-body line. NOTE: the line singularities sit on the body axis, so keep
# the fuselage from extending under a lifting surface (here it ends near x = 3.5, ahead of
# the tail at x = 4) to avoid an unphysically large induced cross-flow on an embedded tail.
fuse = HyperEllipseFuselage(
    radius   = 0.5,
    length   = 4.5,
    x_a      = 0.2,
    x_b      = 0.75,
    d_nose   = -0.2,
    position = [-1.0, 0.0, -0.1],
)

## Assemble the aircraft
#==========================================================================================#

# A moderate chordwise count keeps the trailing-edge panels from becoming microscopic on the
# short-chord tails (which would corrupt the finite-difference surface speeds).
n_chord  = 25
aircraft = (
    wing  = surface_panels(WingMesh(wing,  [12], 16), [12], n_chord),
    htail = surface_panels(WingMesh(htail, [8],  12), [8],  n_chord),
    vtail = surface_panels(WingMesh(vtail, [8],  12), [8],  n_chord),
    fuse  = make_fuselage_line(fuse; n = 20),
)

## Case
#==========================================================================================#
fs = Freestream(
    alpha = 3.0,
    beta  = 0.0,
)

ref = References(
    speed    = 1.0,
    area     = projected_area(wing),
    span     = span(wing),
    chord    = mean_aerodynamic_chord(wing),
    location = mean_aerodynamic_center(wing),
)

## Solve
sys = solve_case(aircraft, fs, ref; wake_length = 100.0)

## Aerodynamic coefficients
#==========================================================================================#

# Nearfield (surface-pressure + fuselage Munk) and farfield (wake-circulation) coefficients,
# printed per component and for the whole aircraft.
print_coefficients(sys, :aircraft; components = true)

# Lift by Kutta–Joukowsky from the shed circulation
@show lift_coefficients(sys)   # per surface
@show lift_coefficient(sys)    # total

# Fuselage slender-body cross-flow: near-zero net lift with a nonzero doublet distribution.
if !isnothing(sys.fuse_doublets)
    λ = sys.fuse_doublets
    @info "Fuselage slender-body cross-flow" λ_range = extrema(λ) net_lift = λ[end] - λ[1]
end

## Stability derivatives
#==========================================================================================#

# Force/moment coefficients and their derivatives w.r.t. angle of attack α and sideslip β
# (stability axes). Look for Cm_α < 0 (longitudinal static stability) and Cn_β > 0
# (directional static stability).
dvs = freestream_derivatives(sys; print = true)
