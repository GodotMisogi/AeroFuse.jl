## Coupled body-panel fuselage + vortex-lattice aircraft (T-tail)
#
# A skin-panelled `HyperEllipseFuselage` (constant-strength source panels, Neumann
# no-penetration) is solved monolithically with a vortex-lattice wing/tail: the body source
# strengths and the lifting-surface circulations share one AIC. This is the 3-D skin-panel
# counterpart to the slender-body `make_fuselage_line` coupling in `vlm_aircraft.jl`.
#
# The lifting surfaces are mounted at the fuselage SURFACE rather than run through the body
# interior: the wing is split into two exposed half-wings whose roots sit on the fuselage
# flanks (a carry-through gap replaces the buried centre section), the vertical tail's root
# sits on the fuselage top, and the horizontal tail is a T-tail carried on the fin tip, clear
# of the body. This is the physically faithful wing-body layout and keeps the vortex legs out
# of the body interior.
#
# Note: side-force symmetry on a symmetric aircraft (CY = 0 at β = 0) comes from the fuselage
# mesh being mirror-symmetric about the x-z plane, which `make_body_panels` enforces (it rounds
# `n_circ` up to an odd value). Mounting the wing on the skin does not by itself null the
# side force — the root vortex then sits against the fuselage flank panels — so the symmetric
# mesh is what makes the ±y loads cancel exactly.
using AeroFuse
using ComponentArrays

## Compute backend
# `nothing` runs the serial host solver. To parallelise the influence-matrix assembly, linear
# solve and O(N²) post-processing, uncomment one backend (its package must be installed):
backend = nothing
# using KernelAbstractions; backend = CPU()          # Multithreaded host (start Julia with `-t auto`)
# using Metal;              backend = MetalBackend() # Apple GPU (Float32 precision)
# using CUDA;               backend = CUDABackend()  # NVIDIA GPU

## Fuselage (paneled skin)
#=========================================================#
fuse = HyperEllipseFuselage(
    radius   = 0.6,               # Max radius (m)
    length   = 6.0,               # Length (m)
    c_nose   = 2,                 # Nose curvature
    c_rear   = 2,                 # Rear curvature
    position = [-1.4, 0.0, 0.0],  # Nose location, so the wing sits mid-body
)

# Fuselage half-width at the wing station and top height at the tail station (skin mounts).
z_top  = 0.45  # Fuselage top at the empennage station (x ≈ 4)

## Lifting surfaces (mounted on the fuselage surface)
#=========================================================#

# Exposed half-wings: roots on the fuselage flanks at y = ±R_side, with a carry-through gap in
# between (built as two `symmetry = false` half-wings, the port one flipped in the x-z plane).
wing_sections = (
    foils     = fill(naca4(2, 4, 1, 2), 2),
    chords    = [1.0, 0.6],
    twists    = [2.0, 0.0],
    spans     = [3.4],
    dihedrals = [5.0],
    sweeps    = [5.0],
)
x_wing = 1.4
wing_r = Wing(; wing_sections..., symmetry = false, flip = false, position = [x_wing, fuse.radius, 0.])
wing_l = Wing(; wing_sections..., symmetry = false, flip = true,  position = [x_wing, -fuse.radius, 0.])

# Reference (gross, symmetric) wing for the non-dimensionalisation and MAC location.
wing_ref = Wing(; wing_sections..., symmetry = true, position = [0., 0., 0.])

# Vertical tail: root on the fuselage top surface.
vtail = Wing(
    foils    = fill(naca4(0, 0, 0, 9), 2),
    chords   = [0.7, 0.42],
    spans    = [1.0],
    sweeps   = [7.97],
    position = [4.0, 0.0, z_top],  # Root at the fuselage crown
    angle    = 90.0,
    axis     = [1.0, 0.0, 0.0],
)

# Horizontal tail: T-tail carried on the fin tip (well clear of the fuselage).
z_fin_tip = z_top + 1.0  # Fin root z + fin span
htail = Wing(
    foils     = fill(naca4(0, 0, 1, 2), 2),
    chords    = [0.7, 0.42],
    spans     = [1.25],
    sweeps    = [6.39],
    position  = [4.15, 0.0, z_fin_tip],  # At the swept-back fin tip
    angle     = -2,
    axis      = [0.0, 1.0, 0.0],
    symmetry  = true,
)

## Meshing
#=========================================================#
wing_r_mesh = WingMesh(wing_r, [24], 12)
wing_l_mesh = WingMesh(wing_l, [24], 12)
htail_mesh  = WingMesh(htail,  [12], 8)
vtail_mesh  = WingMesh(vtail,  [12], 8)

# Assemble the aircraft with the generic `elements(geometry, model)` interface: each component
# is built by passing its geometry together with an element-model specification. The lifting
# surfaces use `Horseshoe()` (swap for `VortexRing()` to get a vortex-ring lattice); the
# fuselage skin uses `SourcePanel()` (constant-strength source panels, Neumann). The source
# panels couple monolithically with the lifting surfaces via the shared AIC (standard Neumann
# rows, no bespoke solve), and the body block's name carries no behaviour — the pressure force
# model is selected by the element type, not the field name. `n_secs`/`n_circ` set the
# axial/circumferential panel density.
aircraft = ComponentVector(
    wing = hcat(
        elements(wing_l_mesh, Horseshoe()), 
        elements(wing_r_mesh, Horseshoe())
    ),
    htail  = elements(htail_mesh,  Horseshoe()),
    vtail  = elements(vtail_mesh,  Horseshoe()),
    body   = elements(fuse, SourcePanel(); n_secs = 20, n_circ = 24),
)

## Case
#=========================================================#
fs = Freestream(
    alpha = 3.0,
    beta  = 0.0,
    omega = [0.0, 0.0, 0.0],
)

ref = References(
    speed     = 150.0,
    density   = 1.225,
    viscosity = 1.5e-5,
    area      = projected_area(wing_ref),
    span      = span(wing_ref),
    chord     = mean_aerodynamic_chord(wing_ref),
    location  = mean_aerodynamic_center(wing_ref),
)

## Solve
#=========================================================#
@time sys = solve_case(
    aircraft, fs, ref;
    backend          = backend, # Compute backend (see top of file)
    print            = true,   # Print aircraft totals
    print_components = true,    # Print per-component (incl. body) breakdown
)

## Aerodynamic coefficients
#=========================================================#
nf = nearfield(sys)
ff = farfield(sys)

nfs = nearfield_coefficients(sys)  # Per-component nearfield coefficients (incl. `body`)

# Body surface pressures and pressure-integrated force
Cps      = body_pressure_coefficients(sys)
body_CFs = body_forces(sys)
@info "Body Cp range" extrema(Cps)

## Stability derivatives (body is fully differentiable through the coupled solve)
# Derivatives re-solve on the serial host path (ForwardDiff), whatever `backend` is set above.
@time dvs = freestream_derivatives(sys;
    axes             = Stability(),
    print            = true,
    print_components = true,
)

## Plotting — surface pressure contours (body skin C_p + lifting-surface loading)
#=========================================================#
# Renders a 3-D scene colouring the fuselage skin by its surface pressure coefficient and the
# lifting surfaces by their load magnitude (see the file for details). Requires CairoMakie.
include("vlm_body_panel_makie_plot.jl")
