## Blown lift: a propeller slipstream washing over a wing
#
# Demonstrates the prescribed actuator-disk slipstream (`PropellerDisk`) coupled to the vortex
# lattice method. The slipstream enters both the boundary condition (raising the circulation)
# and the local Kutta–Joukowsky velocity (the dynamic-pressure boost inside the tube), so a
# powered wing shows augmented lift over the unpowered baseline.
using AeroFuse
using StaticArrays
using LinearAlgebra # For `normalize` in the flap deflection
using Accessors # For `@set` / `setproperties` in the powered-flap section

## Wing
wing = Wing(
    foils     = fill(naca4(2, 4, 1, 2), 2),
    chords    = [1.0, 0.6],
    twists    = [0.0, 0.0],
    spans     = [4.0],
    dihedrals = [5.0],
    sweeps    = [5.0],
    symmetry  = true,
)

wing_mesh  = WingMesh(wing, [35], 12)
wing_rings = make_vortex_rings(wing_mesh)
aircraft   = ComponentVector(wing = wing_rings)

## Freestream and references
fs = Freestream(alpha = 3.0, beta = 0.0, omega = [0.0, 0.0, 0.0])

ref = References(
    speed     = 50.0,
    density   = 1.225,
    viscosity = 1.5e-5,
    area      = projected_area(wing),
    span      = span(wing),
    chord     = mean_aerodynamic_chord(wing),
    location  = mean_aerodynamic_center(wing),
)

## Baseline (unpowered)
sys_baseline = VortexLatticeSystem(aircraft, fs, ref)
nf_baseline  = nearfield(sys_baseline)

## Propellers: a tractor disk ahead of each half-wing, thrust along the freestream (+x).
# `torque` sets the swirl and `sense = ±1` its direction; mirror the sense across the two
# props so their swirl-induced rolling moments cancel.
prop_right = PropellerDisk(center = SVector(-1.0,  1.2, 0.0), axis = SVector(1.0, 0.0, 0.0),
                           radius = 1.5, thrust = 8000.0, torque = 120.0, sense =  1.0)
prop_left  = PropellerDisk(center = SVector(-1.0, -1.2, 0.0), axis = SVector(1.0, 0.0, 0.0),
                           radius = 1.5, thrust = 8000.0, torque = 120.0, sense = -1.0)

## Blown analysis
sys_blown = VortexLatticeSystem(aircraft, fs, ref; slipstream = [prop_right, prop_left])
nf_blown  = nearfield(sys_blown)

## Compare
println("Unpowered : CL = ", round(nf_baseline.CZ, digits = 4), "  Cl = ", round(nf_baseline.Cl, digits = 6))
println("Blown     : CL = ", round(nf_blown.CZ,    digits = 4), "  Cl = ", round(nf_blown.Cl,    digits = 6))
println("Lift augmentation ΔCL = ", round(nf_blown.CZ - nf_baseline.CZ, digits = 4))

## Powered flap (deflected slipstream / jet flap)
#
# A plain wing in an axial slipstream gains lift only linearly with the velocity ratio. The
# large powered-lift gains come from a flap turning the high-momentum jet. `auto_turn` derives
# the jet turning from the immersed wing panels' trailing-edge deflection, so it need not be
# prescribed.

# Deflect trailing-edge flaps over the full span (symmetric → no roll) using Drela's
# small-angle normal update.
δ_flap  = deg2rad(20.0)
nc, ns  = size(wing_rings)
is_flap = zeros(Bool, nc, ns)
is_flap[end-1:end, :] .= true # Rear two chordwise panel rows
flapped = map(wing_rings, is_flap) do ring, is_cs
    is_cs ? (@set ring.normal = normalize(ring.normal + δ_flap * cross(SVector(0.0, 1.0, 0.0), ring.normal))) : ring
end
aircraft_flap = ComponentVector(wing = flapped)

# Same props, but now turning the jet — derived automatically from the immersed flap panels.
turn_props = [ auto_turn(p, flapped, ref) for p in (prop_right, prop_left) ]

CL_flap       = nearfield(VortexLatticeSystem(aircraft_flap, fs, ref)).CZ
CL_flap_axial = nearfield(VortexLatticeSystem(aircraft_flap, fs, ref; slipstream = [prop_right, prop_left])).CZ
CL_flap_jet   = nearfield(VortexLatticeSystem(aircraft_flap, fs, ref; slipstream = turn_props)).CZ

println()
println("Flap 20°, unpowered   : CL = ", round(CL_flap,       digits = 4))
println("Flap 20°, axial jet   : CL = ", round(CL_flap_axial, digits = 4), "  (ΔCL = ", round(CL_flap_axial - CL_flap, digits = 4), ")")
println("Flap 20°, turned jet  : CL = ", round(CL_flap_jet,   digits = 4), "  (ΔCL = ", round(CL_flap_jet   - CL_flap, digits = 4), ")")
