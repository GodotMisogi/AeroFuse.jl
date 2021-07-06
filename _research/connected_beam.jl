using Revise
using AeroMDAO
using LinearAlgebra

## Aerodynamic setup
#=========================#

# Define wing
wing = WingSection(root_foil  = naca4((0,0,1,2)),
                   span       = 1.3,
                   dihedral   = 5.0,
                   sweep_LE   = 0.0,
                   taper      = 1.0,
                   root_chord = 0.314,
                   root_twist = 0.0,
                   tip_twist  = 0.0)
wing_mac    = mean_aerodynamic_center(wing)
wing_plan   = plot_wing(wing)
wing_name   = "Wing"
print_info(wing)

# Mesh
span_num        = 2
chord_num       = 1
xyzs            = chord_coordinates(wing, [span_num], chord_num)
panels, normies = panel_wing(wing, span_num, chord_num);
aircraft        = Dict(wing_name => (panels, normies));

# Set up aerodynamic state
aero_state = VLMState(0., 0., 0., [0.0, 0.0, 0.0], 
                      rho_ref   = 1.225,
                      r_ref     = [ wing_mac[1], 0., 0. ],
                      area_ref  = projected_area(wing), 
                      chord_ref = mean_aerodynamic_chord(wing), 
                      span_ref  = span(wing));

# Test case - Fixed speed
aero_state.speed   = 20.
aero_state.alpha   = deg2rad(1.)
aero_state.rho_ref = 0.98

# Build system with initial guess from aerodynamic-only analysis
aero_system, aero_surfs = solve_case(aircraft, aero_state)
print_coefficients(aero_surfs[1], aero_state);

## Load transfer
#=========================#

# Functions on adjacencies
adjacent_joiner(x1, x2) = [ [ x1[1] ]; x1[2:end] .+ x2[1:end-1]; [ x2[end] ] ]
AeroMDAO.VortexLattice.points(horses) = @. r1(bound_leg_center(horses), horses)[:], r2(bound_leg_center(horses), horses)[:]

# Get variables for structural analysis
forces   = surface_forces(aero_surfs[1]) 
horsies  = horseshoes(aero_system)
bounds  = bound_leg_vector.(horsies)[:]
normies = Ref(aero_state.velocity) .× bounds[:]
r1s, r2s = points(horsies)
Ls       = (norm ∘ bound_leg_vector).(horsies)

# Point values
half_forces = forces[:] / 2
M1s         = @. r1s × half_forces
M2s         = @. r2s × half_forces

# Boundary values
pt_forces   = adjacent_joiner(half_forces, half_forces)
pt_moments  = adjacent_joiner(M1s, M2s)

# Transform everything to principal axes 
streams = @. bounds × normies
glob    = repeat([[1; 0; 0], [0; 1; 0], [0; 0; 1]], 1, 3)
dircos  = [ dot.(repeat([c, n, s], 1, 3), glob)  for (c, n, s) in zip(bounds, normies, streams) ]
F_S     = dircos .* forces[:] # NEEDS POINT FORCES HERE

## Boundary condition, setting F = 0 at center of wing
# n = ceil(Int, length(pt_forces) / 2) - 1
# pt_forces[n] .-= pt_forces[n];

# Assembly
Fx = getindex.(pt_forces, 1)
Fy = getindex.(pt_forces, 2)
Fz = getindex.(pt_forces, 3)
Mx = getindex.(pt_moments, 1)
My = getindex.(pt_moments, 2)
Mz = getindex.(pt_moments, 3)

py = (collect ∘ Iterators.flatten ∘ zip)(Fy, My)
pz = (collect ∘ Iterators.flatten ∘ zip)(Fz, Mz)
px = [ Fx[:]; Mx[:] ]

F  = [ py; pz; px ]

## Structural setup
#=========================#

# Beam properties
E     = 70e9  # Elastic modulus, N/m²
G     = 30e9  # Shear modulus, N/m²
σ_max = 200e6 # Yield stress with factor of safety 2.5, N/m²
ρ     = 3e3   # Density, kg/m³
ν     = 0.3   # Poisson's ratio (UNUSED FOR NOW)
R     = 1e-1  # Outer radius, m
t     = 1e-3  # Thickness, m

# Create material and tubes
aluminum = Material(E, G, σ_max, ρ)
tubes    = Tube.(Ref(aluminum), Ls, R, t)

## "FEM" setup
K = tube_stiffness_matrix(aluminum, tubes)

## Solve system
x = K \ F