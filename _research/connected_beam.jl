using Revise
using AeroMDAO
using LinearAlgebra
using StaticArrays

## Aerodynamic setup
#==========================================================================================#

# Define wing
wing = HalfWingSection(root_foil  = naca4((0,0,1,2)),
                   span       = 1.3,
                   dihedral   = 5.0,
                   sweep_LE   = 20.0,
                   taper      = 0.5,
                   root_chord = 0.314,
                   root_twist = 0.0,
                   tip_twist  = 0.0)
wing_mac    = mean_aerodynamic_center(wing)
wing_plan   = plot_wing(wing)
wing_name   = "Wing"
print_info(wing)

# Mesh
span_num        = 12
chord_num       = 1
panels, normies = panel_wing(wing, span_num, chord_num);
aircraft        = Dict(wing_name => (panels, normies));

# Set up aerodynamic state
aero_state = VLMState(1., 0., 0., [0.0, 0.0, 0.0], 
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

## Load transfer scheme
#==========================================================================================#

# Functions on adjacencies
adjacent_joiner(x1, x2) = [ [ x1[1] ]; x1[2:end] .+ x2[1:end-1]; [ x2[end] ] ]
AeroMDAO.VortexLattice.points(horses) = @. r1(bound_leg_center(horses), horses)[:], r2(bound_leg_center(horses), horses)[:]

# Get variables for structural analysis
forces   = surface_forces(aero_surfs[1]) 
horsies  = horseshoes(aero_system)
r1s, r2s = points(horsies)

# Point values
half_forces = forces[:] / 2
M1s         = @. r1s × half_forces
M2s         = @. r2s × half_forces

pt_forces   = adjacent_joiner(half_forces, half_forces)
pt_moments  = adjacent_joiner(M1s, M2s)

# Transform everything to principal axes
wing_pans = (make_panels ∘ coordinates)(wing)[:]
cs        = normalize.(reduce(vcat, diff.(bound_leg.(wing_pans), dims = 1)))
ns        = normalize.(panel_normal.(wing_pans))
ss        = cs .× ns

gloref = repeat([[1; 0; 0], [0; 1; 0], [0; 0; 1]], 1, 3)
dircos = [ dot.(repeat([c, n, s], 1, 3), gloref) for (c, n, s) in zip(cs, ns, ss) ]

zero_vec = [SVector(0,0,0.)]

Fs  = [ [ zero_vec; pt_forces[2:end]  ] ]
Ms  = [ [ zero_vec; pt_moments[2:end] ] ]
F_S = map(x -> [ x[1] * vec for vec in x[2] ], zip(dircos, Fs))
M_S = map(x -> [ x[1] * vec for vec in x[2] ], zip(dircos, Ms))

## Splitting for Wing case
#==========================================================================================#

function middle_index(x :: AbstractArray)
    n = length(x)
    if n % 2 == 0
        Int(n / 2)
    else
        ceil(Int, n / 2)
    end
end

middle(x :: AbstractVector) = @view x[middle_index(x)]

left_forces   = @views [ pt_forces[1:middle_index(pt_forces)-1]; zero_vec ]
right_forces  = @views [ zero_vec; pt_forces[middle_index(pt_forces)+1:end] ]
left_moments  = @view pt_moments[1:middle_index(pt_moments)]
right_moments = @view pt_moments[middle_index(pt_moments):end]

# F_S = map(x -> [ x[1] * vec for vec in x[2] ], zip(dircos, [left_forces,  right_forces ]))
# M_S = map(x -> [ x[1] * vec for vec in x[2] ], zip(dircos, [left_moments, right_moments]))

## Structural setup
#==========================================================================================#

# Beam properties
E     = 70e9  # Elastic modulus, N/m²
G     = 30e9  # Shear modulus, N/m²
σ_max = 200e6 # Yield stress with factor of safety 2.5, N/m²
ρ     = 3e3   # Density, kg/m³
ν     = 0.3   # Poisson's ratio (UNUSED FOR NOW)
R     = 1e-1  # Outer radius, m
t     = 1e-3  # Thickness, m

Ls    = norm.(cs)

# Create material and tubes
aluminum = Material(E, G, σ_max, ρ)
tubes    = Tube.(Ref(aluminum), Ls, R, t)

# Assemble RHS
function assemble_fem_dynamics(pt_forces, pt_moments)
    Fx = getindex.(pt_forces,  1)
    Fy = getindex.(pt_forces,  2)
    Fz = getindex.(pt_forces,  3)
    Mx = getindex.(pt_moments, 1)
    My = getindex.(pt_moments, 2)
    Mz = getindex.(pt_moments, 3)

    py = (collect ∘ Iterators.flatten ∘ zip)(Fy, My)
    pz = (collect ∘ Iterators.flatten ∘ zip)(Fz, Mz)
    px = [ Fx[:]; Mx[:] ]

    F  = [ py; pz; px ]
end

# FEM_RHS = assemble_fem_dynamics(F_S, M_S)
FEM_RHS = reduce(vcat, assemble_fem_dynamics.(F_S, M_S))

## "FEM" setup
K = tube_stiffness_matrix(aluminum, tubes, fill(span_num, length(tubes)))

## Solve system(s)
xs = K \ FEM_RHS


## Plotting
#==========================================================================================#

using Plots
gr(dpi = 300)

plot(aspect_ratio = 1, camera = (45, 45))
plot!.(plot_panels(panels[:]), color = :black, label = :none)
plot!(wing_plan, color = :blue, label = :none)
plot!()

mid_pts = midpoint.(wing_pans)

quiver!(getindex.(mid_pts, 1)[:], getindex.(mid_pts, 2)[:], getindex.(mid_pts, 3)[:], quiver=(getindex.(cs, 1)[:], getindex.(cs, 2)[:], getindex.(cs, 3)[:]), color = :orange)
quiver!(getindex.(mid_pts, 1)[:], getindex.(mid_pts, 2)[:], getindex.(mid_pts, 3)[:], quiver=(getindex.(ns, 1)[:], getindex.(ns, 2)[:], getindex.(ns, 3)[:]), color = :red)
quiver!(getindex.(mid_pts, 1)[:], getindex.(mid_pts, 2)[:], getindex.(mid_pts, 3)[:], quiver=(getindex.(ss, 1)[:], getindex.(ss, 2)[:], getindex.(ss, 3)[:]), color = :brown)