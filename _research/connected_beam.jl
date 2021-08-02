using Revise
using AeroMDAO
using LinearAlgebra
using StaticArrays
using DataFrames

## Aerodynamic setup
#==========================================================================================#

# Define wing
wing = WingSection(root_foil  = naca4(0,0,1,2),
                   span       = 3.11,
                   dihedral   = 5.0,
                   sweep_LE   = 20.0,
                   taper      = 0.5,
                   root_chord = 0.3,
                   root_twist = 2.0,
                   tip_twist  = 0.0)
wing_mac    = mean_aerodynamic_center(wing)
wing_plan   = plot_wing(wing)
wing_name   = "Wing"
print_info(wing, wing_name)

# Mesh
span_num        = 5
chord_num       = 2
# xyzs            = chop_wing(coordinates(wing), [span_num], chord_num)
panels, normies = panel_wing(wing, span_num, chord_num, spacing = "sine");
aircraft        = Dict(wing_name => (panels, normies));

# Set up aerodynamic state
aero_state = VLMState(1., 0., 0., [0.0, 0.0, 0.0], 
                      rho_ref   = 1.225,
                      r_ref     = [ wing_mac[1], 0., 0. ],
                      area_ref  = projected_area(wing), 
                      chord_ref = mean_aerodynamic_chord(wing), 
                      span_ref  = span(wing));

# Test case - Fixed speed
aero_state.speed   = 22.876
aero_state.alpha   = deg2rad(5.)
aero_state.rho_ref = 0.770816

# Build system with initial guess from aerodynamic-only analysis
aero_system = solve_case(aircraft, aero_state)
aero_surfs  = surfaces(aero_system)
print_coefficients(aero_surfs[1], aero_state);

## Load transfer scheme
#==========================================================================================#

# Functions on adjacencies
adjacent_joiner(x1, x2) = [ [ x1[1] ]; x1[2:end,:] .+ x2[1:end-1,:]; [ x2[end] ] ]
# AeroMDAO.VortexLattice.points(horses) = @. r1(bound_leg_center(horses), horses)[:], r2(bound_leg_center(horses), horses)[:]

# Get variables for structural analysis
forces   = surface_forces(aero_surfs[1]) 
horsies  = horseshoes(aero_surfs[1])
# r1s, r2s = points(horsies)

fem_w       = 0.35
plan_coords = coordinates(wing, span_num, chord_num)
mesh        = reshape(permutedims(reduce(hcat, plan_coords)), (size(plan_coords)..., 3))

# Aerodynamic center locations and forces
half_forces = sum(forces, dims = 1)[:] / 2
ac_mine  = bound_leg_center.(horsies)
ac_pts   = reshape(permutedims(reduce(hcat, ac_mine)), (size(ac_mine)..., 3))

# FEM beam node locations as Matrix
fem_arr  = (1 - fem_w) * mesh[1,:,:] + fem_w * mesh[end,:,:]
M1s = sum(let xs = (ac_pts[i,:,:] - fem_arr[1:end-1,:]); 
              @. SVector(xs[:,1], xs[:,2], xs[:,3]) × forces[i,:] / 2 end
              for i in eachindex(ac_pts[:,1,1]))
M2s = sum(let xs = (ac_pts[i,:,:] - fem_arr[2:end,:]); 
              @. SVector(xs[:,1], xs[:,2], xs[:,3]) × forces[i,:] / 2 end
              for i in eachindex(ac_pts[:,1,1]))

## FEM beam node locations as Vector{SVector{3}}
# fem_pts  = (1 - fem_w) * plan_coords[1,:] + fem_w * plan_coords[end,:]
# r1s, r2s = fem_pts[1:end-1], fem_pts[2:end]
# M1s         = @. r1s × half_forces
# M2s         = @. r2s × half_forces

## Need to check if adjacent_joiner is the right definition...
pt_forces   = adjacent_joiner(half_forces, half_forces)
pt_moments  = adjacent_joiner(M1s, M2s)

# Matrix version
pt_loads = zeros(length(half_forces) + 1, 6)
test_force = permutedims(reduce(hcat, half_forces))
pt_loads[1:end-1,1:3] = test_force
pt_loads[2:end,1:3]  += test_force

test_M1s = permutedims(reduce(hcat, M1s))
test_M2s = permutedims(reduce(hcat, M2s))
pt_loads[1:end-1,4:end] = test_M1s
pt_loads[2:end,4:end]  += test_M2s
pt_loads[end:-1:6,:]

## Transform everything to principal axes
wing_pans = (make_panels ∘ coordinates)(wing)[:]
cs        = diff.(bound_leg.(wing_pans), dims = 1)
n_cs      = normalize.(reduce(vcat, cs))
ns        = normalize.(panel_normal.(wing_pans))
ss        = n_cs .× ns

gloref = repeat([[1.; 0; 0], [0.; 1; 0], [0.; 0; 1]], 1, 3)
dircos = [ dot.(repeat([s, c, n], 1, 3)', gloref) for (c, n, s) in zip(n_cs, ns, ss) ]

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

zero_vec = [SVector(0,0,0.)]

# Ls    = (norm ∘ bound_leg_vector).(horsies)
# Fs  = [ [ zero_vec; pt_forces[2:end]  ] ]
# Ms  = [ [ zero_vec; pt_moments[2:end] ] ]
# F_S = map(x -> [ x[1] * vec for vec in x[2] ], zip(dircos, Fs))
# M_S = map(x -> [ x[1] * vec for vec in x[2] ], zip(dircos, Ms))

Ls    = norm.(diff(fem_pts))

# NEED TO CHECK THE DIRECTIONS
left_forces   = @view pt_forces[1:middle_index(pt_forces)-1]
left_moments  = @view pt_moments[1:middle_index(pt_moments)-1]
right_forces  = @views [ zero_vec; pt_forces[middle_index(pt_forces)+1:end] ]
right_moments = @views [ zero_vec; pt_moments[middle_index(pt_moments)+1:end] ]

F_S = map(x -> [ x[1] * vec for vec in x[2] ], zip(dircos, [left_forces,  right_forces ]))
M_S = map(x -> [ x[1] * vec for vec in x[2] ], zip(dircos, [left_moments, right_moments]))

## Structural setup
#==========================================================================================#

# Beam properties
E     = 85e9  # Elastic modulus, N/m²
G     = 25e9  # Shear modulus, N/m²
σ_max = 350e6 # Yield stress with factor of safety 2.5, N/m²
ρ     = 1.6e3   # Density, kg/m³
ν     = 0.3   # Poisson's ratio (UNUSED FOR NOW)
R     = 1e-2  # Outer radius, m
t     = 8e-3  # Thickness, m

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

Fx = getindex.(pt_forces,  1)
Fy = getindex.(pt_forces,  2) 
Fz = getindex.(pt_forces,  3)
Mx = getindex.(pt_moments, 1)
My = getindex.(pt_moments, 2)
Mz = getindex.(pt_moments, 3)

df = DataFrame([ Fx Fy Fz Mx My Mz ], :auto)
rename!(df, [:Fx, :Fy, :Fz, :Mx, :My, :Mz])

py = (collect ∘ Iterators.flatten ∘ zip)(Fy, My)
pz = (collect ∘ Iterators.flatten ∘ zip)(Fz, Mz)
px = [ Fx[:]; Mx[:] ]

## "FEM" setup
K = tube_stiffness_matrix(aluminum, tubes)
F = reduce(vcat, assemble_fem_dynamics.(F_S, M_S))

## Solve system(s)
xs = K \ F

## Displacement transfer scheme
#==========================================================================================#

function compute_displacements(Δs, n)
    # Assembly
    n1 = 2n                 # dim(Fy + My)
    n2 = n1 + 2n            # dim(Fz + Mz)
    n3 = n2 +  n            # dim(Fx)
    n4 = n3 +  n            # dim(Mx)

    dy = @view Δs[1:2:n1]
    θy = @view Δs[2:2:n1]
    dz = @view Δs[n1+1:2:n2]
    θz = @view Δs[n1+2:2:n2]
    dx = @view Δs[n2+1:n3]
    θx = @view Δs[n3+1:n4]

    dx, θx, dy, θy, dz, θz
end

dx, θx, dy, θy, dz, θz = compute_displacements(xs, length(fem_pts))

ds = SVector.(dx, dy, dz)
θs = SVector.(θx, θy, θz)

## Generate DataFrames
df = DataFrame([ dx θx dy θy dz θz ], :auto)
rename!(df, [:dx, :θx, :dy, :θy, :dz, :θz])

## Plotting
#==========================================================================================#

using Plots
pyplot(dpi = 300)

n_pts = 20
circle3D(r, n) = [ (r*cos(θ), 0, r*sin(θ)) for θ in 0:2π/n:2π ];
circ     = circle3D(R, n_pts) 

beam_pts = collect(zip(tupvector(fem_pts[1:end-1]), tupvector(fem_pts[2:end])))
left_pts = [ [ [ circ_pt .+ pt[1], circ_pt .+ pt[2] ] for circ_pt in circ ] for pt in beam_pts ]

mid_pts = midpoint.(wing_pans)

mid_pans = ac_mine # midpoint.(panels)[:]
mid_xs   = getindex.(mid_pans, 1)
mid_ys   = getindex.(mid_pans, 2)
mid_zs   = getindex.(mid_pans, 3)

CFs = force_coefficient.(reduce(vcat, F_S), dynamic_pressure(aero_state.rho_ref, aero_state.speed), aero_state.area_ref) 
# surface_force_coefficients(aero_surfs[1], aero_state)
Cxs = @. getindex(CFs, 1)
Cys = @. getindex(CFs, 2)
Czs = @. getindex(CFs, 3)

hs_pts = reduce(hcat, points.(collect(eachrow(bound_leg.(horsies)))))
hs_xs  = getindex.(hs_pts, 1)
hs_ys  = getindex.(hs_pts, 2)
hs_zs  = getindex.(hs_pts, 3)

# Plot
b = aero_state.span_ref
plot(camera = (45, 45), 
     xlim = (-b/2, b/2),
    #  ylim = (-b/2, b/2), 
     zlim = (-b/2, b/2)
    )

# Panels
[ plot!(pans, color = :black, label = ifelse(i == 1, "Panels", :none)) for (i, pans) in enumerate(plot_panels(panels[:])) ]

# Planform
plot!(wing_plan, color = :blue, label = "Planform")

# Beams
[ plot!(reduce(vcat, pt), color = :green, label = ifelse(i == 1, "Beams", :none)) for (i, pt) in enumerate(left_pts) ]

# Forces
quiver!(hs_xs[:], hs_ys[:], hs_zs[:], quiver=(Cxs[:], Cys[:], Czs[:]) .* 10, label = "Forces")
quiver!(mid_xs[:], mid_ys[:], mid_zs[:], quiver=(getindex.(forces[:], 1),getindex.(forces[:], 2),getindex.(forces[:], 3)) .* 0.1, label = "Forces")

# Axis system
# quiver!(getindex.(mid_pts, 1)[:], getindex.(mid_pts, 2)[:], getindex.(mid_pts, 3)[:], quiver=(getindex.(n_cs, 1)[:], getindex.(n_cs, 2)[:], getindex.(n_cs, 3)[:]), color = :orange, label = :none)
# quiver!(getindex.(mid_pts, 1)[:], getindex.(mid_pts, 2)[:], getindex.(mid_pts, 3)[:], quiver=(getindex.(ns, 1)[:], getindex.(ns, 2)[:], getindex.(ns, 3)[:]), color = :red, label = :none)
# quiver!(getindex.(mid_pts, 1)[:], getindex.(mid_pts, 2)[:], getindex.(mid_pts, 3)[:], quiver=(getindex.(ss, 1)[:], getindex.(ss, 2)[:], getindex.(ss, 3)[:]), color = :brown, label = :none)

# Nodes
scatter!(ac_pts[:,:,1][:], ac_pts[:,:,2][:], ac_pts[:,:,3][:], label = "Aerodynamic Centers")
# scatter!(getindex.(ac_mine, 1)[:], getindex.(ac_mine, 2)[:], getindex.(ac_mine, 3)[:], label = "ACs (Mine)")
plot!()