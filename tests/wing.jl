## 
using Revise
includet("../src/AeroMDAO.jl")
includet("../src/MathTools.jl")
includet("../src/FoilParametrization.jl")

##
using .FoilParametrization: read_foil, foil_camthick, camthick_foil, cosine_foil, kulfan_CST, naca4
using .AeroMDAO
using .MathTools: linspace, tuparray, tupvector
using DelimitedFiles
using Rotations

## Wing section setup
alpha_u = [0.1, 0.3, 0.2, 0.15, 0.2]
alpha_l = [-0.1, -0.1, -0.1, -0.001, -0.02]
alphas = [alpha_u alpha_l]
dzs = (1e-4, 1e-4)
cst_foil = kulfan_CST(alphas, dzs, 0.2)

num_secs = 3
foils = [ cst_foil for i in 1:num_secs ]
airfoils = Foil.(foils)

wing_chords = [1, 0.5, 0.2]
wing_twists = [0, 0, 0]
wing_spans = [2, 0.5]
wing_dihedrals = [0, 45]
wing_sweeps = [15, 30]

wing_right = HalfWing(airfoils, wing_chords, wing_spans, wing_dihedrals, wing_sweeps, wing_twists)
wing = Wing(wing_right, wing_right)
print_info(wing_right)

## Assembly
ρ = 1.225
uniform = Uniform(10.0, 1.0, 0.0)
V = uniform.mag
S = projected_area(wing_right)
c = mean_aerodynamic_chord(wing_right);

# ## Horseshoe method
# horseshoe_panels = mesh_horseshoes(wing_right, spanwise_panels = 15, chordwise_panels = 1)
# @time forces, moments, stable_force, stable_moment, drag = solve_horseshoes(horseshoe_panels, uniform, ρ)

# # Post-processing
# cl = force_coefficient(stable_force[3], ρ, V, S)
# cdi = force_coefficient(drag, ρ, V, S)
# cm = moment_coefficient(sum(stable_moment), ρ, V, S, c)

# println("Total Force: $stable_force N")
# println("Total Moment: $stable_moment N")
# println("Lift Coefficient: $cl, Drag Coefficient: $cdi, Moment Coefficient: $cm")
# println("Lift-to-Drag Ratio (L/D): $(cl/cdi)")


## Vortex lattice method
camber_panels = mesh_cambers(wing_right, spanwise_panels = 5, chordwise_panels = 20)
@time Γs, camber_panels = solve_vortex_rings(camber_panels, uniform, 10, 2.0, ρ)

## Post-processing
cl = force_coefficient(forces[3], ρ, V, S)
cdi = force_coefficient(drag, ρ, V, S)
cm = moment_coefficient(sum(moments), ρ, V, S, c)

println("Forces: $forces N")
println("Moments: $moments N")
println("Lift Coefficient: $cl, Drag Coefficient: $cdi, Moment Coefficient: $cm")


# ## Panel method: TO DO
# wing_panels = mesh_wing(wing, spanwise_panels = 5, chordwise_panels = 10);

##
using Plots, LaTeXStrings
##
plotlyjs()


## Wing
# horseshoe_collocs = horseshoe_collocation.(horseshoe_panels)[:]
camber_collocs = horseshoe_collocation.(camber_panels)[:]

##
# spans = [ pt[2] for pt in horseshoe_collocs ];
# plot(spans, cl, label = L"C_L")
# plot!(spans, cdi, label = L"C_{D_i}")

##
wing_coords = (tuparray ∘ panel_coords).(wing_panels)[:]
# horseshoe_coords = tuparray.(panel_coords.(horseshoe_panels)[:])
camber_coords = (tuparray ∘ panel_coords).(camber_panels)[:]
# vortex_rings = (tupvector ∘ vortex_ring).(camber_panels)

##
plot(xaxis = "x", yaxis = "y", zaxis = "z", aspect_ratio = :equal, zlim = (-0.5, 5.0), size=(1280, 720))
# plot!.(wing_coords, color = :black, label = :none)
# plot!.(horseshoe_coords, color = :grey, label = :none)
plot!.(camber_coords, color = :grey, label = :none)
# plot!.(vortex_rings, color = :black, label = :none)

# gui();

# scatter!(horseshoe_collocs, c = :grey, markersize = 1, label = "Wing Collocation Points")
scatter!(camber_collocs, c = :grey, markersize = 1, label = "Wing Collocation Points")

##
wing_lines = [ horseshoe_lines(panel, uniform) for panel in horseshoe_panels ]

lines = [ [ tuparray([line.r1'; line.r2']) for line in horseshoe.vortex_lines ] for horseshoe in wing_lines ]
[ [ plot!(line, c = :darkblue, label = :none) for line in horseshoe ] for horseshoe in lines ]

gui();