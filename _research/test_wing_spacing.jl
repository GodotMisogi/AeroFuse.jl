using AeroMDAO

# Define wing
wing = Wing(foils     = Foil.(naca4((2,4,1,2)) for i ∈ 1:3),
            chords    = [0.18, 0.16, 0.08],
            twists    = [2., 0., -2.],
            spans     = [0.5, 0.5],
            dihedrals = [0., 11.3],
            sweep_LEs = [1.14, 8.])

# Define reference values
ρ   = 1.225
ref = [0.25 * mean_aerodynamic_chord(wing), 0., 0.]
Ω   = [0.0, 0.0, 0.0]

# Define freestream condition
uniform = Freestream(10.0, 2.0, 2.0, Ω)

# Evaluate stability case
nf_coeffs, ff_coeffs, CFs, CMs, horseshoe_panels, normals, horseshoes, Γs = solve_case(wing, uniform; rho_ref = ρ, r_ref = ref, span_num = 20, chord_num = 5, spacing = ["cosine", "uniform"])


print_coefficients(nf_coeffs, ff_coeffs)

## Display
horseshoe_coords = plot_panels(horseshoe_panels[:])
wing_coords      = plot_wing(wing);
colpoints        = horseshoe_point.(horseshoe_panels)

# Coordinates
xs = getindex.(colpoints, 1);
ys = getindex.(colpoints, 2);
zs = getindex.(colpoints, 3);

z_limit = 5
plotly(dpi = 300)
plot(xaxis = "x", yaxis = "y", zaxis = "z",
     aspect_ratio = 1, 
     camera = (30, 60),
     zlim = (-0.1, z_limit),
     size = (800, 600))
plot!.(horseshoe_coords, color = :black, label = :none)
scatter!(tupvector(colpoints)[:], marker = 1, color = :black, label = :none)
plot!()