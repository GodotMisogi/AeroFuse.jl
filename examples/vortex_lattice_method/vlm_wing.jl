## Wing analysis case
using AeroMDAO

## Wing section setup
wing_right = HalfWing(foils     = Foil.(fill(naca4((0,0,1,2)), 3)),
                      chords    = [1.0, 0.6, 0.2],
                      twists    = [0.0, 0.0, 0.0],
                      spans     = [5.0, 0.5],
                      dihedrals = [5., 5.],
                      LE_sweeps = [5., 5.]);
wing = Wing(wing_right, wing_right)
wing_mac = mean_aerodynamic_center(wing)
S, b, c = projected_area(wing), span(wing), mean_aerodynamic_chord(wing);
print_info(wing, "Wing")

## Assembly
ρ       = 1.225
x_w     = wing_mac[1]
ref     = [ x_w; 0.; 0. ]
V, α, β = 1.0, 5.0, 0.0
Ω       = [0.0, 0.0, 0.0]
fs      = Freestream(V, α, β, Ω)

## Evaluate case
@time nf_coeffs, ff_coeffs, CFs, CMs, horseshoe_panels, normies, horses, Γs =
solve_case(wing, fs;
           rho_ref   = ρ,
           r_ref     = ref,
           area_ref  = S,
           span_ref  = b,
           chord_ref = c,
           span_num  = [25, 8],
           chord_num = 5,
           viscous   = true, # Only appropriate for α = β = 0, but works for other angles anyway
           x_tr      = [0.3, 0.3],
          #  spacing   = Uniform()
          );

print_coefficients(nf_coeffs, ff_coeffs, "Wing")

## Evaluate case with stability derivatives
@time nf, ff, dvs =
solve_stability_case(wing, fs;
                     rho_ref    = ρ,
                     r_ref      = ref,
                     span_num   = [25, 8],
                     chord_num  = 5,
                     name       = "My Wing",
                     viscous    = true,
                     x_tr       = [0.3, 0.3],
                     print      = false,
                    #  spacing   = [Cosine(), "cosine"]
                    );

#
print_coefficients(nf, ff, "Wing")
print_derivatives(dvs, "Wing")

## Plotting
using StaticArrays
using Plots
gr()

## Streamlines

# Chordwise distribution
# num_points = 50
# max_z = 0.1
# y = span(wing) / 2 - 0.05
# seed = SVector.(fill(-0.1, num_points), fill(y, num_points), range(-max_z, stop = max_z, length = num_points))

# Spanwise distribution
span_points = 20
init        = chop_trailing_edge(wing, span_points)
dx, dy, dz  = 0, 0, 1e-3
seed        = [ init .+ Ref([dx, dy,  dz])  ;
                init .+ Ref([dx, dy, -dz]) ];

distance = 5
num_stream_points = 100
streams = plot_streams(fs, seed, horses, Γs, distance, num_stream_points);

## Display
horseshoe_coords = plot_panels(horseshoe_panels)
wing_coords      = plot_wing(wing);
horseshoe_points = Tuple.(horseshoe_point.(horseshoe_panels))
ys               = getindex.(horseshoe_points[1,:], 2)

## Coordinates
z_limit = 5
plot(xaxis = "x", yaxis = "y", zaxis = "z",
     aspect_ratio = 1,
     camera = (30, 60),
     zlim = (-0.1, z_limit),
     size = (800, 600))
plot!.(horseshoe_coords, color = :black, label = :none)
scatter!(horseshoe_points[:], marker = 1, color = :black, label = :none)
plot!.(streams, color = :green, label = :none)
plot!()

## Spanwise forces
wind_CFs = body_to_wind_axes.(CFs, fs.alpha, fs.beta)
CDis     = @. getindex(wind_CFs, 1)
CYs      = @. getindex(wind_CFs, 2)
CLs      = @. getindex(wind_CFs, 3)

area_scale  = S ./ sum(panel_area, horseshoe_panels, dims = 1)[:]
span_CDis   = sum(CDis, dims = 1)[:] .* area_scale
span_CYs    = sum(CYs,  dims = 1)[:] .* area_scale
span_CLs    = sum(CLs,  dims = 1)[:] .* area_scale
CL_loadings = sum(Γs,   dims = 1)[:] / (0.5 * fs.V * c)

plot_CD = plot(ys, span_CDis, label = :none, ylabel = "CDi")
plot_CY = plot(ys, span_CYs, label = :none, ylabel = "CY")
plot_CL = begin
            plot(ys, span_CLs, label = :none, xlabel = "y", ylabel = "CL")
            plot!(ys, CL_loadings, label = "Normalized", xlabel = "y")
          end
plot(plot_CD, plot_CY, plot_CL, size = (800, 700), layout = (3,1))

## Lift distribution

# Exaggerated CF distribution for plot
hs_pts = Tuple.(bound_leg_center.(horses))[:]

plot(xaxis = "x", yaxis = "y", zaxis = "z",
     aspect_ratio = 1,
     camera = (60, 60),
     zlim = (-0.1, z_limit)
    )
plot!.(horseshoe_coords, color = :gray, label = :none)
# scatter!(cz_pts, zcolor = CLs[:], marker = 2, label = "CL (Exaggerated)")
quiver!(hs_pts, quiver=(CDis[:], CYs[:], CLs[:]) .* 100, label = "Forces (Exaggerated)")
plot!(size = (800, 600))