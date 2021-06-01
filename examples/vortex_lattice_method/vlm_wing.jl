## 
using AeroMDAO

## Wing section setup
wing_foils = Foil.(fill(naca4((0,0,1,2)), 3))
wing_right = HalfWing(wing_foils,
                      [1.0, 0.6, 0.2],
                      [0.0, 0.0, 0.0],
                      [5.0, 0.5],
                      [5., 5.],
                      [5., 5.]);
wing = Wing(wing_right, wing_right)
print_info(wing, "Wing")
S, b, c = projected_area(wing), span(wing), mean_aerodynamic_chord(wing);

## Assembly
ρ 		= 1.225
ref 	= [0.25, 0., 0.]
V, α, β = 1.0, 1.0, 0.0
Ω 		= [0.0, 0.0, 0.0]
fs 		= Freestream(V, α, β, Ω)

## Evaluate case
@time nf_coeffs, ff_coeffs, CFs, CMs, horseshoe_panels, camber_panels, horseshoes, Γs = 
solve_case(wing, fs; 
           rho_ref   = ρ, 
           r_ref     = ref,
           area_ref  = S,
           span_ref  = b,
           chord_ref = c,
           span_num  = [15, 9], 
           chord_num = 6,
           viscous   = true, # Only appropriate for α = β = 0, but works for other angles anyway
           x_tr      = [0.3, 0.3]);

print_coefficients("Wing", nf_coeffs, ff_coeffs)

## Evaluate case with stability derivatives
@time nf, ff, dvs = 
solve_stability_case(wing, fs; 
                     rho_ref    = ρ, 
                     r_ref      = ref, 
                     area_ref   = S, 
                     span_ref   = b, 
                     chord_ref  = c, 
                     span_num   = [15, 9], 
                     chord_num  = 6, 
                     name       = "My Wing",
                     viscous    = true,
                     x_tr       = [0.3, 0.3],
                     print      = false);

#
print_coefficients("Wing", nf, ff)
print_derivatives("Wing", dvs)

## Plotting
using Plots
gr(size = (600, 400), dpi = 300)

## Coordinates
horseshoe_coords 	= plot_panels(horseshoe_panels[:])
camber_coords		= plot_panels(camber_panels[:])
wing_coords 		= plot_wing(wing);

CDis = getindex.(CFs, 1)
CYs	 = getindex.(CFs, 2)
CLs  = getindex.(CFs, 3);
CL_loadings = 2sum(Γs, dims = 1)[:] / (V * b)

colpoints 	= horseshoe_point.(horseshoe_panels)
xs 			= getindex.(colpoints, 1);
ys 			= getindex.(colpoints, 2);
zs 			= getindex.(colpoints, 3);
cl_pts 		= tupvector(SVector.(xs[:], ys[:], zs[:] .+ CLs[:]));

## Streamlines

# Chordwise distirbution
# num_points = 50
# max_z = 0.1
# y = span(wing) / 2 - 0.05
# seed = SVector.(fill(-0.1, num_points), fill(y, num_points), range(-max_z, stop = max_z, length = num_points))

# Spanwise distribution
span_points = 20
init        = trailing_chopper(ifelse(β == 0 && Ω == zeros(3), wing.right, wing), span_points) 
dx, dy, dz  = 0, 0, 1e-3
seed        = [ init .+ Ref([dx, dy, dz])  ; 
                init .+ Ref([dx, dy, -dz]) ];

distance = 5
num_stream_points = 100
streams = plot_streams(fs, seed, horseshoes, Γs, distance, num_stream_points);

## Display
z_limit = 5
plot(xaxis = "x", yaxis = "y", zaxis = "z",
     aspect_ratio = 1, 
     camera = (90,0),
     zlim = (-0.1, z_limit),
     size = (800, 600))
# plot!.(horseshoe_coords, color = :black, label = :none)
plot!.(camber_coords, color = :black, label = :none)
scatter!(tupvector(colpoints)[:], marker = 1, color = :black, label = :none)
plot!.(streams, color = :green, label = :none)
plot!()

## Span forces
plot1 = plot(ys[1,:], sum(CDis, dims = 1)[:], label = :none, ylabel = "CDi")
plot2 = plot(ys[1,:], abs.(sum(CYs, dims = 1)[:]), label = :none, ylabel = "CY")
plot3 = plot(ys[1,:], sum(CLs, dims = 1)[:], label = :none, xlabel = "y", ylabel = "CL")
plot(plot1, plot2, plot3, layout = (3,1))

## Lift distribution
plot(xaxis = "x", yaxis = "y", zaxis = "z",
    aspect_ratio = 1,
    camera = (30, 30),
    zlim = (-0.1, z_limit))
plot!(wing_coords, label = :none)
scatter!(cl_pts, zcolor = CLs[:], marker = 1, label = "CL")
plot!(size = (800, 600))

## Aircraft Assembly
wing_panels = WingSetup(wing, ([20], 10), orientation = AngleAxis{Float64}(0, 1, 0, 0), position = zeros(3));
aircraft = Dict("Wing" => wing_panels)

ac_name 	= "My Wing"

# reset_timer!()
@time data 	= solve_case(aircraft, fs; 
                         rho_ref     = ρ, 
                         r_ref       = ref, 
                         area_ref    = S, 
                         span_ref    = b, 
                         chord_ref   = c, 
                         name        = ac_name, 
                         # print       = true, 
                         #  print_components = true
                        ) 
# print_timer()
#
names = (collect ∘ keys)(data)
component = names[1]
nf_coeffs, ff_coeffs, CFs, CMs, horseshoe_panels, camber_panels, horseshoes, Γs = data[component];

print_coefficients(component, nf_coeffs, ff_coeffs)

## Panelling errors
panel_case(x, y) = solve_case(wing, fs; 
                              rho_ref   = ρ, 
                              r_ref     = ref,
                              area_ref  = S,
                              span_ref  = b,
                              chord_ref = c,
                              span_num  = y, 
                              chord_num = x,
                              viscous   = true,
                              x_tr      = [0.3, 0.3]);

xs = 6:30
ys = xs .* 4

nf_data, ff_data = zeros(5)', zeros(5)'
for (x, y) in zip(xs, ys)
    @time nf_coeffs, ff_coeffs, CFs, CMs, horseshoe_panels, camber_panels, horseshoes, Γs = 
    panel_case(x, y)

    println("Distribution: nx = $x, ny = $y")	
    print_coefficients("Wing", nf_coeffs, ff_coeffs)
    nf_data = [ nf_data; [x; y; nf_coeffs[3:5] ]' ]
    ff_data = [ ff_data; [x; y; ff_coeffs[3:5] ]' ]
end

##
nf_data

##
errors = @. abs((ff_data[3:end,3:end] - ff_data[2:end-1,3:end])/ff_data[2:end-1,3:end])

##
plot(ff_data[2:end-1,2], errors, xlabel = "Spanwise Panels", yscale = :log10)

##
plot(nf_data[2:end,2], nf_data[2:end,3], marker = :dot, label = "Nearfield")
plot!(ff_data[2:end,2], ff_data[2:end,3], marker = :dot, label = "Farfield")