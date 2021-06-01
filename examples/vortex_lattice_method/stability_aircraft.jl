using AeroMDAO
using LinearAlgebra # For norm()

## Trapezoidal lifting surface - consists of one section
TrapezoidalWing(b, δ, Λ, λ, c_root, τ_root, τ_tip, foil_root, foil_tip) =
	HalfWing([ Foil(foil_root), Foil(foil_tip) ], # Foils
			   [c_root, λ * c_root], 			  # Chords
			   [τ_root, τ_tip], 				  # Twists
			   [b],             				  # Span
			   [δ],             				  # Dihedral
			   [Λ])             				  # LE sweep

# Lifting surfaces
wing_right  = TrapezoidalWing(4.0, 0.0, 15.0, 0.4, 2.0, 0.0, -2.0, naca4((2,4,1,2)), naca4((2,4,1,2)))
wing        = Wing(wing_right, wing_right)
wing_mac 	= mean_aerodynamic_center(wing)
wing_pos    = [0., 0., 0.]
wing_plan  	= plot_wing(wing;  
						position = wing_pos)

print_info(wing, "Wing")

htail_right = TrapezoidalWing(1.0, 0.0, 15.0, 0.6, 0.8, 0.0, 0.0, naca4((0,0,1,2)), naca4((0,0,0,9)));
htail		= Wing(htail_right, htail_right)
htail_mac	= mean_aerodynamic_center(htail)
htail_pos	= [5., 0., 0.]
α_h_i		= 0.
htail_plan	= plot_wing(htail;
						position = htail_pos)

print_info(htail, "Horizontal Tail")

vtail		= TrapezoidalWing(0.8, 0.0, 8.0, 0.6, 0.8, 0.0, 0., naca4((0,0,0,9)), naca4((0,0,0,9)))
vtail_mac	= mean_aerodynamic_center(vtail) # NEEDS FIXING FOR ROTATION
vtail_pos	= [5., 0., 0.]
vtail_plan	= plot_wing(vtail; 
						position = vtail_pos,
						angle 	= π/2)
						
print_info(vtail, "Vertical Tail")

## Static stability
S, b, c = projected_area(wing), span(wing), mean_aerodynamic_chord(wing)
x_w = wing_pos + [ wing_mac[1], 0, 0 ]

# Horizontal tail volume coefficient
S_h = projected_area(htail)
x_h = htail_pos + [ htail_mac[1], 0, 0 ]
l_h = norm(x_h - x_w)
V_h = S_h * l_h / (S * c)
println("Horizontal TVC      V_h: $V_h")

# Vertical tail volume coefficient
S_v = projected_area(vtail)
x_v = vtail_pos + [ vtail_mac[1], 0, 0 ]
l_v = norm(x_v - x_w)
V_v = S_v * l_v / (S * b)
println("Vertical TVC        V_v: $V_v")

## Panelling and assembly
wing_panels  = 	panel_wing(wing, [20], 10;
            	           position = wing_pos
                          )
htail_panels =	panel_wing(htail, [10], 5;
				           position	= htail_pos,
				           angle 	= deg2rad(α_h_i),
				           axis 	= [0., 1., 0.]
				          )
vtail_panels = 	panel_wing(vtail, [10], 5;
                           position = vtail_pos,
                           angle    = π/2
						  )

aircraft = Dict("Wing" 			  	=> wing_panels,
				"Horizontal Tail" 	=> htail_panels,
                "Vertical Tail"     => vtail_panels);

## Evaluate case
ac_name = "My Aircraft"
ρ 		= 1.225
ref     = x_w
V, α, β = 1.0, 0.0, 0.0
Ω 		= [0.0, 0.0, 0.0]
fs 	    = Freestream(V, α, β, Ω)

@time dv_data = 
solve_stability_case(aircraft, fs; 
					 rho_ref     = ρ, 			# Reference density
					 r_ref       = ref, 		# Reference point for moments
					 area_ref    = S, 			# Reference area
					 span_ref    = b, 			# Reference span
					 chord_ref   = c, 			# Reference chord
					 name        = ac_name,		# Aircraft name
					 print       = true,		# Prints the results for the entire aircraft
					 print_components = true,	# Prints the results for each component
					);

## Process data
labels = (collect ∘ keys)(dv_data)	# Get aircraft component names from analysis
comp = labels[1]					# Pick your component
nf, ff, dvs = dv_data[comp];		# Get the nearfield, farfield, and stablity derivative coefficients
print_case(dv_data, comp)			# Pretty-print the results

## Aerodynamic quantities
nf_plane, ff_plane, dvs_plane = dv_data[labels[1]]

# Center of pressure
Cm 		= nf_plane[5]	# Moment coefficient
CL 		= nf_plane[3]	# Lift coefficient
cp 		= -c * Cm / CL

# Neutral point
Cm_α 	= dvs_plane[5,1]
CL_α 	= dvs_plane[3,1]
np   	= -c * Cm_α / CL_α

# Locations
x_np 	= [ ref[1] .+ np; zeros(2) ]	# Neutral point
x_cp 	= [ ref[1] .+ cp; zeros(2) ]	# Center of pressure

println("Aerodynamic Center x_ac: $(x_w[1]) m")
println("Neutral Point      x_np: $(x_np[1]) m")
println("Center of Pressure x_cp: $(x_cp[1]) m")

## Plotting everything
using Plots
gr(dpi = 300)

#
z_limit = b
aircraft_plot = 
plot(xaxis = "x", yaxis = "y", zaxis = "z",
	 aspect_ratio = 1, 
	 camera = (30, 60),
     xlim = (0, z_limit),
    #  ylim = (-z_limit/2, z_limit/2),
	 zlim = (-z_limit/2, z_limit/2),
	 size = (1280, 720)
	)

plot!(wing_plan, label = "Wing")
plot!(htail_plan, label = "Horizontal Tail")
plot!(vtail_plan, label = "Vertical Tail")

scatter!(tuple(x_w...), label = "Wing MAC")
scatter!(tuple(x_h...), label = "Horizontal Tail MAC")
scatter!(tuple(x_v...), label = "Vertical Tail MAC")

scatter!(tuple(x_np...), label = "Neutral Point")
scatter!(tuple(x_cp...), label = "Center of Pressure")
# savefig(aircraft_plot, "Aircraft.png")
plot!()