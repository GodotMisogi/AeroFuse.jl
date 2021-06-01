## Panel helpers
#============================================#

# Evaluation of velocities on points
source_velocity(σ1, σ2, panel :: Panel2D, x, y) = panel_velocity(linear_source_velocity_a, linear_source_velocity_b, σ1, σ2, panel, x, y)
vortex_velocity(γ1, γ2, panel :: Panel2D, x, y) = panel_velocity(linear_vortex_velocity_a, linear_vortex_velocity_b, γ1, γ2, panel, x, y)

total_velocity(σ_1j :: Real, σ_2j :: Real, γ_1j :: Real, γ_2j :: Real, panel_j :: Panel2D, x, y) = source_velocity(σ_1j, σ_2j, panel_j, x, y) .+ vortex_velocity(γ_1j, γ_2j, panel_j, x, y)

# Evaluation of velocities on panels
source_velocity(σ_1j, σ_2j, panel_j :: Panel2D, panel_i :: Panel2D) = panel_velocity(linear_source_velocity_a, linear_source_velocity_b, σ_1j, σ_2j, panel_j, panel_i)
vortex_velocity(γ_1j, γ_2j, panel_j :: Panel2D, panel_i :: Panel2D) = panel_velocity(linear_vortex_velocity_a, linear_vortex_velocity_b, γ_1j, γ_2j, panel_j, panel_i)

total_velocity(σ_1j, σ_2j, γ_1j, γ_2j, panel_j :: Panel2D, panel_i :: Panel2D) = source_velocity(σ_1j, σ_2j, panel_j, panel_i) .+ vortex_velocity(γ_1j, γ_2j, panel_j, panel_i)

## Matrix construction
#============================================#

linear_vortex_matrix(panels_1 :: Vector{<: Panel}, panels_2 :: Vector{<: Panel}) = two_point_matrix(linear_vortex_velocity_a, linear_vortex_velocity_b, panels_1, panels_2)
linear_source_matrix(panels_1 :: Vector{<: Panel}, panels_2 :: Vector{<: Panel}) = two_point_matrix(linear_source_velocity_a, linear_source_velocity_b, panels_1, panels_2)

kutta(panels_1 :: Vector{<: Panel}, panels_2 :: Vector{<: Panel}) = [ 1 zeros(length(panels_1) - 1)' 1 zeros(length(panels_2) - length(panels_1))' ]
source_influence_matrix(panels_1 :: Vector{<: Panel}, panels_2 :: Vector{<: Panel}) = [ linear_source_matrix(panels_1, panels_2); kutta(panels_1, panels_2) ]
vortex_influence_matrix(panels_1 :: Vector{<: Panel}, panels_2 :: Vector{<: Panel}) = [ linear_vortex_matrix(panels_1, panels_2); kutta(panels_1, panels_2) ]

# Diagonal cases
kutta(panels :: Vector{<: Panel}) = [ 1 zeros(length(panels) - 1)' 1 ]
source_influence_matrix(panels :: Vector{<: Panel}) = [ linear_source_matrix(panels, panels); kutta(panels) ]
vortex_influence_matrix(panels :: Vector{<: Panel}) = [ linear_vortex_matrix(panels, panels); kutta(panels) ]

# Constant-strength sources
source_velocity(σ1, panel :: Panel2D, x, y) = panel_velocity(constant_source_velocity, σ1, panel, x, y)
source_velocity(σ_j, panel_j :: Panel2D, panel_i :: Panel2D) = panel_velocity(constant_source_velocity, σ_j, panel_j, panel_i)

constant_source_influence_coefficient(panel_j, panel_i) = ifelse(panel_i === panel_j, 0.5, influence_coefficient(constant_source_velocity, panel_normal, panel_j, panel_i))

constant_source_matrix(panels :: Vector{<: Panel2D}) = [ constant_source_influence_coefficient(panel_j, panel_i) for (panel_j, panel_i) in product(panels, panels)]

constant_source_boundary_condition(panels :: Vector{<: Panel2D}, u) = -dot.(Ref(u), panel_normal.(panels))


# Generic influence coefficient computation on panels
#============================================#

influence_coefficient(velocity_func :: F, angle_func, panel_j :: Panel, panel_i :: Panel) where F <: Function = dot(panel_velocity(velocity_func, 1., panel_j, panel_i), angle_func(panel_i))

neumann_influence_matrix(func, angle_func, panel_is :: Vector{<: Panel}, panel_js :: Vector{<: Panel}) = [ influence_coefficient(func, angle_func, panel_j, panel_i) for (panel_i, panel_j) in product(panel_is, panel_js) ]

function two_point_matrix(vel_a :: F, vel_b :: G, panels_1 :: Vector{<: Panel}, panels_2 :: Vector{<: Panel}) where F where G
	N 	  = min(length(panels_1), length(panels_2))
	inf_a = neumann_influence_matrix(vel_a, panel_normal, panels_1, panels_2)
	inf_b = neumann_influence_matrix(vel_b, panel_normal, panels_1, panels_2)
	[ inf_a zeros(N) ] + [ zeros(N) inf_b ]
end

neumann_boundary_condition(panels :: Vector{<: Panel}, u) = [ -dot.(Ref(u), panel_normal.(panels)); 0 ]