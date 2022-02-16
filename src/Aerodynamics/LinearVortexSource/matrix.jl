## Panel helpers
#============================================#

# Evaluation of velocities on points
source_velocity(σ1, σ2, panel :: AbstractPanel2D, x, y) = panel_velocity(linear_source_velocity_a, linear_source_velocity_b, σ1, σ2, panel, x, y)
vortex_velocity(γ1, γ2, panel :: AbstractPanel2D, x, y) = panel_velocity(linear_vortex_velocity_a, linear_vortex_velocity_b, γ1, γ2, panel, x, y)

total_velocity(σ_1j :: Real, σ_2j :: Real, γ_1j :: Real, γ_2j :: Real, panel_j :: AbstractPanel2D, x, y) = source_velocity(σ_1j, σ_2j, panel_j, x, y) .+ vortex_velocity(γ_1j, γ_2j, panel_j, x, y)

# Evaluation of velocities on panels
source_velocity(σ_1j, σ_2j, panel_j :: AbstractPanel2D, panel_i :: AbstractPanel2D) = panel_velocity(linear_source_velocity_a, linear_source_velocity_b, σ_1j, σ_2j, panel_j, panel_i)
vortex_velocity(γ_1j, γ_2j, panel_j :: AbstractPanel2D, panel_i :: AbstractPanel2D) = panel_velocity(linear_vortex_velocity_a, linear_vortex_velocity_b, γ_1j, γ_2j, panel_j, panel_i)

total_velocity(σ_1j, σ_2j, γ_1j, γ_2j, panel_j :: AbstractPanel2D, panel_i :: AbstractPanel2D) = source_velocity(σ_1j, σ_2j, panel_j, panel_i) .+ vortex_velocity(γ_1j, γ_2j, panel_j, panel_i)

# Influence coefficient computation
influence_coefficient(velocity_func :: F, angle_func, panel_j :: AbstractPanel2D, panel_i :: AbstractPanel2D) where F <: Function = dot(panel_velocity(velocity_func, 1., panel_j, panel_i), angle_func(panel_i))

## Matrix construction
#============================================#

function two_point_matrix(vel_a :: F, vel_b :: G, panels_1, panels_2) where F where G
    N     = min(length(panels_1), length(panels_2))
    inf_a = neumann_influence_matrix(vel_a, normal_vector, panels_1, panels_2)
    inf_b = neumann_influence_matrix(vel_b, normal_vector, panels_1, panels_2)
    [ inf_a zeros(N) ] + [ zeros(N) inf_b ]
end

linear_source_matrix(panels_1, panels_2) = two_point_matrix(linear_source_velocity_a, linear_source_velocity_b, panels_1, panels_2)
linear_vortex_matrix(panels_1, panels_2) = two_point_matrix(linear_vortex_velocity_a, linear_vortex_velocity_b, panels_1, panels_2)

kutta_condition(panels_1, panels_2) = [ 1 zeros(length(panels_1) - 1)' 1 zeros(length(panels_2) - length(panels_1))' ]
source_influence_matrix(panels_1, panels_2) = linear_source_matrix(panels_1, panels_2)
vortex_influence_matrix(panels_1, panels_2) = linear_vortex_matrix(panels_1, panels_2)

# Diagonal cases
kutta_condition(panels) = [ 1 zeros(length(panels) - 1)' 1 ]
source_influence_matrix(panels) = [ linear_source_matrix(panels, panels); kutta_condition(panels) ]
vortex_influence_matrix(panels) = [ linear_vortex_matrix(panels, panels); kutta_condition(panels) ]

# Constant-strength sources
source_velocity(σ1, panel :: AbstractPanel2D, x, y) = panel_velocity(constant_source_velocity, σ1, panel, x, y)
source_velocity(σ_j, panel_j :: AbstractPanel2D, panel_i :: AbstractPanel2D) = panel_velocity(constant_source_velocity, σ_j, panel_j, panel_i)

constant_source_influence_coefficient(panel_j, panel_i) = ifelse(panel_i === panel_j, 0.5, influence_coefficient(constant_source_velocity, normal_vector, panel_j, panel_i))

constant_source_matrix(panels) = [ constant_source_influence_coefficient(panel_j, panel_i) for (panel_j, panel_i) in product(panels, panels)]

neumann_influence_matrix(func, angle_func, panel_is, panel_js) = [ influence_coefficient(func, angle_func, panel_j, panel_i) for (panel_i, panel_j) in product(panel_is, panel_js) ]

# Boundary conditions
#============================================#


constant_source_boundary_condition(panels, u) = -map(pan -> dot(u, normal_vector(pan)), panels)

neumann_boundary_condition(panels, u) = [ constant_source_boundary_condition(panels, u); 0 ]