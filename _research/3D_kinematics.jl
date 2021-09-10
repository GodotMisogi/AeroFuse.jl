using Rotations
using StaticArrays
using LinearAlgebra
using BenchmarkTools
using DifferentialEquations
using ModelingToolkit
using DataFrames
using AeroMDAO

struct StateVector{T <: Real}
	position 		:: SVector{3,T}
	euler_angles 	:: SVector{3,T}
	velocity 		:: SVector{3,T}
	rates 			:: SVector{3,T}
end

StateVector(position :: FieldVector{3,T}, euler_angles :: FieldVector{3,T}, velocity :: FieldVector{3,T}, rates :: FieldVector{3,T}) where T <: Real = StateVector{T}(position, euler_angles, velocity, rates)

vector(state :: StateVector{T}) where T <: Real = [ state.position; state.euler_angles; state.velocity; state.rates ]

# Euler angles, body to Earth
body_to_earth(euler_angles) = RotZYX(euler_angles...)
earth_to_body(euler_angles) = body_to_earth(euler_angles)'

# Orientation rates
function euler_rates_to_earth_rates(ψ, θ)
	sψ, cψ = sincos(ψ)
	sθ, cθ = sincos(θ)

	[ cψ * cθ -sψ  0 ;
	  sψ * cθ  cψ  0 ;
	    -sθ	   0   1 ]
end

function euler_rates_to_body_rates(θ, φ)
	sφ, cφ = sincos(φ)
	sθ, cθ = sincos(θ)

	[ 1   0    -sθ   ;
	  0  cφ  sφ * cθ ;
	  0 -sφ  cφ * cθ ]
end

function body_to_earth_rates(Φs)
	ψ, θ, φ = Φs
	sφ, cφ = sincos(φ)
	cθ, tθ = cos(θ), tan(θ)

	[ 1 sφ * tθ cφ * tθ ;
	  0    cφ  	  -sφ	;
	  0	sφ / cθ	cφ / cθ ]
end

function quaternion_matrix(ω)
	ω1, ω2, ω3 = ω ./ 2
	@SMatrix [ 0  -ω1 -ω2 -ω3 ;
			   ω1  0   ω3 -ω2 ; 
			   ω2 -ω3  0   ω1 ;  
			   ω3  ω2 -ω1  0  ]
end

# Kinematic equations
linear_momentum_body_axes(F, g, mass, U, a, Ω) = F + mass * (g - a + Ω × U)
angular_momentum_body_axes(M, I, Ω, α, h) = M - I * α + Ω × (I * Ω + h)

linear_acceleration(F, g, mass, U, Ω) = F / mass - (Ω × U - g)
angular_acceleration(M, I, Ω, h) = M - Ω × (I * Ω + h)

# Equations of motion
function aircraft_eom!(dx, x, p, t)
	r_e		= @view x[1:3] # Unused
	φ_b		= @view x[4:6]
	U_b		= @view x[7:9]
	p_b		= @view x[10:end]

	mass, g, F_b, M_b, I, h	= p(t, U_b, φ_b, p_b)

	dr_dt 	= @view dx[1:3]
	dφ_dt 	= @view dx[4:6]
	du_dt 	= @view dx[7:9]
	dp_dt 	= @view dx[10:end]

	x0 		= zeros(3)

	dr_dt 	.= mul!(x0, body_to_earth(φ_b), U_b)
	dφ_dt 	.= mul!(x0, body_to_earth_rates(φ_b), U_b)
	du_dt 	.= linear_acceleration(F_b, g, mass, U_b, p_b)
	dp_dt 	.= angular_acceleration(M_b, I, p_b, h)
end

# Trapezoidal lifting surface
TrapezoidalWing(b, δ, Λ, λ, c_root, τ_root, τ_tip, foil_root, foil_tip) =
HalfWing([ Foil(foil_root) , 
		   Foil(foil_tip) ], 		# Foils
		   [c_root, λ * c_root], 	# Chords
		   [τ_root, τ_tip], 		# Twists
		   [b],             		# Span
		   [δ],             		# Dihedral
		   [Λ])             		# LE sweep

# Surfaces
wing_right  = TrapezoidalWing(4.0, 0.0, 15.0, 0.4, 2.0, 0.0, -2.0, naca4((2,4,1,2)), naca4((2,4,1,2)))
wing        = Wing(wing_right, wing_right)
wing_mac 	= mean_aerodynamic_center(wing)
wing_pos    = [0., 0., 0.]
wing_plan  	= plot_wing(wing;  
						position = wing_pos)

print_info(wing, "Wing")

S, b, c = projected_area(wing), span(wing), mean_aerodynamic_chord(wing)
x_w = wing_pos + [ wing_mac[1], 0, 0 ]

ρ 		= 1.225
ref     = x_w
V, α, β = 5.0, 1.0, 0.0
Ω 		= @SVector [0., 0., 0.] 	# Angular velocity
fs 	    = Freestream(V, α, β, Ω)

mass = 250 		# Mass 
inertia = @SMatrix [ 1. 0. 0. ;
					 0. 1. 0. ;
					 0. 0. 1. ]  # Inertia matrix

O = zeros(3,3)
M = @SMatrix [ I O O	O	 ;
			   O I O	O	 ;
			   O O I	O	 ;
			   O O O inertia ]

φs 	= @SVector [0., 0., 0.] 	# Euler angles
α 	= @SVector [-1., 0., 1.] 	# Angular acceleration
R_o = @SVector [1., 1., 1.] 	# Position vector of aircraft
g 	= @SVector [0., 0., 1.] 	# Gravitational acceleration, Earth frame
U 	= velocity(fs)				# Freestream velocity

## Differential equations setup
thrust_force(u) = SVector(20u^2 + 10u + 20, 0., 0.)
aero_stability(CFs, CMs, q, S, b, c) = SVector(1., -1, 1) .* CFs .* q * S, CMs .* SVector(b, c, b) .* q * S

function params(t, U, φ, Ω, mass, g, inertia, h, ac, ρ, ref, S, b, c)
	fs = Freestream(U, zeros(3))
	q = dynamic_pressure(ρ, fs.V)

	# Evaluate forces
	F_T = thrust_force(fs.V)
	CFs, CMs = vlm_analysis(ac, fs, ρ, ref, S, b, c)
	F_A, M_A = aero_stability(CFs, CMs, q, S, b, c)

	# Force equilibrium
	F = F_T - F_A

	# M_A = zeros(3)
	
	g_b = earth_to_body(φ) * g

	# free = cartesian_to_freestream(velocity(fs))
	# rates = rate_coefficient(Ω, V, b, c)

	# println("Freestream: $free")
	# println("Gravity: $g_b")
	# println("Rotation: $rates")
	# println("Force coefficients: $CFs")
	# println("Moment coefficients: $CMs")

	mass, g_b, F, M_A, inertia, h
end

cuck(t, U, φs, Ω) = params(t, U, φs, Ω, mass, g, inertia, zeros(3), aircraft, ρ, ref, S, b, c)

##
@time cocker = cuck(1.0, U, φs, Ω)

## DifferentialEquations setup
state 	= StateVector(R_o, φs, U, Ω)
x0 		= (collect ∘ vector)(state)

tspan 	= (0., 5.0)
run 	= ODEFunction(aircraft_eom!, syms = [:xₑ, :yₑ, :zₑ, :φ, :θ, :ψ, :u, :v, :w, :p, :q, :r])
prob	= ODEProblem(run, x0, tspan, cuck, mass_matrix = M)
# de		= modelingtoolkitize(prob)
# jac		= eval(ModelingToolkit.generate_jacobian(de)[2])
# f		= ODEFunction(run, jac = jac, mass_matrix = M)

# prob_jac = ODEProblem(f, x0, tspan, cuck)

##
@time sol = solve(prob, DifferentialEquations.CVODE_BDF(linear_solver=:GMRES), reltol = 1e-6, abstol = 1e-6);
df = DataFrame(sol)

##
using Plots
gr(dpi = 300)
plot(sol, tspan = tspan, layout = (4,3), size = (900, 600))
plot!()