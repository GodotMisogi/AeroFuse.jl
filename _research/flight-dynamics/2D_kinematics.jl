using StaticArrays
using LinearAlgebra
using BenchmarkTools
using DifferentialEquations
using ModelingToolkit
using DataFrames
using NLsolve
using AeroMDAO
using Plots

## Aerodynamics
include("aerodynamics.jl")

# Plotting
# plot(15:90, post_stall_CL.(deg2rad.(15:90), 6, 1.2, deg2rad(15)))
# plot!(15:90, post_stall_CD.(deg2rad.(15:90), 6, 0.5, deg2rad(15), 0.12))

## Propulsion
include("propulsion.jl")

# Constants
V       = 10
κ       = 1.
ρ       = 1.225
A_disk  = 1.0
P       = 1/2 * ρ * A_disk * V^3

# Problem setup
momentum_power_residual(P, T, V_perp, A_disk, ρ, κ) = P - momentum_power(T, V_perp, A_disk, ρ, κ)
solve_thrust!(R, T) = R .= momentum_power_residual(P, T[1], V, A_disk, ρ, κ)
T0   = 10.
# prob = nlsolve(solve_thrust!, [T0])
# prob.zero

## Flight dynamics
struct FlightState2D{T <: Real}
	position 		:: SVector{3,T}
	velocity 		:: SVector{3,T}
end

FlightState2D(position :: FieldVector{3,T}, velocity :: FieldVector{3,T}) where T <: Real = FlightState2D{T}(position, velocity)

# Get state vector
vector(state :: FlightState2D{T}) where T <: Real = [ state.position; state.velocity ]

# Equations of motion
linear_acceleration(F, mass, g) = F / mass + g
angular_acceleration(M, I, h)   = M / I + h

rotation_2D(α) = [  cos(α) sin(α) ;
                   -sin(α) cos(α) ]

inverse_rotation_2D(α) = rotation_2D(-α)

# Compute system of equations
function aircraft_eom_2D!(dx, x, p, t)
    # State variables
	r_e		= @view x[1:2] # Unused
    α       = @view x[3]
	u_b		= @view x[4:5]
    ω		= @view x[6]

    # State derivatives
	dr_dt 	= @view dx[1:2]
    dα_dt   = @view dx[3]
	du_dt 	= @view dx[4:5]
	dω_dt 	= @view dx[6]

    # Parametric inputs
	mass, g, F_b, M_b, I, h	= p(t, r_e, α[1], u_b, ω[1])

    # Equations of motion
	dr_dt  .= inverse_rotation_2D(α[1]) * u_b
	dα_dt  .= ω
	du_dt  .= linear_acceleration(F_b, mass, g)
	dω_dt  .= angular_acceleration(M_b, I, h)
end

# Dimensionalise forces
function thrust_force(V, ρ, α, δ_T)
    κ       = 1.
    A_disk  = 1.0
    P       = 1/2 * ρ * A_disk * V^3

    # Thrust evaluation
    solve_thrust!(R, T) = R .= momentum_power_residual(P, T[1], V, A_disk, ρ, κ)
    T0       = 3 / 2 * ρ * A_disk * V^2
    prob     = nlsolve(solve_thrust!, [T0])
    T_y, T_x = prob.zero[1] .* sincos(α)

    δ_T .* SVector(T_x, T_y)
end

function aerodynamic_dynamics(α, δ_e, δ_f, ρ, V, S, c)
    # Wing
    wing_AR           = 6.
    wing_Cl_α         = 5.6
    wing_Cl_s         = 2.2
    wing_α_s          = 15.
    wing_foil_max_tc  = 0.12
    wing_k            = 1.
    wing_Cd_max       = 2.
    wing_Cd_0         = 1e-3
    wing              = Airfoil(global_cl, global_cd, identity)
    cl_wing           = cl(wing, α, wing_Cl_α, wing_Cl_s, wing_α_s, wing_AR)
    cd_wing           = cd(wing, α, wing_Cd_0, wing_Cd_max, wing_α_s, wing_AR, wing_foil_max_tc, wing_k)

    # Horizontal tail
    htail_AR          = 3.
    htail_Cl_α        = 3.0
    htail_Cl_s        = 1.2
    htail_α_s         = 8.
    htail_foil_max_tc = 0.09
    htail_k           = 1.
    htail_Cd_max      = 1.2
    htail_Cd_0        = 1e-4
    htail             = Airfoil(global_cl, global_cd, identity)
    cl_htail          = cl(htail, α, htail_Cl_α, htail_Cl_s, htail_α_s, htail_AR)
    cd_htail          = cd(htail, α, htail_Cd_0, htail_Cd_max, htail_α_s, htail_AR, htail_foil_max_tc, htail_k)

    CF = (1 + δ_f) .* SVector(cl_wing, cd_wing) .+ (1 + δ_e) .* SVector(cl_htail, cd_htail)
    Cm = -0.1
    q  = dynamic_pressure(ρ, V)

    aerodynamic_dynamics(CF, Cm, q, S, c)
end

aerodynamic_dynamics(CFs, Cm, q, S, c) = CFs .* q * S, Cm .* q * S * c


function params(t, r_e, α, u_b, ω, h, δ_e, δ_f, δ_T, ρ, S, c, mass)
    V = norm(u_b)

	# Evaluate forces
	F_T      = thrust_force(V, ρ, α, δ_T)
    F_A, M_A = aerodynamic_dynamics(α, δ_e, δ_f, ρ, V, S, c)

	# Body force equilibrium
	F = F_T - F_A

    g_b = rotation_2D(α) * [0.0, -9.8]
    inertia = 1.0

	mass, g_b, F, M_A, inertia, h
end

## Setup
ρ    = 1.225
h    = 0.
δ_e  = 0.1
δ_f  = 0.1
δ_T  = 0.8
S    = 1.98
c    = 0.198
mass = 3.

cock(t, r_e, α, u_b, ω) = params(t, r_e, α, u_b, ω, h, δ_e, δ_f, δ_T, ρ, S, c, mass)

r0 = SVector(1.0, 0.0, deg2rad(5.0))
v0 = SVector(1.0, 1.0, 0.1)

state = FlightState2D(r0, v0)
x0 	  = (collect ∘ vector)(state)

## Run case
tspan = (0., 1.0)
run   = ODEFunction(aircraft_eom_2D!, syms = [:xₑ, :yₑ, :α, :u, :v, :ω])
prob  = ODEProblem(run, x0, tspan, cock)
sol   = @time solve(prob, RK4(), reltol = 1e-6, abstol = 1e-6);
df    = DataFrame(sol)

##
gr(dpi = 300)
plot(sol, layout = (2,3))
plot!()