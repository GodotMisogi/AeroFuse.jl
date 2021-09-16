using AeroMDAO
using NLsolve
using DifferentialEquations
using TaylorSeries
using Plots

## Aerodynamic coefficients (convert to struct for generality)
#===========================================================================#

# Lift coefficient
CL_0  = 0.895 # Lift coefficient at α = 0
CL_a  = 5.01  # Lift curve slope (Lift-AoA derivative)
CL_de = 0.722 # Lift-elevator derivative

lift_coefficient(α, δe, CL_0, CL_a, CL_δe) = CL_0 + CL_a * α + CL_δe * δe
lift_coefficient(α, δe) = lift_coefficient(α, δe, CL_0, CL_a, CL_de)

# Drag coefficient
CD₀      = 0.177 # Profile drag coefficient (at α when CL = 0? Hmm, this is inconsistent...)
∂CD_∂α   = 0.232 # Drag-AoA derivative
∂²CD_∂α² = 1.393 # Drag-AoA second derivative

drag_coefficient(α, CD₀, ∂CD_∂α, ∂²CD_∂α²) = CD₀ + ∂CD_∂α * α + ∂²CD_∂α² * α^2
drag_coefficient(α) = drag_coefficient(α, CD₀, ∂CD_∂α, ∂²CD_∂α²)

# Moment coefficient
CM_0  = -0.046 # Moment coefficient at α = 0?
CM_α  = -1.087 # Moment curve slope (Moment-AoA derivative)
CM_δe = -1.88  # Moment-elevator derivative
CM_Q̂  = -7.055 # Moment-dynamic pressure derivative

moment_coefficient(α, δe, q̂, CM_0, CM_α, CM_δe, CM_q̂) = CM_0 + CM_α * α + CM_δe * δe + CM_q̂ * q̂
moment_coefficient(α, δe, Q̂) = moment_coefficient(α, δe, Q̂, CM_0, CM_α, CM_δe, CM_Q̂)

function compute_aerodynamics(α, δe, q̂, q∞, S_ref, c_ref)
    # Compute aerodynamic coefficients from Taylor series approximations
    CD =   drag_coefficient(α)
    CL =   lift_coefficient(α, δe)
    CM = moment_coefficient(α, δe, q̂)

    # Compute dimensional aerodynamic forces and moment
    D =  force(CD, q∞, S_ref)
    L =  force(CL, q∞, S_ref)
    M = moment(CM, q∞, S_ref, c_ref)

    return D, L, M
end

## Equations of motion
#===========================================================================#

# Net forces
translational_forces(T, D, L, W, α, Θ) = [T; 0] + rotation(Θ) * [0; W] + rotation(α) * [D; L]
horizontal_forces(T, D, L, W, α, Θ)    = T - D * cos(α) + L * sin(α) - W * sin(Θ)
vertical_forces(D, L, W, α, Θ)         = -D * sin(α) - L * cos(α) + W * cos(Θ)
longitudinal_moment(M_A, T, Δ_zT)      = M_A - T * Δ_zT

## Trim analysis
#===========================================================================#

function total_dynamics!(R, T, D, L, W, M, α, Θ, Δ_zT)
    R[1] = horizontal_forces(T, D, L, W, α, Θ)
    R[2] = vertical_forces(D, L, W, α, Θ)
    R[3] = longitudinal_moment(M, T, Δ_zT)
end

function trim_equations!(R, x, params)
    # Unpack variables
    T  = x[1] # Thrust
    δe = x[2] # Elevator deflection angle
    Θ  = x[3] # Pitch angle

    # Other parameters
    W, Δ_zT, q∞, S_ref, c_ref, γ = params

    # Compute aerodynamics
    α       = Θ - γ # Angle of attack
    Q̂       = 0.    # Trim definition
    D, L, M = compute_aerodynamics(α, δe, Q̂, q∞, S_ref, c_ref)

    total_dynamics!(R, T, D, L, W, M, α, Θ, Δ_zT)
end

## Integration of 3-DOF EOM
#===========================================================================#

function compute_dynamics(Q, α, δe, T_in, mass, g, ρ, V_ref, S_ref, c_ref)
    # Dynamic pressure for dimensionalisation
    q∞ = dynamic_pressure(ρ, V_ref)

    # Aerodynamic forces and moment
    Q̂       = rate_coefficient(Q, V_ref, c_ref)
    D, L, M = compute_aerodynamics(α, δe, Q̂, q∞, S_ref, c_ref)

    # Propulsive forces
    # (Replace with some velocity-dependent curve with thrust coefficient?)
    T = T_in

    # Weight (Replace with fuel burn computation?)
    W = mass * g

    T, D, L, W, M
end

function longitudinal_equations_of_motion!(dx, x, params, t)
    # State variables at current timestep t
    u = x[1] # x-velocity (Forward speed)
    w = x[2] # z-velocity (Heave speed)
    Q = x[3] # Pitch rate
    Θ = x[6] # Pitch attitude

    δe = x[7] # Elevator deflection angle

    # Parameters
    mass, g, Iyy, Δ_zT, T_in, ρ, S_ref, c_ref = params

    # Compute freestream values
    speed, α = cartesian_to_freestream(u, w)

    # Compute forces and moment
    T, D, L, W, M = compute_dynamics(Q, α, δe, T_in, mass, g, ρ, speed, S_ref, c_ref)

    # Nonlinear longitudinal equations of motion
    #============================================#

    # Translational dynamics
    dx[1]   = -Q * w + horizontal_forces(T, D, L, W, α, Θ) / mass
    dx[2]   =  Q * u + vertical_forces(D, L, W, α, Θ) / mass

    # Longitudinal dynamics
    dx[3]   = longitudinal_moment(M, T, Δ_zT) / Iyy

    # Coordinate transformations to Earth axes
    dx[4:5] = inverse_rotation(Θ) * x[1:2]

    # Pitch rate
    dx[6]   = Q

    # Elevator deflection rate (Could be modified for a controller law later?)
    dx[7]   = 0.

    return dx
end

## Test case
#===========================================================================#

mass  = 7484.4 # Mass, kg
g     = 9.81   # Gravity, m/s²
Iyy   = 84309  # Moment of inertia in xz-plane, kg·m²
S_ref = 32.8   # Wing area, m²
c_ref = 2.29   # Reference chord length, m
Δ_zT  = -0.378 # Thrust vector location from CG

weight = mass * g

## Trim analysis
ρ    = 1.225
V, α = 120 * 0.514, 0.
q∞   = dynamic_pressure(ρ, V)
γ    = 0.

# Initial guess
T_0  = 15000.
Θ_0  = deg2rad(2)
δe_0 = deg2rad(0)

x0     = [ T_0, Θ_0, δe_0 ]
params = [ weight, Δ_zT, q∞, S_ref, c_ref, γ ]

# Closure
trim_equations!(R, x) = trim_equations!(R, x, params)

## Solve system
res_trim = nlsolve(trim_equations!,    # Nonlinear equation
                   x0,                 # Initial guess
                   method = :newton,   # Solution method
                  )

T_s, δe_s, Θ_s = trim_state = res_trim.zero

## ODE integration

# Initial state setup
Δδe            = deg2rad(-1.)
δe_init        = δe_s + Δδe

Θ_init         = Θ_s
α_init         = Θ_init - γ

u_init, w_init = freestream_to_cartesian(V, α_init)

x_e_init       = 0.
z_e_init       = 0.
Q_init         = 0
T_init         = T_s

# ODE initial state
x_init = [ u_init, w_init, Q_init, x_e_init, z_e_init, Θ_init, δe_init ]
ps     = [ mass, g, Iyy, Δ_zT, T_init, ρ, S_ref, c_ref ]
tspan  = (0, 100.)

## ODE setup
ode = ODEFunction(longitudinal_equations_of_motion!,
                  syms = [:u, :w, :Q, :xₑ, :zₑ, :Θ, :δe])

prob = ODEProblem(ode, x_init, tspan, ps)

## ODE solution
sol = solve(prob)

## Sensitivity analysis
# dg(out,u,p,t,i) = (out.=1.0.-u)
# ts = 0:0.5:100
# du0, dp = adjoint_sensitivities(sol, RK4(), dg, ts; sensealg=ZygoteAdjoint())

## Plotting
pyplot()

##
plot(sol, vars = [1,2,3,4,5,6], layout = (2,3))