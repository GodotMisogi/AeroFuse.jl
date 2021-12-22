using Rotations
using StaticArrays, SparseArrays
using LinearAlgebra
using BenchmarkTools
using DifferentialEquations
using SimulationLogs
using DataFrames
using AeroMDAO

struct StateVector{T <: Real}
    position     :: SVector{3,T}
    euler_angles :: SVector{3,T}
    velocity     :: SVector{3,T}
    rates        :: SVector{3,T}
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
        -sθ    0   1 ]
end

function euler_rates_to_body_rates(θ, φ)
    sφ, cφ = sincos(φ)
    sθ, cθ = sincos(θ)

    [ 1    0    -sθ   ;
      0   cφ  sφ * cθ ;
      0  -sφ  cφ * cθ ]
end

function body_to_earth_rates(θ, φ)
    sφ, cφ = sincos(φ)
    cθ, tθ = cos(θ), tan(θ)

    [ 1  sφ * tθ  cφ * tθ ;
      0     cφ     -sφ    ;
      0  sφ / cθ  cφ / cθ ]
end

body_to_earth_rates(Φs) = body_to_earth_rates(Φs[2], Φs[3])

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
angular_acceleration(M, I, Ω, h)      = M - Ω × (I * Ω + h)

# Equations of motion
function aircraft_eom!(dx, x, p, t)
    φ  = @view x[4:6]
    U  = @view x[7:9]

    acc, rot = p(t, x)

    dx[1:3]    = body_to_earth(φ) * U       # mul!(x0, body_to_earth(φ), U)
    dx[4:6]    = body_to_earth_rates(φ) * U # mul!(x0, body_to_earth_rates(φ), U)
    dx[7:9]    = acc
    dx[10:end] = rot

    dx
end

## Surfaces
function vlm_analysis(aircraft, fs, ρ, ref, S, b, c, print = false)
    data =  solve_case(aircraft, fs;
                       rho_ref   = ρ,
                       r_ref     = ref,
                       area_ref  = S,
                       span_ref  = b,
                       chord_ref = c,
                       print     = print);

    # Get data
    nf_coeffs, ff_coeffs = @views data["Aircraft"][1:2]

    # Filter relevant data
    @views [ ff_coeffs; nf_coeffs[4:6] ]
end

# Aircraft definition
include("../vortex-lattice/aircraft_types.jl")

vanilla   = VanillaAirplane()
wing      = vanilla.surfs[1]
wing_mac  = mean_aerodynamic_center(wing)
wing_plan = plot_wing(wing)

print_info(wing, "Wing")

S, b, c = projected_area(wing), span(wing), mean_aerodynamic_chord(wing)
x_w = [ wing_mac[1], 0, 0 ]

## Meshing
wing, htail, vtail = vanilla.surfs

wing_panels, wing_normals   = panel_wing(wing, 8, 6)
htail_panels, htail_normals = panel_wing(htail, 6, 6)
vtail_panels, vtail_normals = panel_wing(vtail, 6, 5)

# Aircraft assembly
aircraft = Dict(
                "Wing"            => Horseshoe.(wing_panels,  wing_normals),
                "Horizontal Tail" => Horseshoe.(htail_panels, htail_normals),
                "Vertical Tail"   => Horseshoe.(vtail_panels, vtail_normals),
               );

## Evaluate case
ρ       = 1.225
ref     = x_w
V, α, β = 5.0, 0.0, 0.0
Ω       = @SVector [0., 0., 0.]     # Angular velocity
fs      = Freestream(V, α, β, Ω)

res = vlm_analysis(aircraft, fs, ρ, ref, S, b, c, true);

## Differential equations setup
thrust_force(u) = SVector(20u^2 + 10u + 20, 0., 0.)
aero_stability(CFs, CMs, q, S, b, c) = SVector(1., -1, 1) .* CFs .* q * S, CMs .* SVector(b, c, b) .* q * S

function compute_6DOF_dynamics(t, x, p)
    r = @view x[1:3]
    φ = @view x[4:6]
    U = @view x[7:9]
    Ω = @view x[10:end]

    ac, mass, g, I, h, ρ, ref, S, b, c = p

    fs = Freestream(U, Ω)
    q  = dynamic_pressure(ρ, fs.speed)

    # Aerodynamics
    @log t
    @log α, β, Ω     = fs.alpha, fs.beta, fs.omega
    @log aero_coeffs = vlm_analysis(ac, fs, ρ, ref, S, b, c)
    
    CFs      = @view aero_coeffs[1:3]
    CMs      = @view aero_coeffs[4:end]
    F_A, M_A = aero_stability(CFs, CMs, q, S, b, c)

    # Propulsion (Gliding for now)
    F_T      = zeros(3) # thrust_force(fs.speed)

    # Force equilibrium
    F = F_T - F_A

    g_b = earth_to_body(φ) * g

    acc = F / mass - (Ω × U - g_b)
    rot = I^(-1) * (M_A - Ω × (I * Ω + h))

    # free = cartesian_to_freestream(velocity(fs))
    # rates = rate_coefficient(Ω, V, b, c)

    # println("Freestream: $free")
    # println("Gravity: $g_b")
    # println("Rotation: $rates")
    # println("Force coefficients: $CFs")
    # println("Moment coefficients: $CMs")

    acc, rot
end

## Test
mass = 250      # Mass
inertia = mass * @SMatrix [ 1. 0. 0. ;
                            0. 1. 0. ;
                            0. 0. 1. ]  # Inertia matrix

O   = zeros(3, 3)
M   = @SMatrix [ I O O    O    ;
                 O I O    O    ;
                 O O I    O    ;
                 O O O inertia ]

φs  = @SVector [0., 0., 0.] # Euler angles
α   = @SVector [0., 0., 0.] # Angular rates
R_o = @SVector [1., 1., 1.] # Position vector of aircraft, Earth frame
g   = @SVector [0., 0., 1.] # Gravitational acceleration, Earth frame
U   = velocity(fs)          # Freestream velocity

## DifferentialEquations setup
state  = StateVector(R_o, φs, U, Ω)
x0     = (collect ∘ vector)(state)
ps     = (aircraft, mass, g, inertia, zeros(3), ρ, ref, S, b, c)


##
tspan       = (0., 1.0)
run         = ODEFunction(aircraft_eom!, syms = [:xₑ, :yₑ, :zₑ, :φ, :θ, :ψ, :u, :v, :w, :p, :q, :r])
prob        = ODEProblem(run, x0, tspan, (t, x) -> compute_6DOF_dynamics(t, x, ps)) # , mass_matrix = M)


# using ModelingToolkit
# de          = modelingtoolkitize(prob)
# jac         = eval(ModelingToolkit.generate_jacobian(de)[2])
# f           = ODEFunction(run, jac = jac, mass_matrix = M)

# prob_jac = ODEProblem(f, x0, tspan, params)

##
@time sol = solve(prob, reltol = 1e-6, abstol = 1e-6);

##
df   = DataFrame(sol)

##
log  = get_log(sol)

##
using Plots
gr(dpi = 300)
plot(sol, tspan = tspan, layout = (4,3), size = (900, 600))
plot!()

##
plot(log.t, reduce(hcat, log.aero_coeffs)', layout = (2, 3), labels = ["CD" "CY" "CL" "Cl" "Cm" "Cn"])