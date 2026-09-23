## Wing planform optimization with SciML framework
using AeroFuse
using Roots
using LinearAlgebra
using ComponentArrays
using DifferentiationInterface
using Zygote

## Elliptic wing planform prediction test
#==========================================================================================#

include("wing_definition.jl")

## Initial guess
n_vars = 32 # Number of spanwise stations
c = 0.125 # Fixed chord
c_w = LinRange(c, c, n_vars) # Constant distribution
CL_tgt = 1.6 # Target lift coefficient
nc = length(c_w)

wing_init = make_wing(c_w)
Sw = projected_area(wing_init) # Reference area

refs = References(
    speed=10.,
    area=Sw,
    span=span(wing_init),
    chord=mean_aerodynamic_chord(wing_init),
    location=mean_aerodynamic_center(wing_init)
)

# Find angle of attack which matches target CL
α0 = find_zero(1.0, Roots.Secant()) do α
    sys = make_case(α, wing_init, refs)
    CL_tgt - get_forces(sys, wing_init).CL
end

## Initial run
sys = make_case(α0, wing_init, refs)
init = get_forces(sys, wing_init)
print_coefficients(sys)

# Common
function get_res(x, sweep_ratio=0.25, ref=refs)
    α = x[1]
    c_w = @view x[2:end]

    # Setup
    wing_mesh = make_wing(c_w, sweep_ratio)
    system = make_case(α, wing_mesh, ref)

    return system, wing_mesh
end

# Objective
function opt_drag(x, p=nothing)
    sys, mesh = get_res(x)

    # Get forces
    res = get_forces(sys, mesh)

    res.CD
end

# Constraints
function con_all(R, x, p)
    c_w = @view x[2:end]
    sys, mesh = get_res(x)

    _, _, CL = farfield(sys)

    R[1] = CL # Lift coefficient
    R[2] = projected_area(mesh) # Area
    R[3:end] = -diff(c_w) # Chord length differences along span

    # @info "Variables": x
    # @info "Constraints:" R

    return nothing
end

## Initial setup and test
x0 = ComponentVector(alpha=α0, chords=c_w)  # Initial guess
cons = ComponentVector( # Constraints
    CL=0.,
    Sw=0.,
    chords=zeros(nc - 1)
)

CD = opt_drag(x0, nothing)
con_all(cons, x0, nothing)