## Wing planform optimization with SciML framework using B-splines
using AeroFuse
using Roots
using LinearAlgebra
using Optimization, OptimizationIpopt
using ComponentArrays
using BSplineKit

import Plots: plot, plot!, scatter!, savefig

## Elliptic wing planform prediction test
#==========================================================================================#

include("wing_definition.jl")

## Initial guess
n_vars = 64 # Number of spanwise evaluation stations
n_cp = 6 # Number of control points for B-spline
c = 0.125 # Fixed chord initial guess
c_cp_init = fill(c, n_cp) # Constant distribution for control points
CL_tgt = 1.6 # Target lift coefficient

# B-spline parameterization
y_cp = LinRange(0, 1, n_cp) # Normalized spanwise control points
y_eval = LinRange(0, 1, n_vars) # Normalized evaluation stations

function evaluate_chords(c_cp, y_cp, y_eval)
    itp = BSplineKit.interpolate(y_cp, c_cp, BSplineOrder(4))
    return itp.(y_eval)
end

c_w = evaluate_chords(c_cp_init, y_cp, y_eval)
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
    c_cp = @view x[2:end]

    # Evaluate chords using B-Spline
    c_w = evaluate_chords(c_cp, y_cp, y_eval)

    # Setup
    wing_mesh = make_wing(c_w, sweep_ratio)
    system = make_case(α, wing_mesh, ref)

    return system, wing_mesh, c_w
end

# Objective
function opt_drag(x, p=nothing)
    sys, mesh, c_w = get_res(x)

    # Get forces
    res = get_forces(sys, mesh)

    res.CD
end

# Constraints
function con_all(R, x, p)
    sys, mesh, c_w = get_res(x)

    _, _, CL = farfield(sys)

    R[1] = CL # Lift coefficient
    R[2] = projected_area(mesh) # Area

    # Apply monotonicity constraint directly to control points
    c_cp = @view x[2:end]
    R[3:end] = -diff(c_cp) # Control point differences along span must be negative/zero

    return nothing
end

## Initial setup and test
x0 = ComponentVector(alpha=α0, chords=c_cp_init)  # Initial guess
cons = ComponentVector( # Constraints
    CL=0.,
    Sw=0.,
    chords=zeros(n_cp - 1)
)

CD = opt_drag(x0, nothing)
con_all(cons, x0, nothing)

# Bounds for control points
lx = [-Inf; zeros(n_cp)]
ux = [Inf; Inf * ones(n_cp)]
lg = [CL_tgt; Sw; zeros(n_cp - 1)]
ug = [CL_tgt; Sw; Inf * ones(n_cp - 1)]

ng = length(x0) # Number of constraints

## Problem construction
optprob = OptimizationFunction(opt_drag, Optimization.AutoForwardDiff(), cons=con_all)
prob = OptimizationProblem(optprob, x0[:], lcons=lg, ucons=ug, lb=lx, ub=ux)

## Choose optimizer
opt = IpoptOptimizer(;
    hessian_approximation="limited-memory"
)

## Solve
@time sol = solve(prob, opt; verbose=true);

## Substitute
xopt = ComponentArray(sol.u, getaxes(x0))
c_w_opt = evaluate_chords(xopt.chords, y_cp, y_eval)
wing_opt = make_wing(c_w_opt)
sys_opt = make_case(xopt.alpha, wing_opt, refs)
opt_forces = get_forces(sys_opt, wing_opt)
print_coefficients(sys_opt)

# Exact solution
y_exact = LinRange(0., 1., n_vars)
x_exact = @. √(1. - y_exact^2) * 4 / π * c

wing_exact = make_wing(x_exact)
sys_exact = make_case(xopt.alpha, wing_exact, refs)
exact = get_forces(sys_exact, wing_exact)
print_coefficients(sys_exact)

## Plotting
#==========================================================================================#

## Plot spanwise loading
ll_init = spanwise_loading(wing_init, sys.reference, surface_coefficients(sys)[1].wing, sys.circulations.wing)
ll_opt = spanwise_loading(wing_opt, sys_opt.reference, surface_coefficients(sys_opt)[1].wing, sys_opt.circulations.wing)
ll_exact = spanwise_loading(wing_exact, sys_exact.reference, surface_coefficients(sys_exact)[1].wing, sys_exact.circulations.wing)

##
using Plots, LaTeXStrings

pgfplotsx() # Needs LaTeX

plt_opt = plot(
    camera=(90, 90),
    legend=:bottom,
    xlabel=L"x,~m",
    guidefontrotation=90.0,
    title=LaTeXString("Planform, \$ S = $(round(Sw; digits = 4)),~C_{L_{req}} = $(round(CL_tgt; digits = 4)) \$"),
    grid=false,
)

# Initial planform
plot!(wing_init.surface,
    lc=:black,
    mc=:black,
    lw=0.8,
    alpha=0.6,
    label="",
)
# Optimized planform
plot!(wing_opt.surface,
    lc=:cornflowerblue,
    mc=:cornflowerblue,
    lw=0.8,
    label="",
)

# Exact solution planform
plot!(wing_exact.surface,
    lc=:green,
    mc=:green,
    alpha=0.6,
    lw=0.8,
    label="",
)

# Planform distribution
plt_plan = plot(
    ylabel=L"Chord Length $c$, $m$",
    legend=:bottom,
    title="Chord Distribution",
    grid=false
)

plot!(
    [-cumsum(wing_init.surface.spans)[end:-1:1]; 0; cumsum(wing_init.surface.spans)],
    [wing_init.surface.chords[end:-1:2]; wing_init.surface.chords],
    lc=:black,
    label=LaTeXString("Initial Wing: \$ (C_{D_i}, C_{D_v}, C_D, C_L) = $(round.([init.CDi;init.CDv;init.CD;init.CL]; digits = 4)) \$"),
)
plot!(
    [-cumsum(wing_opt.surface.spans)[end:-1:1]; 0; cumsum(wing_opt.surface.spans)],
    [wing_opt.surface.chords[end:-1:2]; wing_opt.surface.chords],
    lc=:cornflowerblue,
    label=LaTeXString("Optimized Wing: \$ (C_{D_i}, C_{D_v}, C_D, C_L) = $(round.([opt_forces.CDi;opt_forces.CDv;opt_forces.CD;opt_forces.CL]; digits = 4)) \$"),
)
plot!(
    [-cumsum(wing_exact.surface.spans)[end:-1:1]; 0; cumsum(wing_exact.surface.spans)],
    [wing_exact.surface.chords[end:-1:2]; wing_exact.surface.chords],
    lc=:green,
    label=LaTeXString("Inviscid Optimum: \$ (C_{D_i}, C_{D_v}, C_D, C_L) = $(round.([exact.CDi;exact.CDv;exact.CD;exact.CL]; digits = 6)) \$"),
)

plt_CL = plot(
    title="Lift Distribution",
    ylabel=L"C_L",
    xlabel=L"Spanwise Location $y$, $m$",
    grid=false,
)
plot!(ll_init[:, 1], ll_init[:, 5],
    lc=:black, label=""
)
plot!(ll_opt[:, 1], ll_opt[:, 5],
    lc=:cornflowerblue, label=""
)
plot!(ll_exact[:, 1], ll_exact[:, 5],
    lc=:green, label=""
)

l = @layout [a; b; c]
plt_wing = plot(
    plt_opt,
    plt_plan,
    plt_CL,
    layout=l,
    size=(700, 700)
)

mkpath("plots")
savefig(plt_wing, "plots/SciMLBSplineWingOptimization.pdf")
