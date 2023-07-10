## Wing planform optimization with SciML framework
using AeroFuse
using Roots
using LinearAlgebra
using Optimization, OptimizationMOI
using Ipopt

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
    speed     = 10.,
    area      = Sw,
    span      = span(wing_init),
    chord     = mean_aerodynamic_chord(wing_init),
    location  = mean_aerodynamic_center(wing_init)
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
function get_res(x, chord_ratio = 0.25, ref=  refs)
    α = x[1]
    c_w = @view x[2:end]

    # Setup
    wing_mesh = make_wing(c_w, chord_ratio)
    system = make_case(α, wing_mesh, ref)

    return system, wing_mesh
end

# Objective
function opt_drag(x, p = nothing) 
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
x0 = ComponentVector(alpha = α0, chords = c_w)  # Initial guess
cons = ComponentVector( # Constraints
    CL = 0.,
    Sw = 0.,
    chords = zeros(nc - 1)
)

CD = opt_drag(x0, nothing)
con_all(cons, x0, nothing)

# Bounds
lx = [ -Inf; zeros(nc) ]
ux = [ Inf; Inf * ones(nc) ] 
lg = [ CL_tgt; Sw; zeros(nc - 1) ]
ug = [ CL_tgt; Sw; Inf * ones(nc - 1) ]

ng = length(x0) # Number of constraints

## Problem construction
optprob = OptimizationFunction(opt_drag, Optimization.AutoForwardDiff(), cons = con_all)
prob = OptimizationProblem(optprob, x0, lcons = lg, ucons = ug, lb = lx, ub = ux)

## Choose optimizer
opt = OptimizationMOI.MOI.OptimizerWithAttributes(
    # NLopt.Optimizer,
    # "algorithm" => :LN_NELDERMEAD,
    Ipopt.Optimizer, 
    "hessian_approximation" => "limited-memory"
)

# opt = IPNewton()

## Solve
@time sol = solve(prob, opt);

## Substitute
xopt = sol.u
wing_opt = make_wing(xopt[2:end])
sys_opt = make_case(xopt[1], wing_opt, refs)
opt = get_forces(sys_opt, wing_opt)
print_coefficients(sys_opt)

# Exact solution
y_exact = LinRange(0., 1., n_vars)
x_exact = @. √(1. - y_exact^2) * 4 / π * c

wing_exact = make_wing(x_exact)
sys_exact = make_case(xopt[1], wing_exact, refs)
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
# gr()
# plotlyjs()
##
plt_opt = plot(
    camera = (90, 90),
    legend = :bottom,
    xlabel = L"x,~m",
    guidefontrotation = 90.0,
    title = LaTeXString("Planform, \$ S = $(round(Sw; digits = 4)),~C_{L_{req}} = $(round(CL_tgt; digits = 4)) \$"),
    grid = false,
    # aspect_ratio = 1,
    # zlim = (-0.5, 0.5) .* span(wing_init)
)

# Initial planform
plot!(wing_init.surface,
    lc = :black, 
    mc = :black,
    lw = 0.8, 
    alpha = 0.6, 
    label = "",
)
# plot!(sys, wing_init.surface, 
#     span = 10,
#     lw = 0.3,
#     alpha = 0.5,
#     lc = :grey
# )

# Optimized planform
plot!(wing_opt, 
    lc = :cornflowerblue, 
    mc = :cornflowerblue, 
    lw = 0.8, 
    label = "",
)
# plot!(sys_opt, wing_init.surface, 
#     span = 10, 
#     lw = 0.1,
#     lc = :cornflowerblue,
#     alpha = 0.5,
#     label = ""
# )

# Exact solution planform
plot!(wing_exact.surface, 
    lc = :green, 
    mc = :green, 
    alpha = 0.6,
    lw = 0.8, 
    label = "",
)

# Planform distribution
plt_plan = plot(
    ylabel = L"Chord Length $c$, $m$", 
    legend = :bottom,
    title = "Chord Distribution",
    grid = false
)

plot!(
    [ -cumsum(wing_init.surface.spans)[end:-1:1]; 0; cumsum(wing_init.surface.spans) ], 
    [ wing_init.surface.chords[end:-1:2]; wing_init.surface.chords ], 
    lc = :black,
    label = LaTeXString("Initial Wing: \$ (C_{D_i}, C_{D_v}, C_D, C_L) = $(round.([init.CDi;init.CDv;init.CD;init.CL]; digits = 4)) \$"),
)
plot!(
    [ -cumsum(wing_opt.surface.spans)[end:-1:1]; 0; cumsum(wing_opt.surface.spans) ], 
    [ wing_opt.surface.chords[end:-1:2]; wing_opt.surface.chords ],
    lc = :cornflowerblue, 
    label = LaTeXString("Optimized Wing: \$ (C_{D_i}, C_{D_v}, C_D, C_L) = $(round.([opt.CDi;opt.CDv;opt.CD;opt.CL]; digits = 4)) \$"),
)
plot!(
    [ -cumsum(wing_exact.surface.spans)[end:-1:1]; 0; cumsum(wing_exact.surface.spans) ], 
    [ wing_exact.surface.chords[end:-1:2]; wing_exact.surface.chords ],
    lc = :green, 
    label = LaTeXString("Inviscid Optimum: \$ (C_{D_i}, C_{D_v}, C_D, C_L) = $(round.([exact.CDi;exact.CDv;exact.CD;exact.CL]; digits = 6)) \$"),
)

plt_CL = plot(
    title = "Lift Distribution",
    ylabel = L"C_L",
    xlabel = L"Spanwise Location $y$, $m$",
    grid = false,
)
plot!(ll_init[:,1], ll_init[:,5],
    lc = :black, label = ""
)
plot!(ll_opt[:,1], ll_opt[:,5], 
    lc = :cornflowerblue, label = ""
)
plot!(ll_exact[:,1], ll_exact[:,5], 
    lc = :green, label = ""
)

l = @layout [ a; b; c ]
plt_wing = plot(
    plt_opt, 
    plt_plan, 
    plt_CL, 
    layout = l, 
    size = (700,700)
)

##
savefig(plt_wing, "plots/SciMLWingOptimization.pdf")