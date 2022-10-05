## Optimisation tests
using AeroFuse
using Roots
using LinearAlgebra

## Using IPOPT with SNOW
using SNOW # Helps set up optimization problems
using Ipopt # The optimizer

## Elliptic wing planform prediction test
#==========================================================================================#

# Reference values
const refs = References(
    speed     = 10.,
    area      = 0.25,
    span      = 2.0,
    chord     = 0.125,
    location  = [0.03125, 0., 0.]
)

# Geometry
function make_wing(xc, w = 0.25)
    # Planform geometry
    n_vars = length(xc)
    bs = fill(1., n_vars - 1)
    wing = Wing(
        foils     = fill(naca4(0,0,1,2), n_vars),
        chords    = xc, # Design variables
        spans     = bs / sum(bs), # Normalizing halfspan to 1
        w_sweep   = w, # Quarter-chord
        symmetry  = true
    )

    # Meshing
    wing_mesh = WingMesh(
        wing, fill(2, n_vars - 1), 1, 
        span_spacing = Uniform(),
    );

    return wing_mesh
end

# VLM analysis
function make_case(α, wing_mesh, refs)
    aircraft = ComponentVector(wing = make_horseshoes(wing_mesh))

    # Freestream conditions
    fs = Freestream(alpha = α)  # Design variable: Angle of attack

    # Solve system
    return solve_case(aircraft, fs, refs)
end

# Aerodynamic forces
function get_forces(system, wing_mesh)
    # Evaluate aerodynamic coefficients
    CDi, _, CL = farfield(system)

    # Calculate equivalent flat-plate skin-friction drag
    # CDv = profile_drag_coefficient(wing_mesh, 1.0, system.reference)

    # Calculate local-dissipation/local-friction drag
    CVs = norm.(surface_velocities(system)).wing
    CDv = profile_drag_coefficient(wing_mesh, 0.98, CVs, system.reference)

    return (CDi = CDi, CDv = CDv, CD = CDi + CDv, CL = CL)
end

## Objective and constraints
function optimize_drag!(g, x, w = 0.25)
    α = x[1]
    xc = @view x[2:end]
    wing_mesh = make_wing(xc, w)
    system = make_case(α, wing_mesh, refs)

    res = get_forces(system, wing_mesh)

    # Objective
    f = res.CD

    # Constraints
    g[1] = res.CL # Lift coefficient
    g[2] = projected_area(wing_mesh) # Area
    g[3:end] = @views xc[1:end-1] - xc[2:end] # Chord length differences along span

    return f
end

## Initial guess
n_vars = 60 # Number of spanwise stations
c = 0.125 # Fixed chord
xs = LinRange(c, c, n_vars) # Constant distribution
CL_tgt = 1.2 # Target lift coefficient

wing_init = make_wing(xs)

# Find angle of attack which matches target CL
α0 = find_zero(1.0, Roots.Secant()) do α
    sys = make_case(α, wing_init, refs)
    CL_tgt - get_forces(sys, wing_init).CL
end

## Initial setup and test
x0 = ComponentVector(alpha = α0, chords = xs)  # Initial guess
ng = length(x0) # Number of constraints
cons = zeros(ng)

CD = optimize_drag!(cons, x0)

# Bounds
lx = [ -Inf; zero(xs) ]
ux = [ Inf; Inf * ones(length(xs)) ] 
lg = [ CL_tgt; Sw; zeros(ng - 2) ]
ug = [ CL_tgt; Sw; Inf * zeros(ng - 2) ]

# Sparsity
# sp = SparsePattern(ForwardAD(), optimize_drag!, ng, lx, ux)
opts = Options(
    solver = IPOPT(),
    derivatives = ForwardAD(),
    # sparsity = sp
)

## Run optimization
@time xopt, fopt, info = minimize(optimize_drag!, x0, ng, lx, ux, lg, ug, opts)

## Substitute
wing_opt = make_wing(xopt[2:end])
sys_opt = make_case(xopt[1], wing_opt, refs)
print_coefficients(sys_opt)
opt = get_forces(sys_opt, wing_opt)

## Comparison with initial
Sw = projected_area(wing_init)
sys = make_case(xopt[1], wing_init, refs)
print_coefficients(sys)
init = get_forces(sys, wing_init)

## Plotting
#==========================================================================================#

# Exact solution
ys = LinRange(0., 1., n_vars)
xs = @. √(1. - ys^2) * 4 / π * c

wing_exact = make_wing(xs)
sys_exact = make_case(xopt[1], wing_exact, refs)
print_coefficients(sys_exact)
exact = get_forces(sys_exact, wing_exact)

## Plot spanwise loading
ll_init = spanwise_loading(wing_init, surface_coefficients(sys)[1].wing, sys.reference.area)
ll_opt = spanwise_loading(wing_opt, surface_coefficients(sys_opt)[1].wing, sys_opt.reference.area)

CL_init = vec(sum(sys.circulations.wing, dims = 1)) / (0.5 * sys.reference.speed * sys.reference.chord)
CL_loads = vec(sum(sys_opt.circulations.wing, dims = 1)) / (0.5 * sys_opt.reference.speed * sys_opt.reference.chord)

ll_exact = spanwise_loading(wing_exact, surface_coefficients(sys_exact)[1].wing, sys_exact.reference.area)
CL_exact = vec(sum(sys_exact.circulations.wing, dims = 1)) / (0.5 * sys_exact.reference.speed * sys_exact.reference.chord)

##
using Plots, LaTeXStrings

pgfplotsx() # Needs LaTeX
# gr()
# plotlyjs()

plt_opt = plot(
    camera = (90, 90),
    legend = :bottom,
    xlabel = L"x,~m",
    guidefontrotation = 90.0,
    title = LaTeXString("Planform, \$ S = $(round(Sw; digits = 4)),~C_L = $(round(CL_tgt; digits = 4)) \$"),
    grid = false,
    # aspect_ratio = 1,
    # zlim = (-0.5, 0.5) .* span(wing_init)
)

# Initial planform
plot!(wing_init.surface,
    lc = :black, 
    mc = :black,
    lw = 0.8, 
    # alpha = 0.5, 
    label = "",
)
# plot!(sys, wing_init.surface, 
#     span = 10,
#     lw = 0.3,
#     alpha = 0.5,
#     lc = :grey
# )

# Optimized planform
plot!(wing_opt.surface, 
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
    label = LaTeXString("Initial Wing: \$ (C_{D_i}, C_{D_v}, C_D) = $(round.([init.CDi;init.CDv;init.CD]; digits = 4)) \$"),
)
plot!(
    [ -cumsum(wing_opt.surface.spans)[end:-1:1]; 0; cumsum(wing_opt.surface.spans) ], 
    [ wing_opt.surface.chords[end:-1:2]; wing_opt.surface.chords ],
    lc = :cornflowerblue, 
    label = LaTeXString("Optimized Wing: \$ (C_{D_i}, C_{D_v}, C_D) = $(round.([opt.CDi;opt.CDv;init.CD]; digits = 4)) \$"),
)


plot!(ys, xs, label = LaTeXString("Inviscid Optimum: \$ (C_{D_i}, C_{D_v}, C_D) = $(round.([exact.CDi;exact.CDv;exact.CD]; digits = 4)) \$"))

plt_CL = plot(
    title = "Lift Distribution",
    ylabel = L"C_L",
    xlabel = L"Spanwise Location $y$, $m$",
    grid = false,
)
plot!(ll_init[:,1], CL_init,
    lc = :black, label = ""
)
plot!(ll_opt[:,1], CL_loads, 
    lc = :cornflowerblue, label = ""
)
plot!(ll_exact[:,1], CL_exact, 
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

#
savefig(plt_wing, "plots/WingOptimization.pdf")