## VLM inverse analysis by nonlinear solution
using AeroFuse
using LinearAlgebra

## Surfaces
#================================================#

function make_surfaces()
    ## Wing
    wing = Wing(
        foils     = fill(naca4(2,4,1,2), 2),
        chords    = [1.0, 0.6],
        twists    = [2.0, 0.0],
        spans     = [4.0],
        dihedrals = [5.],
        sweeps    = [5.],
        symmetry  = true
    )

    ## Horizontal tail
    htail = Wing(
        foils     = fill(naca4(0,0,1,2), 2),
        chords    = [0.7, 0.42],
        twists    = [0.0, 0.0],
        spans     = [1.25],
        dihedrals = [0.],
        sweeps    = [6.39],
        position  = [4., 0, -0.1],
        angle     = -2.,
        axis      = [0., 1., 0.],
        symmetry  = true
    )

    ## Vertical tail
    vtail = Wing(
        foils     = fill(naca4(0,0,0,9), 2),
        chords    = [0.7, 0.42],
        twists    = [0.0, 0.0],
        spans     = [1.0],
        dihedrals = [0.],
        sweeps    = [7.97],
        position  = [4., 0, 0],
        angle     = 90.,
        axis      = [1., 0., 0.]
    )

    ## Meshing
    #================================================#

    wing_mesh  = WingMesh(wing, [24], 12, span_spacing = Cosine())
    htail_mesh = WingMesh(htail, [12], 6, span_spacing = Cosine())
    vtail_mesh = WingMesh(vtail, [12], 6, span_spacing = Cosine())

    return wing_mesh, htail_mesh, vtail_mesh
end

wing, htail, vtail = make_surfaces()

# Assemble
aircraft = ComponentArray(
    wing  = make_horseshoes(wing),
    htail = make_horseshoes(htail),
    vtail = make_horseshoes(vtail)
);

## Plot
using Plots
# plotlyjs(dpi = 300)

p1 = plot(
    aspect_ratio = 1,
    zlim = (-0.5, 0.5) .* span(wing)
)

plot!(wing, label = "Wing")
plot!(htail, label = "Horizontal Tail")
plot!(vtail, label = "Vertical Tail")

## VLM analysis
#================================================#

fs  = Freestream(
    alpha = 0.0, 
    beta  = 0.0, 
    omega = [0., 0., 0.]
)

ref = References(
    speed     = 150.0, 
    density   = 1.225,
    viscosity = 1.5e-5,
    area      = projected_area(wing),
    span      = span(wing),
    chord     = mean_aerodynamic_chord(wing),
    location  = mean_aerodynamic_center(wing)
)

# Run initial for compilation
@time sys = solve_case(aircraft, fs, ref)

print_coefficients(sys)


## Nonlinear analyses
#================================================#

## Specified CY and CL
function solve_angles!(R, x, CL_ref = 0.5, CY_ref = -0.01)
    fs_n = Freestream(alpha = x[1], beta = x[2])
    sys = solve_case(aircraft, fs_n, ref)
    _, CY, CL = farfield(sys)

    R[1] = CL - CL_ref
    R[2] = CY - CY_ref

    return R
end

# Initial guess
x0 = ComponentArray(
    alpha = 1.,
    beta = 1.
)

## Check derivative functionality
using ForwardDiff

R = similar(x0)
R_x = R .* R'

ForwardDiff.jacobian!(R_x, solve_angles!, R, x0)

## Solve nonlinear system
using NLsolve

res = nlsolve(
    solve_angles!,
    x0,
    method = :newton,
    autodiff = :forward,
    show_trace = true
)

## Check results
fs_n = Freestream(
    alpha = res.zero.alpha, 
    beta = res.zero.beta
)

sys = solve_case(aircraft, fs_n, ref)

print_coefficients(sys)

## Plot
plot!(p1, sys, wing.surface, dist = 10.)
plot!(p1, sys, htail.surface, dist = 2.)
plot!(p1, sys, vtail.surface, dist = 2.)