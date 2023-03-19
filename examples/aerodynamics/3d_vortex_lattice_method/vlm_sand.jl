##
using AeroFuse
using JuMP
using Ipopt
using NLsolve

## Wing section setup
wing = Wing(
    foils     = fill(naca4((2,4,1,2)), 3),
    chords    = [2.0, 1.6, 0.2],
    twists    = [0.0, 0.0, 0.0],
    spans     = [5., 0.6],
    dihedrals = [5., 5.],
    sweeps    = [20.,20.],
    w_sweep   = 0.25, # Quarter-chord sweep
    symmetry  = true,
    # flip      = true
)

##
x_w, y_w, z_w = wing_mac = mean_aerodynamic_center(wing)
print_info(wing, "Wing")

## Meshing and assembly
wing_mesh = WingMesh(wing, [24,12], 6, 
                     span_spacing = Cosine()
                    );
aircraft = ComponentVector(wing = make_horseshoes(wing_mesh))

# Freestream conditions
fs  = Freestream(
    alpha = 2.0, # deg
    beta  = 2.0, # deg
    omega = [0.,0.,0.]
)

# Reference values
ref = References(
    speed     = 150, # m/s
    density   = 1.225, # kg/m³
    viscosity = 1.5e-5, # ???
    area      = projected_area(wing), # m²
    span      = span(wing), # m
    chord     = mean_aerodynamic_chord(wing), # m
    location  = mean_aerodynamic_center(wing) # m
)

## Solve system
system = VortexLatticeSystem(aircraft, fs, ref)

## NLsolve
R = similar(system.circulations)

##
AeroFuse.solve_nonlinear!(R, x) = solve_nonlinear!(R, system.vortices, x ./ ref.speed, velocity(fs) / ref.speed, fs.omega)

##
xz = nlsolve(
    solve_nonlinear!,
    system.circulations,
    autodiff = :forward,
    show_trace = true,
)

n = length(system.circulations)

## THIS IS THE IMPLICIT FUNCTION THEOREM APPLIED, R(x) = R(x, u(x))
obtain_root(x) = nlsolve(u -> solve_nonlinear!(R, x, u), system.circulations).zero # Rootfinding method (doesn't matter which)


## ChainRules adjoint
## Evaluate
x0 = [1.2, -2.0]
obtain_root(x0)
loss_functional(obtain_root(x0))

## Define adjoint rule
using ChainRulesCore

function ChainRulesCore.rrule(::typeof(obtain_root), x)
    # Solve residual equation 
    # R(x,u) =(by implicit function theorem)= u -> R(u(x)) = 0
    u_class = obtain_root(x)

    # Define pullback
    function obtain_root_pullback(∂J_∂u)
        # Evaluate gradients ∂R/∂x, ∂R/∂u
        ∂R_∂x, ∂R_∂u = ForwardDiff.jacobian!(residual, x, u_class)

        # Adjoint variable:
        # λᵀ = ∂J/∂R = (∂R/∂u)^(-1) * ∂J/∂u
        λᵀ = ∂R_∂u' \ ∂J_∂u

        # Loss sensitivity to input:
        # ∂J/∂x = -λ^T * ∂R/∂x
        dJ_dx = -λᵀ' * ∂R_∂x

        # Loss sensitivity to residual equation:
        # ∂J/∂R = 0
        ∂J_∂R = NoTangent()

        return ∂J_∂R, dJ_dx'
    end

    return u_class, obtain_root_pullback
end




## Create JuMP model, using Ipopt as the solver
model = Model(optimizer_with_attributes(Ipopt.Optimizer))

@variables(model, begin
    Γs[j = 1:n] == system.circulations[j]
end)

register(model, :R_A, n, solve_nonlinear; autodiff = true)

# @NLconstraint(model, [j = 1:n], solve_nonlinear(Γs)[i] == 0)

##
optimize!(model)