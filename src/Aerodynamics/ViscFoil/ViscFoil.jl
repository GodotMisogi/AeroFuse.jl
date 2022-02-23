module ViscFoil

import ..MathTools: midpair_map, forward_difference, forward_sum, weighted_vector

import ..Laplace: Uniform2D, velocity

import ..PanelGeometry: AbstractPanel, AbstractPanel2D, Panel2D, WakePanel2D, distance, panel_length, panel_angle, tangent_vector, panel_points, panel_location, collocation_point

import ..DoubletSource: doublet_matrix, source_matrix, source_strengths, boundary_vector, solve_linear, lift_coefficient

import ..LinearVortexSource: vortex_influence_matrix, source_influence_matrix, linear_vortex_matrix, linear_source_matrix, neumann_boundary_condition

using NLsolve
using LineSearches: BackTracking
using ForwardDiff
using LinearAlgebra

include("closure_relations.jl")

reynolds_number(U, L, ν)    = U * L / ν
δ_star(m, U_e)              = m / U_e             # δ* = m/U_e
shape_parameter(m, U_e, θ)  = δ_star(m, U_e) / θ  # H = δ*/θ

include("boundary_layer.jl")

## Inviscid doublet-source solution and edge velocities
#============================================#

function solve_inviscid_doublets(panels, wakes, u)
    # Inviscid solution
    all_panels      = [ panels; wakes ]
    foil_doublets   = doublet_matrix(panels, panels)
    all_sources     = source_matrix(panels, all_panels)
    A_inv           = foil_doublets^(-1)
    P               = -A_inv * all_sources

    # Mass defect block computation
    all_lengths = panel_length.(all_panels)         #  N + N_w
    B           = source_matrix(wakes, all_panels)  # (N_w, N + N_w)
    A_w         = doublet_matrix(wakes, panels)     # (N_w, N)
    dKds        = -midpair_map(-, defect_block(P, A_w, B)) ./ all_lengths
    D           = dKds ./ all_lengths'

    # Inviscid velocities
    φs          = A_inv * boundary_vector(panels, u, (0., 0.))
    q0s, Δrs    = tangential_velocities(panels, φs, u, false)
    U_invs      = [ q0s; Ref(u) .+ tangent_vector.(wakes) ]

    U_invs, D, all_lengths
end

# Inviscid linear vortex solution and edge velocities
#============================================#

function solve_inviscid_vortices(panels, wakes, u)
    # Inviscid solution
    all_panels      = [ panels; wakes ]
    foil_vortices   = vortex_influence_matrix(panels)
    all_sources     = source_influence_matrix(panels, all_panels)
    A_inv           = foil_vortices^(-1)
    P               = -A_inv * all_sources

    # Mass defect block computation
    all_lengths     = panel_length.(all_panels)                 # Foil and Wake: N + N_w
    B               = linear_source_matrix(wakes, all_panels)
    A_w             = linear_vortex_matrix(wakes, panels)
    C               = defect_block(P, A_w, B)[:,1:end-1]
    D               = C ./ all_lengths'

    # Inviscid velocities
    U_inv_f     = A_inv * neumann_boundary_condition(panels, u) # Foil: N + 1
    U_inv_w     = neumann_boundary_condition(wakes, u)          # Wake: N_w + 1
    U_invs      = [ U_inv_f[1:end-1]; U_inv_w ]                 # Foil and wake: N + N_w + 2

    U_invs, D, all_lengths
end

## Defect formulation
#============================================#

defect_block(P, Aw, B)         = [ P; Aw * P + B ]
edge_velocities(ms, D, U_invs) = U_invs + D * ms

## Newton setup
#============================================#

function solve_system(x_in, panels, all_lengths, C, U_invs, n_crit)
    # Number of nodes
    n  = length(all_lengths)

    # Segregating into inputs
    ms = x_in[1:n+1] # N + N_w + 1
    θs = x_in[n+2:2n+1] # N + N_w
    Ts = x_in[2n+2:end] # N + N_w

    # Forward differencing mass defects
    Δms  = forward_difference(ms) # N + N_w + 1 -> N + N_w

    # println(size(U_invs), size(ms), size(θs), size(Ts))
    # println(size(Δms), size(all_lengths))

    # Compute edge velocities
    σs      = Δms ./ all_lengths
    Us      = abs.(edge_velocities(σs, C, U_invs))
    δ_stars = @. abs(ms[1:end-1] / Us[1:end-1])     # N + Nw
    ΔXs     = all_lengths

    # Constants
    a = 330
    ν = 1.5e-5

    # Trailing edge initial conditions for θ and δ*

    # Initial wake shear coefficient
    # c_τ_wake = c_τ

    # Evaluate residuals
    # R1, R2, R3 = 
    
    R = boundary_layer_finite_difference(Us[1:end-1], δ_stars, θs, Ts, ΔXs, panels, n_crit, a, ν)[:]

    # [ R1; R2; R3 ]
end

function solve_viscous_case(panels, wakes, uniform :: Uniform2D)
    u = velocity(uniform)
    α = uniform.angle

    # Evaluate inviscid edge velocities and mass defect matrix
    #======================================================#

    ## Doublet-source panel method schema
    # U_invs, D, all_lengths = solve_inviscid_doublets(panels, wakes, u)

    ## Linear vortex-source panel method schema
    U_invs, D, all_lengths = solve_inviscid_vortices(panels, wakes, u)

    @show U_invs

    # Ugly tagging and setup
    #======================================================#

    # tags    = [ fill(typeof(panels[1]), length(panels)) ; # Foil panels
                # fill(typeof(wakes[1] ), length(wakes) ) ] # Wake panels
    n_crit  = 9 # Critical amplification ratio for Tollmien-Schlichting waves to transition

    # Initialize variables
    num_nodes   = length(panels) + length(wakes)                        # N + N_w
    ms          = U_invs ./ [ all_lengths; all_lengths[end] ] .* 0.01   # N + N_w + 1
    θs          = fill(0.1, num_nodes)                                  # N + N_w
    ns          = fill(6.0, num_nodes)                                  # N + N_w
    x0          = [ ms; θs; ns ] # Is there any better way to guess than the heuristic?

    # Solving system
    #======================================================#

    # R  = zeros(num_nodes * num_vars + 1)
    system(x)   = solve_system(x, [panels; wakes], all_lengths, D, U_invs, n_crit)
    state       = nlsolve(system,
                          x0,
                          iterations = 20,
                        #   autodiff = :forward,
                          show_trace = true,
                          # extended_trace = true,
                          store_trace = true,
                          method = :newton,
                          linesearch = BackTracking())
end

end

# if writer
#  open("state_variables.txt", "w") do io
#  write(io, "Iteration $n \n")
#  write(io, "Solving system... \n")
#  write(io, "ms: $ms \n")
#  write(io, "θs: $θs \n")
#  write(io, "Ue/V∞: $(ues ./ norm(u)) \n")
#  write(io, "δ*: $δ_stars \n")
#  write(io, "H: $Hs \n")
#  write(io, "H*: $H_stars \n")
#  write(io, "Cfs: $cfs \n")
#  end;
# end

# NOTE: central differencing has different endpoints for all_panels
# foil_lengths = panel_length.(panels)  # N  -> N
# dP  = -midpair_map(-, P) # (N, N + N_w)  -> (N, N + N_w)

# Constructing forward differences for mass defect
# foil_lengths = distance.(panels[1:end-1], panels[2:end])  # N  -> N - 1
# all_lengths  = distance.(all_panels[1:end-1], all_panels[2:end])  # N + N_w  -> N + N_w - 1
# dP   = P[2:end,:] .- P[1:end-1,:]  # (N, N + N_w)  -> (N - 1, N + N_w)
