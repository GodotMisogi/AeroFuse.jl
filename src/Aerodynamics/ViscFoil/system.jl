## NEWTON SYSTEM SETUP
#============================================#

function boundary_layer_finite_difference(Us, δss, θs, Ts, ΔXs, panels, ñ_crit, a, ν)

    R = zeros(length(panels), 3)
    evaluate_residuals!(R, bl_sys)


    R
end

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

    ## Linear vorticity-stream panel method schema
    

    @show U_invs

    # Ugly tagging and initialization
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
