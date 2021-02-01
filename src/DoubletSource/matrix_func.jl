"""
doublet_matrix(panels_1, panels_2)

Creates the matrix of doublet potential influence coefficients between pairs of panels_1 and panels_2.
"""
doublet_matrix(panels_1 :: AbstractVector{<: Panel2D}, panels_2 :: AbstractVector{<: Panel2D}) = map(x -> ifelse(x[1] === x[2], 0.5, doublet_influence(x[2], x[1])), Iterators.product(panels_1, panels_2))

"""
    kutta_condition(panels)

Creates the vector describing Morino's Kutta condition given Panel2Ds.
"""
kutta_condition(panels :: AbstractVector{<: Panel2D}) = [ 1; -1; zeros(length(panels) - 4); 1; -1 ]

"""
    wake_vector(wake_panel, panels)

Creates the vector of doublet potential influence coefficients from the wake on the panels given the wake panel and the array of Panel2Ds.
"""
wake_vector(woke_panel :: Panel2D, panels :: AbstractVector{<: Panel2D}) = doublet_influence.(Ref(woke_panel), panels)


"""
    influence_matrix(panels, wake_panel)

Assembles the Aerodynamic Influence Coefficient matrix consisting of the doublet matrix, wake vector, Kutta condition given Panel2Ds and the wake panel.
"""
influence_matrix(panels :: AbstractVector{<: Panel2D}, woke_panel :: Panel2D) = 
    vcat(hcat(doublet_matrix(panels, panels), wake_vector(woke_panel, panels)), hcat(kutta_condition(panels)', 0))
                            #  0                ]

"""
    source_matrix(panels_1, panels_2)

Creates the matrix of source potential influence coefficients between pairs of `panels_1` and `panels_2`.
"""
source_matrix(panels_1 :: AbstractVector{<: Panel2D}, panels_2 :: AbstractVector{<: Panel2D}) = [ source_influence(panel_j, panel_i) for panel_i ∈ panels_1, panel_j ∈ panels_2 ]

"""
    source_strengths(panels, freestream)

Creates the vector of source strengths for the Dirichlet boundary condition ``\\sigma = \\vec U_{\\infty} \\cdot \\hat{n}`` given Panel2Ds and a Uniform2D.
"""
source_strengths(panels :: AbstractVector{<: Panel2D}, u) = let u_ref = Ref(u); @. dot(u_ref, panel_normal(panels)) end

"""
    boundary_vector(panels, u)

Creates the vector for the boundary condition of the problem given an array of Panel2Ds and velocity ``u``.
"""
boundary_vector(panels :: AbstractVector{<: Panel2D}, u) = [ - source_matrix(panels, panels) * source_strengths(panels, u); 0 ]

"""
    solve_strengths(panels, u, bound)

Solve the system of equations ``[AIC][\\phi] = [\\hat{U} \\cdot \\vec{n}]`` condition given the array of Panel2Ds, a velocity ``u``, and an optional bound for the length of the wake.
"""
function solve_strengths(panels :: AbstractVector{<: Panel2D}, u, bound = 1e3) 
    AIC  = influence_matrix(panels, wake_panel(panels, bound))
    boco = boundary_vector(panels, u)

    AIC \ boco 
end

"""
    lift_coefficient(panels, φs, u)

Computes the lift coefficient given the array of Panel2Ds, the associated doublet strengths φs, and the velocity ``u``.
"""
function lift_coefficient(panels :: AbstractVector{<: Panel2D}, φs, u)
    speed = norm(u)
    u_ref = Ref(u)
    
    Δφs   = -midpair_map(-, φs[1:end-1])
    Δrs   = midpair_map(panel_dist, panels)
    vels  = @. panel_velocity(Δφs, Δrs, u_ref, panel_tangent(panels))
    cps   = @. pressure_coefficient(speed, vels)
    cls   = @. lift_coefficient(cps, Δrs / 2, panel_angle(panels))
    
    cl    = sum(cls)
end