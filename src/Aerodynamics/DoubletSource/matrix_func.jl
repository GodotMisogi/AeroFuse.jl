"""
doublet_matrix(panels_1, panels_2)

Create the matrix of doublet potential influence coefficients between pairs of panels_1 and panels_2.
"""
doublet_matrix(panels_1, panels_2) = [ doublet_influence(panel_j, panel_i) for panel_i in panels_1, panel_j in panels_2 ]

"""
    kutta_condition(panels)

Create the vector describing Morino's Kutta condition given Panel2Ds.
"""
kutta_condition(panels) = [ 1 zeros(length(panels) - 2)' -1 ]

"""
    wake_vector(wake_panel, panels)

Create the vector of doublet potential influence coefficients from the wake on the panels given the wake panel and the array of Panel2Ds.
"""
wake_vector(woke_panel :: AbstractPanel2D, panels) = doublet_influence.(Ref(woke_panel), panels)

"""
    influence_matrix(panels, wake_panel)

Assemble the Aerodynamic Influence Coefficient matrix consisting of the doublet matrix, wake vector, Kutta condition given Panel2Ds and the wake panel.
"""
influence_matrix(panels, woke_panel :: AbstractPanel2D) =
    [ doublet_matrix(panels, panels)  wake_vector(woke_panel, panels) ;
        kutta_condition(panels)                      1.               ]

"""
    source_matrix(panels_1, panels_2)

Create the matrix of source potential influence coefficients between pairs of `panels_1` and `panels_2`.
"""
source_matrix(panels_1, panels_2) = [ source_influence(panel_j, panel_i) for (panel_i, panel_j) in product(panels_1, panels_2) ]

"""
    source_strengths(panels, freestream)

Create the vector of source strengths for the Dirichlet boundary condition ``\\sigma = \\vec U_{\\infty} \\cdot \\hat{n}`` given Panel2Ds and a Uniform2D.
"""
source_strengths(panels, u) = dot.(Ref(u), panel_normal.(panels))

"""
    boundary_vector(panels, u)

Create the vector for the boundary condition of the problem given an array of Panel2Ds and velocity ``u``.
"""
boundary_vector(panels, u) = [ - source_matrix(panels, panels) * source_strengths(panels, u); 0 ]

boundary_vector(colpoints, u, r_te) = [ dot.(colpoints, Ref(u)); dot(u, r_te) ]

# boundary_vector(panels, u, r_te) = [ dot.(collocation_point.(panels), Ref(u)); dot(u, r_te) ]

function boundary_vector(panels :: Vector{<: AbstractPanel2D}, wakes :: Vector{<: AbstractPanel2D}, u) 
    source_panels = [ panels; wakes ]
    [ - source_matrix(panels, source_panels) * source_strengths(source_panels, u); 0 ]
end

"""
    solve_strengths(panels, u, sources, bound)

Solve the system of equations ``[AIC][\\phi] = [\\vec{U} \\cdot \\hat{n}] - B[\\sigma]`` condition given the array of Panel2Ds, a velocity ``\\vec U``, a condition whether to disable source terms (``\\sigma = 0``), and an optional named bound for the length of the wake.
"""
function solve_strengths(panels, u, α, r_te, sources :: Bool; bound = 1e2)
    # Wake
    woke_panel  = wake_panel(panels, bound, α)
    woke_vector = wake_vector(woke_panel, panels)
    woke_matrix = [ -woke_vector zeros(length(panels), length(panels) -2) woke_vector ]
    
    # AIC
    AIC 	= doublet_matrix(panels, panels) + woke_matrix
    boco 	= dot.(collocation_point.(panels), Ref(u)) - woke_vector * dot(u, r_te)
    
    # AIC
    # AIC 	= influence_matrix(panels, woke_panel)
    # boco 	= boundary_vector(ifelse(sources, panels, collocation_point.(panels)), u, r_te) - [ woke_vector; 0 ] .* dot(u, r_te) 

    AIC \ boco 
end

"""
    tangential_velocities(panels, φs, u, sources :: Bool)

Compute the tangential velocities and panel distances given the array of `Panel2D`s, their associated doublet strengths ``\\phi``s, the velocity ``u``, and a condition whether to disable source terms (``\\sigma = 0``).
"""
function tangential_velocities(panels, φs, u, sources :: Bool)
    # Δrs   = midpair_map(panel_dist, panels)
    # Δφs   = -midpair_map(-, φs[1:end-1])

    Δrs   	 = @. panel_dist(panels[2:end], panels[1:end-1])
    Δφs   	 = @. φs[1:end-1] - φs[2:end] 
    tangents = @. panel_tangent(panels[2:end])

    vels  = ifelse(sources, panel_velocity.(Δφs, Δrs, Ref(u), tangents), Δφs ./ Δrs)

    vels, Δrs
end

## WAKE VERSIONS
#==========================#

"""
    solve_strengths(panels, u, wakes, bound)

Solve the system of equations ``[AIC][\\phi] = [\\vec{U} \\cdot \\hat{n}] - B[\\sigma]`` condition given the array of Panel2Ds, a velocity ``\\vec U``, a vector of wake `Panel2D`s, and an optional named bound for the length of the wake.
"""
function solve_strengths(panels, u, α, wakes; bound = 1e2) 
    AIC  = influence_matrix(panels, wake_panel(panels, bound, α))
    boco = boundary_vector(panels, wakes, u)

    AIC \ boco 
end