
"""
doublet_matrix(panels_1, panels_2)

Creates the matrix of doublet potential influence coefficients between pairs of panels_1 and panels_2.
"""
doublet_matrix(panels_1 :: AbstractVector{<: Panel2D}, panels_2 :: AbstractVector{<: Panel2D}) = [ panel_i === panel_j ? 0.5 : doublet_influence(panel_j, panel_i) for panel_i ∈ panels_1, panel_j ∈ panels_2 ]

"""
source_matrix(panels_1, panels_2)

Creates the matrix of source potential influence coefficients between pairs of panels_1 and panels_2.
"""
source_matrix(panels_1 :: AbstractVector{<: Panel2D}, panels_2 :: AbstractVector{<: Panel2D}) = [ source_influence(panel_j, panel_i) for panel_i ∈ panels_1, panel_j ∈ panels_2 ]

"""
    source_strengths(panels, freestream)

Creates the vector of source strengths for the Dirichlet boundary condition ``\\sigma = \\vec U_{\\infty} \\cdot \\hat{n}`` given Panel2Ds and a Uniform2D.
"""
source_strengths(panels :: AbstractVector{<: Panel2D}, u) = dot.(Ref(u), panel_normal.(panels))

"""
kutta_condition(panels)

Creates the vector describing Morino's Kutta condition given Panel2Ds.
"""
kutta_condition(panels :: AbstractVector{<: Panel2D}) = [ 1; -1; zeros(length(panels) - 4); 1; -1]

"""
wake_vector(panels, bound)

Creates the vector of doublet potential influence coefficients from the wake on the panels, with an option for the length of the wake.
"""
function wake_vector(panels :: AbstractVector{<: Panel2D}, bound = 1e3)
lastx, lasty = point2(panels[end])
woke_panel = Panel2D{eltype(lastx)}((lastx, lasty), (bound * lastx, lasty))

[ doublet_potential(woke_panel, 1., pt...) for pt ∈ collocation_point.(panels) ]
end

"""
influence_matrix(panels)

Assembles the Aerodynamic Influence Coefficient matrix consisting of the doublet matrix, wake vector, Kutta condiition given Panel2Ds.
"""
function influence_matrix(panels :: AbstractVector{<: Panel2D}) 
@timeit "Doublet Matrix" dub_mat = doublet_matrix(panels, panels)
@timeit "Wake Vector" wake_vec   = wake_vector(panels)
@timeit "Kutta Condition" kutta  = kutta_condition(panels)'
@timeit "Matrix Assembly" mat    = [ dub_mat   wake_vec;
                                     kutta        0.   ]
end


"""
boundary_vector(panels, freestream)

Assembles the vector for the boundary condition of the problem given Panel2Ds and a Uniform2D.
"""
boundary_vector(panels :: AbstractVector{<: Panel2D}, u) = [ boundary_condition(panels, u); 0. ]