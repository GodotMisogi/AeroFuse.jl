## Visualization
#===========================================================================#

function doublet_velocity(strength, panel_j :: Panel2D, panel_i :: Panel2D)
    xp, yp = transform_panel(panel_j, panel_i)
    u, v = doublet_velocity(strength, xp, yp, 0., panel_length(panel_j))

    inverse_rotation(u, v, -panel_angle(panel_i))
end

stream_velocity(x, φs, panels, V) = V .+ sum(doublet_velocity(φ_j, p_j, p_i) for (ψ_j, p_j) in zip(φs, panels) for p_i in panels)

function streamlines(point, V, φs, panels, length, num_steps :: Integer)
    streamlines = fill(point, num_steps)
    cuck(x) = stream_velocity(x, φs, panels, V)
end