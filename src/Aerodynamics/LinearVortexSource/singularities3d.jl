"""
    d_ij(i :: Int64, j :: Int64, local_panel :: AbstractPanel)

Compute the distance from i-th to j-th point of the panel.
"""
function d_ij(i :: Int64, j :: Int64, local_panel :: AbstractPanel)
    local_coordinates = panel_coordinates(local_panel)
    return norm(local_coordinates[j] - local_coordinates[i])
end

"""
    m_ij(i :: Int64, j :: Int64, local_panel :: AbstractPanel)

Compute the slope of the line segment consisting i-th and j-th point of the panel.
"""
function m_ij(i :: Int64, j :: Int64, local_panel :: AbstractPanel)
    local_coordinates = panel_coordinates(local_panel)
    x,y = local_coordinates[j] - local_coordinates[i]
    if abs(x) < 1e-12
        return Inf
    else
        return y/x
    end
end

"""
    r_k(k :: Int64, local_panel :: AbstractPanel, local_point :: Point3D)

Compute the distance from i-th point of the panel to the arbitary point.
"""
function r_k(k :: Int64, local_panel :: AbstractPanel, local_point :: Point3D)
    local_coordinates = panel_coordinates(local_panel)
    return norm(local_point - local_coordinates[k])
end

"""
    e_k(k :: Int64, local_panel :: AbstractPanel, local_point :: Point3D)

Compute (x - xk)^2 + z^2.
"""
function e_k(k :: Int64, local_panel :: AbstractPanel, local_point :: Point3D)
    local_coordinates = panel_coordinates(local_panel)
    return norm([local_point.x - local_coordinates[k].x, local_point.z])
end

"""
    h_k(k :: Int64, local_panel :: AbstractPanel, local_point :: Point3D)

Compute (x - xk) * (y - yk).
"""
function h_k(k :: Int64, local_panel :: AbstractPanel, local_point :: Point3D)
    local_coordinates = panel_coordinates(local_panel)
    return (local_point.x - local_coordinates[k].x) * (local_point.y - local_coordinates[k].y)
end

# function constant_quadrilateral_source_velocity()

# end