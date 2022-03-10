"""
    d_ij(i :: Int64, j :: Int64, local_coordinates)

Compute the distance from i-th to j-th point of the panel.
"""
d_ij(i :: Int64, j :: Int64, local_coordinates) = norm(local_coordinates[j] - local_coordinates[i])


"""
    m_ij(i :: Int64, j :: Int64, local_coordinates)

Compute the slope of the line segment consisting i-th and j-th point of the panel.
"""
function m_ij(i :: Int64, j :: Int64, local_coordinates)
    x, y = local_coordinates[j] - local_coordinates[i]
    return y/x
end


"""
    r_k(k :: Int64, local_coordinates, local_point)

Compute the distance from i-th point of the panel to the arbitary point.
"""
r_k(k :: Int64, local_coordinates, local_point) = norm(local_point - local_coordinates[k])


"""
    r_k(k :: Int64, local_panel :: AbstractPanel3D, local_point)

Compute the distance from i-th point of the panel to the arbitary point.
"""
r_k(k :: Int64, local_panel :: AbstractPanel3D, local_point) = r_k(k, panel_coordinates(local_panel), local_point)


"""
    k :: Int64, local_coordinates, local_point

Compute (x - xk)^2 + z^2.
"""
e_k(k :: Int64, local_coordinates, local_point) = norm([local_point.x - local_coordinates[k].x, local_point.z])


"""
    k :: Int64, local_coordinates, local_point

Compute (x - xk) * (y - yk).
"""
h_k(k :: Int64, local_coordinates, local_point) = (local_point.x - local_coordinates[k].x) * (local_point.y - local_coordinates[k].y)


"""
    const_quad_source_u(i, j, x, y, xi, yi, dij, ri)

Helper function: term u for constant_quadrilateral_source_velocity().
"""
const_quad_source_u(i, j, yi, dij, ri) = (yi[j] -yi[i]) / dij[i] * log( (ri[i] + ri[j] - dij[i]) / (ri[i] + ri[j] + dij[i]) )


"""
    const_quad_source_v(i, j, x, y, xi, yi, dij, ri)

Helper function: term v for constant_quadrilateral_source_velocity().
"""
const_quad_source_v(i, j, xi, dij, ri) = (xi[i] -xi[j]) / dij[i] * log( (ri[i] + ri[j] - dij[i]) / (ri[i] + ri[j] + dij[i]) )


"""
    constant_quadrilateral_source_velocity(σ, local_panel :: AbstractPanel, local_point)

    Compute the panel velocity (û,v̂,ŵ) in local panel (x̂,ŷ,ẑ) direction of a constant source.
"""
function constant_quadrilateral_source_velocity(σ, local_panel :: AbstractPanel3D, local_point :: Point3D)
    # A panel of type `Panel3D` has 4 points of type `SVector{3, Any}` whose coordinates are in the format of
    # [
    #     [x1 y1 0]
    #     [x2 y2 0]
    #     [x3 y3 0]
    #     [x4 y4 0]
    # ]
    coord = panel_coordinates(local_panel)

    # Check whether the panel is expressed in local coordinates. All z-coordinates must be zeros.
    local_zs = zs(local_panel)
    if norm(local_zs) >= 1e-10
        AssertionError("Panel is not in its local coordinate!")
    end

    xi = xs(local_panel)
    yi = ys(local_panel)

    x, y, z = local_point.x, local_point.y, local_point.z

    ri = [
        r_k(1, coord, local_point),
        r_k(2, coord, local_point),
        r_k(3, coord, local_point),
        r_k(4, coord, local_point)
    ]

    dij = [
        d_ij(1, 2, coord)
        d_ij(2, 3, coord)
        d_ij(3, 4, coord)
        d_ij(4, 1, coord)
    ]

    mij = [
        m_ij(1, 2, coord)
        m_ij(2, 3, coord)
        m_ij(3, 4, coord)
        m_ij(4, 1, coord)
    ]

    ei = [
        e_k(1, coord, local_point)
        e_k(2, coord, local_point)
        e_k(3, coord, local_point)
        e_k(4, coord, local_point)
    ]
    
    hi = [
        h_k(1, coord, local_point)
        h_k(2, coord, local_point)
        h_k(3, coord, local_point)
        h_k(4, coord, local_point)
    ]

    # Results from textbook differs from a minus sign
    u = -σ / 4π * (
        const_quad_source_u(1, 2, yi, dij, ri) +
        const_quad_source_u(2, 3, yi, dij, ri) +
        const_quad_source_u(3, 4, yi, dij, ri) +
        const_quad_source_u(4, 1, yi, dij, ri)
    )

    v = -σ / 4π * (
        const_quad_source_v(1, 2, yi, dij, ri) +
        const_quad_source_v(2, 3, yi, dij, ri) +
        const_quad_source_v(3, 4, yi, dij, ri) +
        const_quad_source_v(4, 1, yi, dij, ri)
    )

    w = -σ / 4π * (
        const_quad_source_phi_term2(1, 2, z, mij, ei, hi, ri) + 
        const_quad_source_phi_term2(2, 3, z, mij, ei, hi, ri) +
        const_quad_source_phi_term2(3, 4, z, mij, ei, hi, ri) +
        const_quad_source_phi_term2(4, 1, z, mij, ei, hi, ri)
    )

    @. return SVector(u, v, w)
end


"""
    quadrilateral_panel_area(panel :: AbstractPanel3D)

Compute the area of a quadrilateral panel.
"""
function quadrilateral_panel_area(panel :: AbstractPanel3D)
    coord = panel_coordinates(panel)

    d12 = d_ij(1, 2, coord)
    d23 = d_ij(2, 3, coord)
    d34 = d_ij(3, 4, coord)
    d41 = d_ij(4, 1, coord)
    d13 = d_ij(1, 3, coord)

    s1 = 0.5 * (d12 + d23 + d13)
    s2 = 0.5 * (d34 + d41 + d13)

    A1 = sqrt(s1 * (s1 - d12) * (s1 - d23) * (s1 - d13))
    A2 = sqrt(s2 * (s2 - d34) * (s2 - d41) * (s2 - d13))

    return A1 + A2
end


"""
    constant_quadrilateral_source_velocity_farfield(σ, local_panel :: AbstractPanel, local_point)

Compute the panel velocity (û,v̂,ŵ) in local panel (x̂,ŷ,ẑ) direction of a constant source at FARFIELD.
"""
function constant_quadrilateral_source_velocity_farfield(σ, local_panel :: AbstractPanel3D, local_point :: Point3D)
    local_zs = zs(local_panel)
    if norm(local_zs) >= 1e-10
        AssertionError("Panel is not in its local coordinate!")
    end

    A = quadrilateral_panel_area(local_panel)
    centroid = midpoint(local_panel) # It is better to use centroid

    u = ( σ * A * (local_point.x - centroid.x) ) / ( 4π * norm(local_point - centroid)^3 )
    v = ( σ * A * (local_point.y - centroid.y) ) / ( 4π * norm(local_point - centroid)^3 )
    w = ( σ * A * (local_point.z - centroid.z) ) / ( 4π * norm(local_point - centroid)^3 )

    @. return SVector(u, v, w)
end


"""
    const_quad_source_phi_term1(i, j, x, y, xi, yi, dij, ri)

Helper function: term 1 for constant_quadrilateral_source_potential().
"""
const_quad_source_phi_term1(i, j, x, y, xi, yi, dij, ri) = ( (x - xi[i]) * (yi[j] - yi[i]) - (y - yi[i])(xi[j] - xi[i]) ) / dij[i] * log( (ri[i] + ri[j] + dij[i]) / (ri[i] + ri[j] - dij[i]) )

"""
    const_quad_source_phi_term2(i, j, z, mij, ei, hi, ri)
Helper function: term 2 for constant_quadrilateral_source_potential().
"""
const_quad_source_phi_term2(i, j, z, mij, ei, hi, ri) = atan( (mij[i] * ei[i] - hi[i]) / (z * ri[i]) ) - atan( (mij[i] * ei[j] - hi[j]) / (z * ri[j]) )


"""
    constant_quadrilateral_source_potential(σ, local_panel :: AbstractPanel, local_point)

Compute the flow potential in local panel (x̂,ŷ,ẑ) direction of a constant source.
"""
function constant_quadrilateral_source_potential(μ, local_panel :: AbstractPanel3D, local_point:: Point3D)
    coord = panel_coordinates(local_panel)

    # Check whether the panel is expressed in local coordinates. All z-coordinates must be zeros.
    local_zs = zs(local_panel)
    if norm(local_zs) >= 1e-10
        AssertionError("Panel is not in its local coordinate!")
    end

    xi = xs(local_panel)
    yi = ys(local_panel)

    x, y, z = local_point.x, local_point.y, local_point.z

    ri = [
        r_k(1, coord, local_point),
        r_k(2, coord, local_point),
        r_k(3, coord, local_point),
        r_k(4, coord, local_point)
    ]

    dij = [
        d_ij(1, 2, coord)
        d_ij(2, 3, coord)
        d_ij(3, 4, coord)
        d_ij(4, 1, coord)
    ]

    mij = [
        m_ij(1, 2, coord)
        m_ij(2, 3, coord)
        m_ij(3, 4, coord)
        m_ij(4, 1, coord)
    ]

    ei = [
        e_k(1, coord, local_point)
        e_k(2, coord, local_point)
        e_k(3, coord, local_point)
        e_k(4, coord, local_point)
    ]

    
    hi = [
        h_k(1, coord, local_point)
        h_k(2, coord, local_point)
        h_k(3, coord, local_point)
        h_k(4, coord, local_point)
    ]

    ϕ = -σ / 4π * (
        const_quad_source_phi_term1(1, 2, x, y, xi, yi, dij, ri) + 
        const_quad_source_phi_term1(2, 3, x, y, xi, yi, dij, ri) +
        const_quad_source_phi_term1(3, 4, x, y, xi, yi, dij, ri) +
        const_quad_source_phi_term1(4, 1, x, y, xi, yi, dij, ri) -
        abs(z) * (
            const_quad_source_phi_term2(1, 2, z, mij, ei, hi, ri) + 
            const_quad_source_phi_term2(2, 3, z, mij, ei, hi, ri) +
            const_quad_source_phi_term2(3, 4, z, mij, ei, hi, ri) +
            const_quad_source_phi_term2(4, 1, z, mij, ei, hi, ri)
        )
    )

    return ϕ
end


"""
    const_quad_doublet_den(i, j, ri, xi, yi, x, y, z)

Helper function: denominator for constant_quadrilateral_doublet().
"""
const_quad_doublet_den(i, j, ri, xi, yi, x, y, z) = ri[i] * ri[j] * ( ri[i] * ri[j] - ( (x - xi[i]) * (x - xi[j]) + (y - yi[i]) * (y - yi[j]) + z^2 ) )


"""
    const_quad_doublet_num_u(i, j, ri, yi, z)

Helper function: numerator of û for constant_quadrilateral_doublet().
"""
const_quad_doublet_num_u(i, j, ri, yi, z) = z * (yi[i] - yi[j]) * (ri[i] + ri[j])


"""
    const_quad_doublet_num_u(i, j, ri, yi, z)

Helper function: numerator of v̂ for constant_quadrilateral_doublet().
"""
const_quad_doublet_num_v(i, j, ri, xi, z) = z * (xi[j] - xi[i]) * (ri[i] + ri[j])


"""
    const_quad_doublet_num_u(i, j, ri, yi, z)

Helper function: numerator of ŵ for constant_quadrilateral_doublet().
"""
const_quad_doublet_num_w(i, j, ri, xi, yi, x, y) = (ri[i] + ri[j]) * ((x - xi[j]) * (y - yi[i]) - (x - xi[i]) * (y - yi[j]))


"""
    constant_quadrilateral_doublet_velocity(σ, local_panel :: AbstractPanel, local_point)

Compute the panel velocity (û,v̂,ŵ) in local panel (x̂,ŷ,ẑ) direction of a constant doublet panel.
"""
function constant_quadrilateral_doublet_velocity(μ, local_panel :: AbstractPanel3D, local_point :: Point3D)
    local_zs = zs(local_panel)
    if norm(local_zs) >= 1e-10
        AssertionError("Panel is not in its local coordinate!")
    end

    ri = SVector(
        r_k(1, local_panel, local_point), 
        r_k(2, local_panel, local_point), 
        r_k(3, local_panel, local_point),
        r_k(4, local_panel, local_point)
        )
    xi = xs(local_panel)
    yi = ys(local_panel)
    x, y, z = local_point.x, local_point.y, local_point.z

    u = μ / 4π * (
        const_quad_doublet_num_u(1, 2, ri, yi, z) / const_quad_doublet_den(1, 2, ri, xi, yi, x, y, z) +
        const_quad_doublet_num_u(2, 3, ri, yi, z) / const_quad_doublet_den(2, 3, ri, xi, yi, x, y, z) +
        const_quad_doublet_num_u(3, 4, ri, yi, z) / const_quad_doublet_den(3, 4, ri, xi, yi, x, y, z) +
        const_quad_doublet_num_u(4, 1, ri, yi, z) / const_quad_doublet_den(4, 1, ri, xi, yi, x, y, z)
    )

    v = μ / 4π * (
        const_quad_doublet_num_v(1, 2, ri, xi, z) / const_quad_doublet_den(1, 2, ri, xi, yi, x, y, z) +
        const_quad_doublet_num_v(2, 3, ri, xi, z) / const_quad_doublet_den(2, 3, ri, xi, yi, x, y, z) +
        const_quad_doublet_num_v(3, 4, ri, xi, z) / const_quad_doublet_den(3, 4, ri, xi, yi, x, y, z) +
        const_quad_doublet_num_v(4, 1, ri, xi, z) / const_quad_doublet_den(4, 1, ri, xi, yi, x, y, z)
    )

    w = μ / 4π * (
        const_quad_doublet_num_w(1, 2, ri, xi, yi, x, y) / const_quad_doublet_den(1, 2, ri, xi, yi, x, y, z) +
        const_quad_doublet_num_w(2, 3, ri, xi, yi, x, y) / const_quad_doublet_den(2, 3, ri, xi, yi, x, y, z) +
        const_quad_doublet_num_w(3, 4, ri, xi, yi, x, y) / const_quad_doublet_den(3, 4, ri, xi, yi, x, y, z) +
        const_quad_doublet_num_w(4, 1, ri, xi, yi, x, y) / const_quad_doublet_den(4, 1, ri, xi, yi, x, y, z)
    )

    @. return SVector(u, v, w)
end

"""
    constant_quadrilateral_doublet_velocity_farfield(σ, local_panel :: AbstractPanel, local_point)

Compute the panel velocity (û,v̂,ŵ) in local panel (x̂,ŷ,ẑ) direction of a constant doublet panel at FARFIELD.
"""
function constant_quadrilateral_doublet_velocity_farfield(μ, local_panel :: AbstractPanel3D, local_point:: Point3D)
    local_zs = zs(local_panel)
    if norm(local_zs) >= 1e-10
        AssertionError("Panel is not in its local coordinate!")
    end

    A = quadrilateral_panel_area(local_panel)
    centroid = midpoint(local_panel) # It is better to use centroid
    x, y, z = local_point.x, local_point.y, local_point.z
    x0, y0 = centroid.x, centroid.y
    u = ( 3μ * A * (x - x0) * z ) / ( 4π * norm(local_point - centroid)^5 )
    v = ( 3μ * A * (y - y0) * z ) / ( 4π * norm(local_point - centroid)^5 )
    w = ( -μ * A * ( (x - x0)^2 + (y - y0)^2 - 2 * z^2 ) ) / ( 4π * norm(local_point - centroid)^5 )

    @. return SVector(u, v, w)
end


"""
    constant_quadrilateral_doublet_potential(σ, local_panel :: AbstractPanel, local_point)

Compute the flow potential in local panel (x̂,ŷ,ẑ) direction of a constant doublet panel.
"""
function constant_quadrilateral_doublet_potential(μ, local_panel :: AbstractPanel3D, local_point:: Point3D)
    coord = panel_coordinates(local_panel)

    # Check whether the panel is expressed in local coordinates. All z-coordinates must be zeros.
    local_zs = zs(local_panel)
    if norm(local_zs) >= 1e-10
        AssertionError("Panel is not in its local coordinate!")
    end

    z = local_point.z

    ri = [
        r_k(1, coord, local_point),
        r_k(2, coord, local_point),
        r_k(3, coord, local_point),
        r_k(4, coord, local_point)
    ]

    mij = [
        m_ij(1, 2, coord)
        m_ij(2, 3, coord)
        m_ij(3, 4, coord)
        m_ij(4, 1, coord)
    ]

    ei = [
        e_k(1, coord, local_point)
        e_k(2, coord, local_point)
        e_k(3, coord, local_point)
        e_k(4, coord, local_point)
    ]

    
    hi = [
        h_k(1, coord, local_point)
        h_k(2, coord, local_point)
        h_k(3, coord, local_point)
        h_k(4, coord, local_point)
    ]

    ϕ = μ / 4π * (
        const_quad_source_phi_term2(1, 2, z, mij, ei, hi, ri) +
        const_quad_source_phi_term2(2, 3, z, mij, ei, hi, ri) +
        const_quad_source_phi_term2(3, 4, z, mij, ei, hi, ri) +
        const_quad_source_phi_term2(4, 1, z, mij, ei, hi, ri)
    )

    return ϕ
end