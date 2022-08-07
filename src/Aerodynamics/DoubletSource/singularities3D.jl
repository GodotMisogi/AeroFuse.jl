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

"""
function check_panel_status(panel :: AbstractPanel3D, point :: Point3D, transformation_error = 1e-7)
	if prod([panel.p1.z; panel.p2.z; panel.p3.z; panel.p4.z] .>= transformation_error)
		@warn "Panel is not in local coordinate system! Check AIC matrix assembly!"
		@info "Transform to local coordinate system now."
		return transform_panel(panel, point)
	else
		return panel, point
	end
end

"""
    quadrilateral_source_velocity(σ, local_panel :: AbstractPanel, local_point)

    Compute the panel velocity (û,v̂,ŵ) in local panel (x̂,ŷ,ẑ) direction of a constant source.
"""
function quadrilateral_source_velocity(σ, local_panel :: AbstractPanel3D, local_point :: Point3D)
    # A panel of type `Panel3D` has 4 points of type `SVector{3, Any}` whose coordinates are in the format of
    # [
    #     [x1 y1 0]
    #     [x2 y2 0]
    #     [x3 y3 0]
    #     [x4 y4 0]
    # ]

	const_quad_source_u(i, j, yi, dij, ri) = (yi[j] -yi[i]) / dij[i] * log( (ri[i] + ri[j] - dij[i]) / (ri[i] + ri[j] + dij[i]) )
	const_quad_source_v(i, j, xi, dij, ri) = (xi[i] -xi[j]) / dij[i] * log( (ri[i] + ri[j] - dij[i]) / (ri[i] + ri[j] + dij[i]) )

    coord = panel_coordinates(local_panel)

    # Check whether the panel is expressed in local coordinates. All z-coordinates must be zeros.
    local_panel, local_point = check_panel_status(local_panel, local_point)
    z = local_point.z

	xi = xs(local_panel)
    yi = ys(local_panel)
    ri  = @SVector [r_k(i, coord, local_point)  for i=1:4]
    ei  = @SVector [e_k(i, coord, local_point)  for i=1:4]
    hi  = @SVector [h_k(i, coord, local_point)  for i=1:4]
    dij = @SVector [d_ij(i, i%4+1, coord)       for i=1:4]
    mij = @SVector [m_ij(i, i%4+1, coord)       for i=1:4]

    # Results from textbook differs from a minus sign
    u = σ / 4π * sum(
        [const_quad_source_u(i, i%4+1, yi, dij, ri) for i=1:4]
    )

    v = σ / 4π * sum(
        [const_quad_source_v(i, i%4+1, xi, dij, ri) for i=1:4]
    )

    w = σ / 4π * sum(
        [const_quad_source_phi_term2(i, i%4+1, z, mij, ei, hi, ri) for i=1:4]
    )

    return SVector(u, v, w)
end


"""
    quadrilateral_source_velocity_farfield(σ, local_panel :: AbstractPanel, local_point)

Compute the panel velocity (û,v̂,ŵ) in local panel (x̂,ŷ,ẑ) direction of a constant source at FARFIELD.
"""
function quadrilateral_source_velocity_farfield(σ, local_panel :: AbstractPanel3D, local_point :: Point3D)
    check_panel_status(local_panel, local_point)

    A = panel_area(local_panel)
    centroid = midpoint(local_panel) # It is better to use centroid

    u = ( σ * A * (local_point.x - centroid.x) ) / ( 4π * norm(local_point - centroid)^3 )
    v = ( σ * A * (local_point.y - centroid.y) ) / ( 4π * norm(local_point - centroid)^3 )
    w = ( σ * A * (local_point.z - centroid.z) ) / ( 4π * norm(local_point - centroid)^3 )

    return SVector(u, v, w)
end

const_quad_source_phi_term1(i::Int64, j::Int64, x, y, xi, yi, dij, ri) = ( (x - xi[i]) * (yi[j] - yi[i]) - (y - yi[i]) * (xi[j] - xi[i]) ) / dij[i] * log( (ri[i] + ri[j] + dij[i]) / (ri[i] + ri[j] - dij[i]) )
const_quad_source_phi_term2(i::Int64, j::Int64, z, mij, ei, hi, ri) = atan( (mij[i] * ei[i] - hi[i]) / (z * ri[i]) ) - atan( (mij[i] * ei[j] - hi[j]) / (z * ri[j]) )

"""
    quadrilateral_source_potential(σ, local_panel :: AbstractPanel, local_point)

Compute the flow potential in local panel (x̂,ŷ,ẑ) direction of a constant source.
"""
function quadrilateral_source_potential(σ, local_panel :: AbstractPanel3D, local_point:: Point3D)
    coord = panel_coordinates(local_panel)

    # Check whether the panel is expressed in local coordinates. All z-coordinates must be zeros.
    check_panel_status(local_panel, local_point)

    xi = xs(local_panel)
    yi = ys(local_panel)

    x, y, z = local_point.x, local_point.y, local_point.z

    ri  = @SVector [r_k(i, coord, local_point)  for i=1:4]
    ei  = @SVector [e_k(i, coord, local_point)  for i=1:4]
    hi  = @SVector [h_k(i, coord, local_point)  for i=1:4]
    dij = @SVector [d_ij(i, i%4+1, coord)       for i=1:4]
    mij = @SVector [m_ij(i, i%4+1, coord)       for i=1:4]

    return -σ / 4π * (
        sum(
            [const_quad_source_phi_term1(i, i%4+1, x, y, xi, yi, dij, ri) for i=1:4]
        ) -
        abs(z) * sum(
            [const_quad_source_phi_term2(i, i%4+1, z, mij, ei, hi, ri) for i=1:4]
        )
    )
end


# """
#     const_quad_doublet_den(i, j, ri, xi, yi, x, y, z)

# Helper function: denominator for quadrilateral_doublet().
# """
# const_quad_doublet_den(i, j, ri, xi, yi, x, y, z) = ri[i] * ri[j] * ( ri[i] * ri[j] - ( (x - xi[i]) * (x - xi[j]) + (y - yi[i]) * (y - yi[j]) + z^2 ) )


# """
#     const_quad_doublet_num_u(i, j, ri, yi, z)

# Helper function: numerator of û for quadrilateral_doublet().
# """
# const_quad_doublet_num(i, j, ri, yi, z) = z * (yi[i] - yi[j]) * (ri[i] + ri[j])


# """
#     const_quad_doublet_num_u(i, j, ri, yi, z)

# Helper function: numerator of ŵ for quadrilateral_doublet().
# """
# const_quad_doublet_num_w(i, j, ri, xi, yi, x, y) = (ri[i] + ri[j]) * ((x - xi[j]) * (y - yi[i]) - (x - xi[i]) * (y - yi[j]))


"""
    quadrilateral_doublet_velocity(σ, local_panel :: AbstractPanel, local_point)

Compute the panel velocity (û,v̂,ŵ) in local panel (x̂,ŷ,ẑ) direction of a constant doublet panel.
"""
function quadrilateral_doublet_velocity(μ, local_panel :: AbstractPanel3D, local_point :: Point3D)
    # Check whether the panel is expressed in local coordinates. All z-coordinates must be zeros.
    check_panel_status(local_panel, local_point)

    coord = panel_coordinates(local_panel)

    rvs = @SVector [local_point - coord[i] for i=1:4]
    rs = norm.(rvs)
    ds = @SVector [d_ij(i, i%4+1, coord) for i=1:4]
    
    return μ / 4π * sum(
        [ ( rvs[i] × rvs[i%4+1] * (rs[i] + rs[i%4+1]) ) / ( (rs[i] * rs[i%4+1]) * (rs[i] * rs[i%4+1] + rvs[i] ⋅ rvs[i%4+1]) + 0.005 * ds[i] ) for i=1:4 ]
    )

end

"""
    quadrilateral_doublet_velocity_farfield(σ, local_panel :: AbstractPanel, local_point)

Compute the panel velocity (û,v̂,ŵ) in local panel (x̂,ŷ,ẑ) direction of a constant doublet panel at FARFIELD.
"""
function quadrilateral_doublet_velocity_farfield(μ, local_panel :: AbstractPanel3D, local_point:: Point3D)
    check_panel_status(local_panel, local_point)

    A = panel_area(local_panel)
    centroid = midpoint(local_panel) # It is better to use centroid
    x, y, z = local_point.x, local_point.y, local_point.z
    x0, y0 = centroid.x, centroid.y
    u = ( 3μ * A * (x - x0) * z ) / ( 4π * norm(local_point - centroid)^5 )
    v = ( 3μ * A * (y - y0) * z ) / ( 4π * norm(local_point - centroid)^5 )
    w = ( -μ * A * ( (x - x0)^2 + (y - y0)^2 - 2 * z^2 ) ) / ( 4π * norm(local_point - centroid)^5 )

    return SVector(u, v, w)
end


"""
    quadrilateral_doublet_potential(σ, local_panel :: AbstractPanel, local_point)

Compute the flow potential in local panel (x̂,ŷ,ẑ) direction of a constant doublet panel.
"""
function quadrilateral_doublet_potential(μ, local_panel :: AbstractPanel3D, local_point:: Point3D)
    coord = panel_coordinates(local_panel)

    # Check whether the panel is expressed in local coordinates. All z-coordinates must be zeros.
    check_panel_status(local_panel, local_point)

    z = local_point.z

    ri  = @SVector [r_k(i, coord, local_point)  for i=1:4]
    ei  = @SVector [e_k(i, coord, local_point)  for i=1:4]
    hi  = @SVector [h_k(i, coord, local_point)  for i=1:4]
    mij = @SVector [m_ij(i, i%4+1, coord)       for i=1:4]
    
    return μ / 4π * sum(
        [const_quad_source_phi_term2(i, i%4+1, z, mij, ei, hi, ri) for i=1:4]
    )
end