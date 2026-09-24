const_quad_source_phi_term1(x, y, dij, xi, yi, xj, yj, ri, rj) = ((x - xi) * (yj - yi) - (y - yi) * (xj - xi)) / dij * log((ri + rj + dij) / (ri + rj - dij))

const_quad_source_phi_term2(z, mij, ei, hi, ri, ej, hj, rj) = atan((mij * ei - hi) / (z * ri)) - atan((mij * ej - hj) / (z * rj))

function potential_processing(p, p_i, p_j)
    # Dereferencing
    x, y, z   = p
    xi, yi, _ = p_i
    xj, yj, _ = p_j

    # Distances
    rj  = norm(p - p_j)
    ej  = (x - xj)^2 + z^2 
    hj  = (x - xj) * (y - yj)
    dij = √((xj - xi)^2 + (yj - yi)^2)

    # Local gradient
    mij = (yj - yi) / (xj - xi)

    mij, dij, ej, hj, rj
end

# Axis permutation into the kernel's local frame. Static and allocation-free (with the corner
# tuple below) so the kernel also compiles for GPU backends.
const DOUBLET_AXIS_PERMUTATION = SMatrix{3,3}(0, 1, 0, 1, 0, 0, 0, 0, -1)

function quadrilateral_doublet_potential(μ, panel :: AbstractPanel3D, point)
    # Local coordinate system transformation
    T = get_transformation(panel, DOUBLET_AXIS_PERMUTATION)
    panel, p = T(panel), T(point)
    pans = (panel.p1, panel.p2, panel.p3, panel.p4)

    x, y, z = p

    ε = eltype(p)(1e-10)
    if abs(z) <= ε
        return (abs(x) <= ε && abs(y) <= ε) ? μ / 2 : zero(μ)
    else
        # LESS READABLE, BUT CORRECT

        # Initial condition for first point
        xi, yi, _ = p_i = pans[1]
        ri  = norm(p - p_i)
        ei  = (x - xi)^2 + z^2
        hi  = (x - xi) * (y - yi)

        inf = zero(μ)
        for i = 1:4
            # Next point (could be more efficient)
            j = i % 4 + 1

            mij, _, ej, hj, rj = potential_processing(p, pans[i], pans[j])

            # Evaluate influence term
            inf += const_quad_source_phi_term2(z, mij, ei, hi, ri, ej, hj, rj)

            # Now set to the new point's values
            ri = rj
            ei = ej
            hi = hj
        end

        return -μ / 4 / π * inf
    end
end

function quadrilateral_source_potential(σ, panel :: AbstractPanel3D, p)
    pans = panel_coordinates(panel)

    x, y, z = p
    xi, yi, _ = p_i = pans[1]
    ri = norm(p - p_i)
    ei = (x - xi)^2 + z^2
    hi = (x - xi) * (y - yi)

    inf = zero(σ)
    for i = 1:4
        # Next point (could be more efficient)
        j = i % 4 + 1

        @views xi, yi, _ = pans[i]
        @views xj, yj, _ = pans[j]

        mij, dij, ej, hj, rj = potential_processing(p, pans[i], pans[j])

        # Evaluate influence term
        inf += const_quad_source_phi_term1(x, y, dij, xi, yi, xj, yj, ri, rj) - abs(z) * const_quad_source_phi_term2(z, mij, ei, hi, ri, ej, hj, rj)

        # Now set to the new point's values
        ri = rj
        ei = ej
        hi = hj
    end

    return -σ / 4π * inf
end