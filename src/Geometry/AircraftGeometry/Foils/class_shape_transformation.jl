## Class shape transformation method
#==========================================================================================#

# Basic shape function
function shape_function(x, basis_func, coeffs, coeff_LE = 0)
    n     = length(coeffs)
    terms = basis_func.(x, n - 1, 0:n-1)
    dot(coeffs, terms) + coeff_LE * (x^0.5) * (1 - x)^(n - 0.5)
end

# Computing coordinates
CST_coordinates(class_func, basis_func, x, alphas, dz, coeff_LE, args...) = class_func(x) * shape_function(x, basis_func, alphas, coeff_LE) + x * dz

## Bernstein basis
#==========================================================================================#

bernstein_class(x, N1, N2) = x^N1 * (1 - x)^N2
bernstein_basis(x, n, k)   = binomial(n, k) * bernstein_class(x, k, n - k)

"""
    kulfan_CST(alpha_u, alpha_l,
               (Δz_u, Δz_l) = (0., 0.),
               (LE_u, LE_l) = (0., 0.),
               n            = 40)

Define a cosine-spaced foil with ``2n`` points using the Class Shape Transformation method on a Bernstein polynomial basis for the upper and lower coordinates.

The foil is defined by arrays of coefficients ``(α_u,~ α_l)`` for the upper and lower surfaces (not necessarily of the same lengths), trailing-edge displacement values ``(Δz_u,~ Δz_l)``, and coefficients for leading edge modifications on the upper and lower surfaces at the nose.
"""
function kulfan_CST(alpha_u, alpha_l, (dz_u, dz_l) = (0., 0.), (LE_u, LE_l) = (0., 0.), n :: Integer = 40, N1 = 0.5, N2 = 1.)
    # Cosine spacing for airfoil of unit chord length
    xs = cosine_spacing(0.5, 1, n)

    # λ-function for Bernstein polynomials
    bernie(x, alphas, dz, LE) = CST_coordinates(y -> bernstein_class(y, N1, N2), bernstein_basis, x, alphas, dz, LE)

    # Upper and lower surface generation
    upper_surf = [ bernie(x, alpha_u, dz_u, LE_u) for x ∈ xs ]
    lower_surf = [ bernie(x, alpha_l, dz_l, LE_l) for x ∈ xs ]

    # Counter-clockwise ordering
    @views Foil([ xs[end:-1:2] upper_surf[end:-1:2] ;
                  xs           lower_surf           ], "Kulfan CST")
end

"""
    camber_CST(α_c, α_t,
               (Δz_u, Δz_l) :: NTuple{2, Real},
               coeff_LE = 0.,
               n :: Integer = 40)

Define a cosine-spaced foil with ``2n`` points using the Class Shape Transformation method on a Bernstein polynomial basis for the camber and thickness coordinates.

The foil is defined by arrays of coefficients ``(α_c,~ α_t)`` for the upper and lower surfaces, trailing-edge spacing values ``(Δz_u,~Δz_l)``, and a coefficient for the leading edge modifications at the nose.
"""
function camber_CST(α_cam, α_thicc, dz_thicc = 0., coeff_LE = 0, n :: Integer = 40, N1 = 0.5, N2 = 1.)
    # Cosine spacing for airfoil of unit chord length
    xs = cosine_spacing(0.5, 1, n)

    # λ-function for Bernstein polynomials
    bernie(x, αs, dz = 0.) = CST_coordinates(y -> bernstein_class(y, N1, N2), bernstein_basis, x, αs, dz, coeff_LE)

    # Upper and lower surface generation
    cam   = [ bernie(x, α_cam) for x ∈ xs ]
    thicc = [ bernie(x, α_thicc, dz_thicc) for x ∈ xs ]

    Foil(camber_thickness_to_coordinates(xs, cam, thicc), "Camber-Thickness CST")
end

"""
    coordinates_to_CST(coords, num_dvs)

Convert coordinates to a specified number of CST variables by performing a least-squares solution.
"""
function coordinates_to_CST(coords, num_dvs)
    xs       = @views coords[:,1]
    S_matrix = reduce(hcat, @. bernstein_class(xs, 0.5, 1.0) * bernstein_basis(xs, num_dvs - 1, i) for i in 0:num_dvs - 1)
    alphas   = @views S_matrix \ coords[:,2]
end

"""
    camber_thickness_to_CST(coords, num_dvs)

Convert camber-thickness coordinates to a specified number of CST variables by performing a least-squares solution.
"""
function camber_thickness_to_CST(coords, num_dvs)
    xs, camber, thickness = (columns ∘ coordinates_to_camber_thickness)(coords)

    alpha_cam   = coordinates_to_CST([ xs camber ], num_dvs)
    alpha_thick = coordinates_to_CST([ xs thickness ], num_dvs)

    alpha_cam, alpha_thick
end