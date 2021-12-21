## Foil type
#==========================================================================================#

abstract type AbstractFoil end

"""
    Foil(coords)

Airfoil structure consisting of foil coordinates as an array of points. Should be in Selig format for compatibility with other AeroMDAO tools.
"""
struct Foil{T <: Real} <: AbstractFoil
    x      :: Vector{T}
    y      :: Vector{T}
    name   :: String
end

coordinates(foil :: Foil) = [ foil.x foil.y ]

Foil(coords, name = "Unnamed") where T <: Real = Foil{T}(getindex.(coords, 1), getindex.(coords, 2), name)
Foil(coords :: AbstractMatrix{T}, name = "Unnamed") where T <: Real = @views Foil{T}(coords[:,1], coords[:,2], name)

arc_length(foil :: Foil) = let c = coordinates(foil); @views norm(c[2:end,:] .- c[1:end-1,:]) end

"""
    scale_foil(foil :: Foil, chord)

Scale the coordinates of a Foil, usually to some chord length.
"""
scale_foil(foil :: Foil, chord) = chord .* coordinates(foil)

"""
    cosine_foil(foil :: Foil, num :: Integer)

Return a Foil with cosine spacing for a given number of points.
"""
cosine_foil(foil :: Foil, num :: Integer) = cosine_foil(coordinates(foil), num)

"""
    camber_thickness(foil :: Foil, num :: Integer)

Compute the camber-thickness distribution of a Foil with cosine spacing.
"""
camber_thickness(foil :: Foil, num :: Integer = 40) = foil_camthick(cosine_foil(coordinates(foil)), num + 1)

function max_thickness_to_chord_ratio_location(coords)
    @views xs, thiccs = coords[:,1], coords[:,3]
    max_thick_arg = argmax(thiccs)
    @views xs[max_thick_arg], thiccs[max_thick_arg]
end

## Foil processing
#==========================================================================================#

"""
    read_foil(path :: String; header = true)

Read a '.dat' file consisting of 2D coordinates, for an airfoil as an array of `SVector`s, with an optional argument to skip the header.
"""
function read_foil(path :: String; header = true)
    coords = readdlm(path, skipstart = header ? 1 : 0)
    # @views @SMatrix [ coords[:,1] coords[:,2] ]
end

function split_foil(coords)
    # display(coords)
    for (i, ((xp, yp), (x, y), (xn, yn))) ∈ (enumerate ∘ adj3 ∘ collect ∘ eachrow)(coords)
        if x < xp && x < xn
            return ifelse(slope(x, y, xp, yp) >= slope(x, y, xn, yn), splitat(i, coords), splitat(i, coords[end:-1:1,:]))
        end
    end
    (coords, [])
end

function paneller(foil :: Foil, num_panels :: Integer)
    coords = cosine_foil(coordinates(foil), Int(ceil(num_panels / 2)))
    vecs   = SVector.(coords[:,1], coords[:,2])
    @views Panel2D.(vecs[2:end,:], vecs[1:end-1,:])[end:-1:1]
end

"""
Discretises a foil profile into panels by projecting the x-coordinates of a circle onto the geometry.
"""
function cosine_foil(coords, n :: Integer = 40)
    upper, lower = split_foil(coords)
    n_upper = @views [ upper       ;
                       lower[1,:]' ] # Append leading edge point from lower to upper

    upper_cos = @views cosine_interp(n_upper[end:-1:1,:], n)
    lower_cos = cosine_interp(lower, n)

    @views [ upper_cos[end:-1:2,:] ;
             lower_cos             ]
end

## Class shape transformation method
#==========================================================================================#

# Basic shape function
function shape_function(x, basis_func, coeffs, coeff_LE = 0)
    n     = length(coeffs)
    terms = basis_func.(x, n - 1, 0:n-1)
    dot(coeffs, terms) + coeff_LE * (x^0.5) * (1 - x)^(n - 0.5)
end

# Computing coordinates
cst_coords(class_func, basis_func, x, alphas, dz, coeff_LE, args...) = class_func(x) * shape_function(x, basis_func, alphas, coeff_LE) + x * dz

## Bernstein basis
#==========================================================================================#

bernstein_class(x, N1, N2) = x^N1 * (1 - x)^N2
bernstein_basis(x, n, k)   = binomial(n, k) * bernstein_class(x, k, n - k)

"""
    kulfan_CST(alpha_u, alpha_l,
               (dz_u, dz_l) :: NTuple{2, Real},
               coeff_LE = 0,
               n :: Integer = 40)

Define a cosine-spaced airfoil with ``2n`` points using the Class Shape Transformation method on a Bernstein polynomial basis with arrays of coefficients ``(\\alpha_u,~ \\alpha_l)`` for the upper and lower surfaces, trailing-edge spacing values ``(\\Delta z_u,~ \\Delta z_l)``, and support for leading edge modifications.
"""
function kulfan_CST(alpha_u, alpha_l, (dz_u, dz_l) = (0., 0.), (LE_u, LE_l) = (0., 0.), n :: Integer = 40, N1 = 0.5, N2 = 1.)
    # Cosine spacing for airfoil of unit chord length
    xs = cosine_spacing(0.5, 1, n)

    # λ-function for Bernstein polynomials
    bernie(x, alphas, dz, LE) = cst_coords(y -> bernstein_class(y, N1, N2), bernstein_basis, x, alphas, dz, LE)

    # Upper and lower surface generation
    upper_surf = [ bernie(x, alpha_u, dz_u, LE_u) for x ∈ xs ]
    lower_surf = [ bernie(x, alpha_l, dz_l, LE_l) for x ∈ xs ]

    # Counter-clockwise ordering
    @views [ xs[end:-1:2] upper_surf[end:-1:2] ;
             xs           lower_surf           ]
end

function camber_CST(alpha_cam, alpha_thicc, dz_thicc = 0., coeff_LE = 0, n :: Integer = 40, N1 = 0.5, N2 = 1.)
    # Cosine spacing for airfoil of unit chord length
    xs = cosine_spacing(0.5, 1, n)

    # λ-function for Bernstein polynomials
    bernie(x, alphas, dz = 0.) = cst_coords(y -> bernstein_class(y, N1, N2), bernstein_basis, x, alphas, dz, coeff_LE)

    # Upper and lower surface generation
    cam   = [ bernie(x, alpha_cam) for x ∈ xs ]
    thicc = [ bernie(x, alpha_thicc, dz_thicc) for x ∈ xs ]

    camber_thickness_to_coordinates(xs, cam, thicc)
end

function coords_to_CST(coords, num_dvs)
    xs       = coords[:,1]
    S_matrix = reduce(hcat, @. bernstein_class(xs, 0.5, 1.0) * bernstein_basis(xs, num_dvs - 1, i) for i in 0:num_dvs - 1)
    alphas   = S_matrix \ coords[:,2]
end

function camthick_to_CST(coords, num_dvs)
    xs, camber, thickness = (columns ∘ foil_camthick)(coords)

    alpha_cam   = coords_to_CST([ xs camber ], num_dvs)
    alpha_thick = coords_to_CST([ xs thickness ], num_dvs)

    alpha_cam, alpha_thick
end

## Camber-thickness representation
#==========================================================================================#

"""
Converts an airfoil to its camber-thickness representation in cosine spacing.
"""
function foil_camthick(coords, num :: Integer = 40)
    upper, lower = split_foil(cosine_foil(coords, num))

    xs, y_LE         = lower[:,1], lower[1,2]            # Getting abscissa and leading edge ordinate
    y_upper, y_lower = upper[end:-1:1,2], lower[2:end,2] # Excluding leading edge point

    camber       = [ y_LE; (y_upper + y_lower) / 2 ]
    thickness    = [ 0; y_upper - y_lower          ]

    [ xs camber thickness ]
end

"""
    camber_thickness_to_coordinates(xs, camber, thickness)

Converts the camber-thickness representation to coordinates given the ``x``-locations and their corresponding camber and thickness values.
"""
function camber_thickness_to_coordinates(xs, camber, thickness)
    coords = @views [ [xs camber + thickness / 2][end:-1:2,:];
                       xs camber - thickness / 2             ]
    # @views SVector.(coords[:,1], coords[:,2])
end

camber_thickness_to_coordinates(coords) = @views camber_thickness_to_coordinates(coords[:,1], coords[:,2], coords[:,3])

"""
    camber_coordinates(coords :: Array{2, <: Real})

Provides the camber coordinates on the ``x``-``z`` plane at ``y = 0``.
"""
camber_coordinates(coords) = @views [ coords[:,1] zeros(length(coords[:,1])) coords[:,2] ]

"""
    thickness_coordinates(coords :: Array{2, <: Real})

Provides the thickness coordinates on the ``x``-``z`` plane at ``y = 0``.
"""
thickness_coordinates(coords) = @views [ coords[:,1] zeros(length(coords[:,1])) coords[:,3] ]


## NACA Parametrisation
#==========================================================================================#

# NACA 4-digit parameter functions
naca4_thickness(t_by_c, xc, sharp_trailing_edge :: Bool) = 5 * t_by_c * (0.2969 * √xc - 0.1260 * xc - 0.3516 * xc^2 + 0.2843 * xc^3 - (ifelse(sharp_trailing_edge, 0.1036, 0.1015) * xc^4))
naca4_camberline(pos, cam, xc) = ifelse(xc < pos, (cam / pos^2) * xc * (2 * pos - xc), cam / (1 - pos)^2 * ( (1 - 2 * pos) + 2 * pos * xc - xc^2) )
naca4_gradient(pos, cam, xc) = atan(2 * cam / (ifelse(xc < pos, pos^2, (1 - pos)^2)) * (pos - xc))

"""
    naca4(digits :: NTuple{4, <: Real}, n :: Integer; sharp_trailing_edge :: Bool)

Defines a NACA 4-digit series profile.
"""
function naca4(digits :: NTuple{4, <: Real}, n :: Integer = 40; sharp_trailing_edge = true)
    # Camber
    cam = digits[1] / 100
    # Position
    pos = digits[2] / 10
    # Thickness-to-chord ratio
    t_by_c = (10 * digits[3] + digits[4]) / 100

    # Cosine spacing
    xs = cosine_spacing(0.5, 1.0, n)

    # Thickness distribution
    thickness = naca4_thickness.(Ref(t_by_c), xs, Ref(sharp_trailing_edge))
    if pos == 0 || cam == 0
        x_upper = xs
        y_upper = thickness
        x_lower = xs
        y_lower = -thickness
    else
        # Compute camberline
        camber  = naca4_camberline.(Ref(pos), Ref(cam), xs)
        # Compute gradients
        grads   = naca4_gradient.(Ref(pos), Ref(cam), xs)
        # Upper surface
        x_upper = @. xs - thickness * sin(grads)
        y_upper = @. camber + thickness * cos(grads)
        # Lower surface
        x_lower = @. xs + thickness * sin(grads)
        y_lower = @. camber - thickness * cos(grads)
    end
    coords = [ [x_upper y_upper][end:-1:2,:];
                x_lower y_lower             ]
end

naca4(a, b, c, d, n = 40; sharp_trailing_edge = true) = naca4((a,b,c,d), n; sharp_trailing_edge = sharp_trailing_edge)
