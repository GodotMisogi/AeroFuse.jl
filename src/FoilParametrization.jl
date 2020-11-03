module FoilParametrization

include("MathTools.jl")
using Base.Math
using Base.Iterators
using Interpolations
using DelimitedFiles
using .MathTools: slope, splitat, adj3


#-------------FOIL PROCESSING------------------#

"""
Reads a '.dat' file consisting of 2D coordinates, for an airfoil.
"""
function read_foil(path :: String; header = true)
    readdlm(path, skipstart = header ? 1 : 0)
    # Point2D{Float64}.(f[:,1], f[:,2])
end

function split_foil(coords)
    for (i, (xp, yp, x, y, xn, yn)) ∈ (enumerate ∘ eachrow ∘ adj3)(coords)
        if x < xp && x < xn
            if slope(x, y, xp, yp) >= slope(x, y, xn, yn)
                return splitat(i, coords)
            else
                return splitat(i, coords[end:-1:1,:])
            end
        end
    end
    (coords, [])
end

"""
Provides the projections to the x-axis for a circle of given diameter and center.
"""
cosine_dist(x_center :: Real, diameter :: Real, n :: Integer = 40) = x_center .+ (diameter / 2) * cos.(range(-π, stop = 0, length = n))

function cosine_interp(coords :: Array{<:Real, 2}, n :: Integer = 40)
    xs, ys = coords[:,1], coords[:,2]

    d = maximum(xs) - minimum(xs)
    x_center = (maximum(xs) + minimum(xs)) / 2
    x_circ = cosine_dist(x_center, d, n)
    
    itp_circ = LinearInterpolation(xs, ys)
    y_circ = itp_circ(x_circ)

    [ x_circ y_circ ]
end

"""
Discretises a foil profile into panels by projecting the x-coordinates of a circle onto the geometry.
"""
function cosine_foil(coords :: Array{<: Real, 2}; n :: Integer = 40)
    upper, lower = split_foil(coords)
    upper = [upper; lower[1,:]'] # Append leading edge from lower to upper
    upper_cos, lower_cos = cosine_interp(upper[end:-1:1,:], n), cosine_interp(lower, n)

    [ upper_cos[end:-1:2,:]; 
      lower_cos ]
end

#-------------------CST METHOD--------------------#

# Basic shape function
function shape_function(x :: Real, basis_func :: Function, coeffs :: Array{<: Real, 1}, coeff_LE :: Real = 0)
    n = length(coeffs)
    terms = [ basis_func(x, n, i) for i in 0:n-1 ]
    sum(coeffs .* terms) + coeff_LE * (x^0.5) * (1 - x)^(n - 0.5)
end

# Computing coordinates
cst_coords(class_func :: Function, basis_func :: Function, x :: Real, alphas :: Array{<: Real, 1}, dz :: Real, coeff_LE :: Real = 0) = class_func(x) * shape_function(x, basis_func, alphas, coeff_LE) + x * dz

#--------------BERNSTEIN BASIS-----------------#

"""
Bernstein basis for class function.
"""
bernstein_class(x, N1 = 0.5, N2 = 1) = x^N1 * (1 - x)^N2 

"""
Bernstein basis element.
"""
bernstein_basis(x, n, k) = binomial(n, k) * bernstein_class(x, k, n - k)

"""
Defines a cosine-spaced airfoil using the Class Shape Transformation method on a Bernstein polynomial basis, with support for leading edge modifications.
"""
function kulfan_CST(alphas :: Array{<: Real, 2}, (dz_u, dz_l), coeff_LE :: Real = 0, num_points :: Integer = 40)

    # Cosine spacing for airfoil of unit chord length
    xs = cosine_dist(0.5, 1, num_points)

    # λ-function for Bernstein polynomials
    bernie = (x, alphas, dz) -> cst_coords(bernstein_class, bernstein_basis, x, alphas, dz, coeff_LE)

    # Upper and lower surface generation
    upper_surf = [ bernie(x, alphas[:,1], dz_u) for x ∈ xs ]
    lower_surf = [ bernie(x, alphas[:,2], dz_l) for x ∈ xs ]

    # Counter-clockwise ordering
    [ [xs upper_surf][end:-1:2,:]; 
       xs lower_surf ]
end

#--------------CAMBER-THICKNESS REPRESENTATION----------------#


"""
Converts an airfoil to its camber-thickness representation in cosine spacing.
"""
function foil_camthick(coords :: Array{<: Real, 2}, num :: Integer = 40)
    upper, lower = split_foil(cosine_foil(coords, n = num))

    xs, y_LE = lower[:,1], lower[1,2]   # Getting abscissa and leading edge ordinate
    y_upper, y_lower = upper[end:-1:1,2], lower[2:end,2] # Excluding leading edge point

    camber = [y_LE; (y_upper .+ y_lower) / 2]
    thickness = [0; y_upper .- y_lower]

    [ xs camber thickness ]
end

"""
Converts the camber-thickness representation to a Foil.
"""
camthick_foil(xs, camber, thickness) = [ [xs camber .+ thickness / 2][end:-1:1,:]; xs camber .- thickness / 2 ]

#--------NACA PARAMETRIZATION------------#

# NACA 4-digit parameter functions
naca4_thickness(t_by_c, xc, sharp_trailing_edge :: Bool) = 5 * t_by_c * (0.2969 * √xc - 0.1260 * xc - 0.3516 * xc^2 + 0.2843 * xc^3 - (sharp_trailing_edge ? 0.1036 : 0.1015) * xc^4)
naca4_camberline(pos, cam, xc) = xc < pos ? (cam / pos^2) * xc * (2 * pos - xc) : cam / (1 - pos)^2 * ( (1 - 2 * pos) + 2 * pos * xc - xc^2)
naca4_gradient(pos, cam, xc) = atan(2 * cam / (xc < pos ? pos^2 : (1 - pos)^2) * (pos - xc))

function naca4(digits :: Tuple, n :: Integer = 40; sharp_trailing_edge :: Bool = false)
    
    # Camber
    cam = digits[1] / 100
    # Position
    pos = digits[2] / 10
    # Thickness-to-chord ratio
    t_by_c = (10 * digits[3] + digits[4]) / 100

    # Cosine spacing
    xs = cosine_dist(0.5, 1.0, n)

    # Thickness distribution
    thickness = [ naca4_thickness(t_by_c, xc, sharp_trailing_edge) for xc in xs ]
    
    if pos == 0 || cam == 0
        x_upper = xs
        y_upper = thickness
        x_lower = xs
        y_lower = -thickness
    else
        # Compute camberline
        camber = [ naca4_camberline(pos, cam, xc) for xc in xs ]
        # Compute gradients
        gradients = [ naca4_gradient(pos, cam, xc) for xc in xs ]
        # Upper surface
        x_upper = xs .- thickness .* sin.(gradients) 
        y_upper = camber .+ thickness .* cos.(gradients)
        # Lower surface
        x_lower = xs .+ thickness .* sin.(gradients) 
        y_lower = camber .- thickness .* cos.(gradients)
    end
    [ [x_upper y_upper][end:-1:2,:]; 
       x_lower y_lower ]
end

end