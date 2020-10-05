module FoilParametrization

include("AeroMDAO.jl")
using Base.Math
using Base.Iterators
using Interpolations
using DelimitedFiles
using .AeroMDAO: Foil

# Copying NumPy's linspace function
linspace(min, max, step) = min:(max - min)/step:max

# Convert 2D array to list of tuples
arraytolist(xs) = (collect ∘ zip)([ xs[:,n] for n in 1:length(xs[1,:])]...)

#-------------HASKELL MASTER RACE--------------#

span(pred, iter) = (takewhile(pred, iter), dropwhile(pred, iter))
splitat(n, xs) = (xs[1:n,:], xs[n+1:end,:])  
lisa(pred, iter) = span(!pred, iter)

# Zieg Heil!

#------------FOIL PROCESSING------------#

"""
Reads a '.dat' file consisting of 2D coordinates, for an airfoil.
"""
function read_foil(path :: String)
    readdlm(path, skipstart = 1)
    # Point2D{Float64}.(f[:,1], f[:,2])
end

slope(x1, y1, x2, y2) = (y2 - y1)/(x2 - x1)

function split_foil(coords :: Array{<:Real, 2})
    cods = arraytolist(coords) # Convert to list of tuples
    for (i, ((xp, yp), (x, y), (xn, yn))) ∈ enumerate(zip(cods[1:end-2], cods[2:end-1], cods[3:end]))
        if x < xp && x < xn
            if slope(x, y, xp, yp) >= slope(x, y, xn, yn)
                return splitat(i, coords)
            else
                return splitat(i, reverse(coords))
            end
        end
    end
    return (coords, [])
    # span(((xp, yp), (x, y), (xn, yn)) -> x < xp && x < xn, zip(cods[1:end-2], cods[2:end-1], cods[3:end]))
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
    [x_circ y_circ]
end

"""
Discretises a geometry consisting of x and y coordinates into panels by projecting the x-coordinate of a circle onto the geometry.
"""
function cosine_foil(airfoil :: Array{<: Real, 2}, n :: Integer = 40)
    upper, lower = split_foil(airfoil)
    upper = [upper; lower[1,:]'] # Append leading edge from lower to upper
    upper_cos, lower_cos = cosine_interp(reverse(upper, dims=1), n), cosine_interp(lower, n)
    [reverse(upper_cos[2:end,:], dims=1); lower_cos]
end

#-------------------CST METHOD--------------------#

abstract type ClassShapeTransformation end

"""
Defines a general class shape transformation to represent a 2D geometric shape given a class and basis function representation.
"""
struct CSTBase <: ClassShapeTransformation
    class_func :: Function
    basis_func :: Function
end

# Basic shape function
function shape_function(cst :: CSTBase, x :: Real, coeffs :: Array{<: Real, 1}, coeff_LE :: Real = 0)
    n = length(coeffs)
    terms = [ cst.basis_func(x, n, i) for i in 0:n ]
    (sum ∘ map)(*, coeffs, terms) + coeff_LE * (x^0.5) * (1 - x)^(n - 0.5)
end

# Computing coordinates
function cst_coords(cst :: CSTBase, x :: Real, alphas :: Array{<: Real, 1}, dz :: Real, coeff_LE :: Real = 0)
    y = cst.class_func(x) * shape_function(cst, x, alphas, coeff_LE) + x * dz
    (x, y)
end

# Basis function setup
function airfoil_basis(coords :: Array{<: Real, 2}, basis_func :: Function, vec_size :: Integer = 40)
    
    num_points = len(coords[:, 0])
    xs = coords[:,0]
    ys = coords[:,1]
    t = linspace(0.0, 1.0, vec_size)

    basis = [ basis_func(num_points - 1, i, t) for i in 0:n ]

    xvals, yvals = xs .* basis, ys .* basis
    [xvals yvals]
end

#----------BERNSTEIN BASIS--------------#

# Bernstein basis for class function
bernstein_class(x, N1 = 0.5, N2 = 1) = x^N1 * (1 - x)^N2 

# Bernstein basis element
bernstein_basis(x, n, k) = binomial(n, k) * bernstein_class(x, k, n - k)

"""Defines a cosine-spaced airfoil using the Class Shape Transformation method on a Bernstein polynomial basis, with support for leading edge modifications.
"""
function airfoil_CST(alpha_u :: Array{<: Real}, alpha_l  :: Array{<: Real}, dz_u :: Real, dz_l :: Real, coeff_LE :: Real = 0, num_points :: Int = 40)
    # Instantiate airfoil method using Bernstein class function and basis by default
    temp_airfoil = CSTBase(bernstein_class, bernstein_basis)

    # Cosine spacing for airfoil of unit chord length
    xs = cosine_dist(0.5, 1, num_points)

    # Upper and lower surface generation
    upper_surf = cst_coords.(temp_airfoil, xs, alpha_u, dz_u, coeff_LE)
    lower_surf = cst_coords.(temp_airfoil, xs, alpha_l, dz_l, coeff_LE)

    # Zipping into coordinates and concatenating 
    # Following counter-clockwise ordering from trailing edge]
    [reverse([xs upper_surf], dims = 1); xs lower_surf]
end

#--------------CAMBER-THICKNESS REPRESENTATION----------------#

# Camber-thickness to coordinates and inverse transformations
"""
Converts coordinates to camber-thickness representation.
"""
function coords_to_camthick(coords :: Array{<: Real, 2})
    upper, lower = (split_foil ∘ cosine_foil)(coords)
    camber = (reverse(upper[:,1]) .+ lower[2:end,1]) / 2
    thickness = reverse(upper[:,1]) .- lower[2:end,1]
    [lower[:,1] [lower[1,2]; camber] [0; thickness]]
end

"""
Converts camber-thickness representation to coordinates.
"""
camthick_to_coords(xs, camber, thickness) = [reverse([xs camber .+ thickness / 2], dims = 1); xs camber .- thickness / 2]

#--------NACA PARAMETRIZATION------------#

function naca4(digits :: Tuple, n :: Integer = 40, closed_te = false, split = false)
    
    # Camber
    cam = digits[1] / 100
    # Position
    pos = digits[2] / 10
    # Thickness-to-chord ratio
    t_by_c = (10 * digits[3] + digits[4]) / 100

    # Cosine spacing
    xs = cosine_dist(0.5, 1.0, n)

    # Thickness distribution
    thickness = [ 5 * t_by_c * (0.2969 * √xc - 0.1260 * xc - 0.3516 * xc^2 + 0.2843 * xc^3 - (closed_te ? 0.1036 : 0.1015) * xc^4) for xc in xs ]
    
    if pos == 0.0 || cam == 0.0
        x_upper = xs
        y_upper = thickness
        x_lower = xs
        y_lower = -thickness
    else
        # Compute camberline
        camline = [ xc < pos ? (cam / pos^2) * xc * (2 * pos - xc) : cam / (1 - pos)^2 * ( (1 - 2 * pos) + 2 * pos * xc - xc^2) for xc in xs ]
        # Compute gradients
        gradients = [ atan(2 * cam / (xc < pos ? pos^2 : (1 - pos)^2) * (pos - xc)) for xc in xs ]
        # Upper surface
        x_upper = xs .- thickness .* sin.(gradients) 
        y_upper = camline .+ thickness .* cos.(gradients)
        # Lower surface
        x_lower = xs .+ thickness .* sin.(gradients) 
        y_lower = camline .- thickness .* cos.(gradients)
    end
    [reverse([x_upper y_upper], dims = 1); x_lower y_lower]
end

end