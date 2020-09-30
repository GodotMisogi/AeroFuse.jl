module FoilParametrization

include("AeroMDAO.jl")
using Base.Math
using Base.Iterators
using .AeroMDAO: cosine_dist

"""
Defines a general class shape transformation to represent a 2D geometric shape given a class and basis function representation.
"""
abstract type ClassShapeTransformation end

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
        # println(x_upper)
        y_upper = camline .+ thickness .* cos.(gradients)
        # Lower surface
        x_lower = xs .+ thickness .* sin.(gradients) 
        y_lower = camline .- thickness .* cos.(gradients)
    end
    [reverse([x_upper y_upper], dims = 1); x_lower y_lower]
end

end