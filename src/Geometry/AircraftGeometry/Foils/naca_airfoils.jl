## NACA Parametrisations
#==========================================================================================#

# NACA 4-digit parameter functions
naca4_thickness(t_by_c, xc, sharp_trailing_edge :: Bool) = 5 * t_by_c * (0.2969 * √xc - 0.1260 * xc - 0.3516 * xc^2 + 0.2843 * xc^3 - (ifelse(sharp_trailing_edge, 0.1036, 0.1015) * xc^4))

naca4_camberline(pos, cam, xc) = ifelse(
                                        xc < pos, 
                                        (cam / pos^2) * xc * (2 * pos - xc), 
                                        cam / (1 - pos)^2 * ( (1 - 2 * pos) + 2 * pos * xc - xc^2) 
                                       )

naca4_gradient(pos, cam, xc) = atan(2 * cam / (ifelse(xc < pos, pos^2, (1 - pos)^2)) * (pos - xc))

"""
    naca4_coordinates(digits :: NTuple{4, <: Real}, n :: Integer, sharp_trailing_edge :: Bool)

Generate the coordinates of a NACA 4-digit series profile with a specified number of points, and a Boolean flag to specify a sharp or blunt trailing edge.
"""
function naca4_coordinates(digits :: NTuple{4, <: Real}, n :: Integer, sharp_trailing_edge :: Bool)
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

"""
    naca4(digits :: NTuple{4, <: Real}, n :: Integer = 40; sharp_trailing_edge :: Bool)

Generate a `Foil` of a NACA 4-digit series profile with a specified number of points (40 by default), and a named option to specify a sharp or blunt trailing edge.
"""
naca4(digits :: NTuple{4, <: Real}, n = 40; sharp_trailing_edge = true) = Foil(naca4_coordinates(digits, n, sharp_trailing_edge), string("NACA ", digits...))

naca4(a, b, c, d, n = 40; sharp_trailing_edge = true) = naca4((a,b,c,d), n; sharp_trailing_edge = sharp_trailing_edge)
