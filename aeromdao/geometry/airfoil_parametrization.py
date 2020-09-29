import numpy as np
import operator as op
from functools import reduce

# Haskellian functions
def zip_with(func, *args): return [ func(*entry) for entry in zip(*args) ]

###---------------------------------------------------------###

def comb(n, r):
    r = min(r, n-r)
    numer = reduce(op.mul, range(n, n-r, -1), 1)
    denom = reduce(op.mul, range(1, r+1), 1)
    return numer // denom

# File I/O
def read_foil(airfoil_path):
    # Import airfoil coordinates
    with open(airfoil_path,'r') as f:
        next(f)
        coords = np.array([ list(map(float,line.split())) for line in f.readlines() ])

    return coords

def write_foil(coords, foil_name, filename='foil.dat'):
    np.savetxt(filename, coords, fmt='%.6f', delimiter='    ', header=foil_name, encoding='UTF-8')

###---------------------------------------------------------###

# Processing coordinates

def slope(c1, c2): return (c2[1] - c1[1])/(c2[0] - c1[0]) # m = (y₂ - y₁)/(x₂ - x₁)

# Splitting airfoil at origin into upper and lower surfaces
def split_foil(coords):
    for (i, (coord_prev, coord, coord_next)) in enumerate(zip(coords[:-2], coords[1:-1], coords[2:])):
        if coord < coord_prev and coord < coord_next:
            i += 1
            if slope(coord, coord_prev) >= slope(coord, coord_next): # Anticlockwise ordering
                return coords[:i+1], coords[i+1:]
            else: # Clockwise ordering
                return coords[i:], coords[:i]
    raise Exception("Turning point of leading edge not found.")

# Cosine spacing functions
def cosine_dist(x_center, diameter, n=40):
    "Provides the x-locations for a cosine spacing with a number of points."
    x_circ = np.cos(np.linspace(-np.pi, 0, n))
    x_scaled = (diameter/2)*x_circ + x_center

    return x_scaled

def cosine_interp(coords, **kwargs):
    "Interpolates a list of coordinates to a cosine spacing normalized to [0,1]."
    # Sorting and segregating coordinates
    coords = coords[coords[:,0].argsort()]
    x = coords[:,0] 
    y = coords[:,1]

    d = (x.max() - x.min())   # Diameter of circle/chord length
    x_c = (x.max() + x.min())/2    # Centre of circle
    x_interp = cosine_dist(x_c, d, **kwargs) # Scale x_interp
    y_interp = np.interp(x_interp, x, y)

    return np.column_stack((x_interp, y_interp))
    
def cosine_airfoil(upper, lower, **kwargs):
    "Interpolates an airfoil consisting of an upper and lower surface to a cosine spacing."

    cos_u = cosine_interp(upper, **kwargs)
    cos_l = cosine_interp(lower, **kwargs)

    return np.append(cos_l[::-1], cos_u, axis=0)

# Camber-thickness to coordinates and inverse transformations
def coords_to_camthick(coords, **kwargs):
    "Converts coordinates to camber-thickness representation."
    upper, lower = split_foil(cosine_airfoil(*split_foil(coords), **kwargs))
    n_upper, n_lower = upper[upper[:,0].argsort()], lower[lower[:,0].argsort()]
    camber, thickness = np.zeros(len(n_upper)), np.zeros(len(n_upper))
    camber[:] = (n_upper[:,1] + n_lower[:,1])/2
    thickness[:] = n_upper[:,1] - n_lower[:,1]

    return upper[:,0], camber, thickness

def camthick_to_coords(xs, camber, thickness):
    "Converts camber-thickness representation to coordinates."
    upper = np.column_stack((xs, camber + thickness/2))
    lower = np.column_stack((xs, camber - thickness/2))
    return np.append(upper[::-1], lower, axis=0)

###---------------------------------------------------------###

# CST Method setup

# Bernstein basis for class function
def bernstein_class(x, N1=0.5, N2=1): return x**N1 * (1 - x)**N2 

# Bernstein basis element
def bernstein_basis(x, n, i): return comb(n, i) * bernstein_class(x, N1=i, N2=n-i)

# Basis function setup
def airfoil_basis(coords, basis_func, vec_size=200, **kwargs):
    
    num_points = len(coords[:, 0])
    xs = coords[:,0]
    ys = coords[:,1]
    t = np.linspace(0.0, 1.0, vec_size)

    basis = np.array([ basis_func(num_points-1, i, t) for i in range(num_points) ])

    xvals = np.dot(xs, basis)
    yvals = np.dot(ys, basis)

    return xvals, yvals

# Class Shape Transformation (CST) Method
class CSTMethod():
    """
    Defines a general class shape transformation to represent a 2D geometric shape given a class and basis function representation.
    """
    def __init__(self, class_func, basis_func):
        self.class_func = class_func
        self.basis_func = basis_func

    # Shape function
    def shape_function(self, x, coeffs):
        n = len(coeffs)
        terms = [ self.basis_func(x, n, i) for i in range(n) ]
        return sum(zip_with(lambda x,y: x * y, coeffs, terms))

    # Computing z-coordinate
    def z(self, x, alphas, dz, **kwargs): 
        return self.class_func(x, **kwargs) * self.shape_function(x, alphas) + x * dz


# CST Leading Edge Modification for airfoils
class CSTAirfoilLEM(CSTMethod):
    
    def __init__(self, class_func, basis_func, coeffLE):
        super().__init__(class_func, basis_func)
        self.coeffLE = coeffLE 

    # Leading edge modification
    def shape_function(self, x, coeffs):
        n = len(coeffs)
        terms = [ self.basis_func(x, n, i) for i in range(n) ]

        return (sum(zip_with(lambda x,y: x * y, coeffs, terms))
                + self.coeffLE * (x**0.5) * (1 - x)**(n - 0.5))

def airfoil_CST(alpha_u, alpha_l, dz_u, dz_l, le=0, num_points=51):
    "Defines a cosine-spaced airfoil using the Class Shape Transformation method, with support for leading edge modifications."

    # Instantiate airfoil method using Bernstein class function and basis by default
    temp_airfoil = CSTAirfoilLEM(bernstein_class, bernstein_basis, le)

    # Cosine spacing for airfoil of unit chord length
    xs = cosine_dist(0.5, 1, num_points)

    # Upper and lower surface generation
    upper_surf = [ temp_airfoil.z(x, alpha_u, dz_u) for x in xs ]
    lower_surf = [ temp_airfoil.z(x, alpha_l, dz_l) for x in xs ]

    # Zipping into coordinates and concatenating 
    # Following counter-clockwise ordering from trailing edge
    upper = np.column_stack((xs, upper_surf))[::-1]
    lower = np.column_stack((xs, lower_surf))
    foil = np.concatenate((upper[:-1], lower))

    return foil

def coords_to_CST(coords):
    num_dvs = len(coords[:,0]) # Number of design variables
    
    S_matrix = np.zeros((num_dvs, num_dvs))
    for i in range(num_dvs):
        S_matrix[:,i] = bernstein_class(coords[:,0])*bernstein_basis(coords[:,0], num_dvs, i)

    alphas = np.linalg.solve(S_matrix, coords[:,1]) 
    
    return alphas

###---------------------------------------------------------###

# Reynolds number calculator
def reynolds_number(chord_length, freestream_velocity, rho = 1.225, mu = 1.7894e-5, water=False):
    if water:
        rho = 998.2
        mu = 8.899e-4
    return rho * freestream_velocity * chord_length / mu


# NACA 4-digit series airfoil generator
def NACA4(digits, c, n, closed_te = False, split = False):
    
    # Airfoil characteristics
    # Camber
    m = digits[0] / 100
    # Position
    p = digits[1] / 10
    # Thickness-to-chord ratio
    t_by_c = (10 * digits[2] + digits[3]) / 100

    # Cosine spacing
    xs = cosine_dist(c/2, c, n)
    x_upper, x_lower, y_upper, y_lower = np.zeros(len(xs)), np.zeros(len(xs)),np.zeros(len(xs)), np.zeros(len(xs))

    # Thickness distribution
    thickness = np.array([ 5 * t_by_c * c * (0.2969 * np.sqrt(xc / c) - 0.126 * xc / c - 0.3516 * (xc / c)**2 + 0.2843 * (xc / c)**3 - (0.1036 if closed_te else 0.1015) * (xc / c)**4) for xc in xs ])

    if p == 0.0 or m == 0.0:
        x_upper = xs
        y_upper = thickness
        x_lower = xs
        y_lower = -thickness
    else:
        # Compute m
        mline = [ (m / p**2) * xc * (2 * p - xc / c) if xc < p * c else (m * (c - xc) / (1 - p)**2) * (1 + xc / c - 2 * p) for xc in xs ]
        # Compute gradients
        gradients = [ np.arctan((2 * m / p**2) * (p - xc / c)) if xc < p * c else np.arctan((2 * m / (1 - p)**2) * (p - xc / c)) for xc in xs ] 
        # Compute surfaces
        x_upper = xs - thickness + np.sin(gradients)
        y_upper = mline + thickness + np.cos(gradients)
        x_lower = xs + thickness + np.sin(gradients)
        y_lower = mline - thickness + np.cos(gradients)

    if split:
        upper = np.column_stack((x_upper, y_upper))
        lower = np.column_stack((x_lower, y_lower))

        return upper[::-1], lower

    else:
        X, Y = np.append(x_upper[::-1], x_lower), np.append(y_upper[::-1], y_lower)

        return np.column_stack((X, Y))