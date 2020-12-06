include("foil_parametrization.jl")
# using .FoilParametrization: read_foil, cosine_foil, foil_camthick

#-------------------------AIRFOIL----------------------#

"""
Airfoil structure consisting of foil coordinates as an array of points.
"""
struct Foil <: Aircraft
    coords :: Array{<: Real, 2} # The foil profile as an array of coordinates, must be in Selig format.
end

"""
Scales the coordinates of a Foil, usually to some chord length.
"""
scale_foil(foil :: Foil, chord) = chord * foil.coords

"""
Translates the coordinates of a Foil by (x, y, z).
"""
shift_foil(foil :: Foil, x, y, z) = [ x y z ] .+ foil.coords

"""
Returns a Foil with cosine spacing for a given number of points. 
"""
cut_foil(foil :: Foil, num) = Foil(cosine_foil(foil.coords, n = num))

"""
Computes the camber-thickness distribution of a Foil with cosine spacing..
"""
camber_thickness(foil :: Foil, num) = Foil(foil_camthick(cosine_foil(foil.coords), num + 1))

"""
Projects a Foil onto the x-z plane at y = 0.
"""
coordinates(foil :: Foil) = [ foil.coords[:,1] (zeros ∘ length)(foil.coords[:,1]) foil.coords[:,2] ]