include("../src/AeroMDAO.jl")

using .AeroMDAO: read_foil, Point3D, Point2D, WingSection, Wing, projected_area, Foil, cosine_foil

foilpath = "airfoil_database/a18.dat"

# Wing section setup
num_secs = 5
xs = zeros(num_secs)
ys = 0:2:num_secs
zs = zeros(num_secs)

chords = repeat([1.0], num_secs)    # Chord lengths
twists = zeros(num_secs)            # Twists
coords = read_foil(foilpath)
foils = [ coords for n in 1:num_secs ]
airfoils = Foil.(foils, chords) # Airfoils
# println(airfoils)

foil = Foil(coords, 1.0)
println(cosine_foil(foil))


# secs = [ WingSection(x...) for x in zip(locs, chords, twists, foils) ]

# wing_loc = Point3D{Float64}(0.0, 0.0, 0.0)
# wing = Wing(wing_loc, secs)

# println(projected_area(wing))