include("../src/AeroMDAO.jl")

using .AeroMDAO: read_foil, Point3D, Point2D, WingSection, Wing

foilpath = "airfoil_database/a18.dat"

coords = let f = read_foil(foilpath); [ Point2D{Float64}(x...) for x in zip(f[:,1], f[:,2]) ] end
num_secs = 5
xs = zeros(num_secs)
ys = range(0, stop = 2, length = num_secs)
zs = zeros(num_secs)

locs = Point3D{Float64} 
chords = repeat([1.0], num_secs)
twists = zeros(num_secs)
foils = repeat(coords, num_secs)

secs = [ WingSection(x...) for x in zip(locs, chords, twists, foils) ]

wing = Wing(Point3D{Float64}(0.0, 0.0, 0.0), secs)

println(projected_area(wing))