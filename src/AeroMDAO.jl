module AeroMDAO

using LinearAlgebra
import Base: *, +
using Base.Iterators
using Interpolations
using Statistics
using DelimitedFiles

#--------------------HACKS----------------------#

# "Lenses" to access subfields on lists of objects
|>(obj, fields :: Array{Symbol}) = foldl(getproperty, fields, init = obj)
|>(list_objs :: Array{T}, fields :: Array{Symbol}) where {T <: Any} = list_objs .|> [fields]

# Convert 2D array to list of tuples
arraytolist(xs) = (collect ∘ zip)([ xs[:,n] for n in 1:length(xs[1,:])]...)

# Copying NumPy's linspace function
linspace(min, max, step) = min:(max - min)/step:max

#-------------HASKELL MASTER RACE--------------#

span(pred, iter) = (takewhile(pred, iter), dropwhile(pred, iter))
splitat(n, xs) = (xs[1:n,:], xs[n:end,:])  
lisa(pred, iter) = span(!pred, iter)

# Zieg Heil!

#----------------VECTOR SPACES?---------------#

# Convert struct entries to lists
structtolist(x) = [ getproperty(x, name) for name ∈ (fieldnames ∘ typeof)(x) ]

abstract type Point end

*(scale :: Real, point :: Point) = typeof(point)(scale * structtolist(point)...)
+(p1 :: Point, p2 :: Point) = typeof(p1)(structtolist(p1) .+ structtolist(p2)...)

struct Point2D <: Point
    x :: Real; y :: Real;
end

struct Point3D <: Point
    x :: Real; y :: Real; z :: Real;
end

#----------------PANEL METHODS------------------#

abstract type Panels end

struct Panel <: Panels
    loc :: Array{Point}
end

"""
Evaluates the midpoint of the coordinates of a panel.
"""
# midpoint(panel :: Panel) = 


abstract type Aircraft end

#----------------AIRFOIL----------------------#

"""
Airfoil structure consisting of foil coordinates as an array of points, a chord length, and twist angle in degrees.
"""
struct Foil <: Aircraft
    foil :: Union{Array{<: Real, 2}, Array{Point}} # The foil profile as an array of coordinates, must be in Selig format.
end

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
            i += 1
            if slope(x, y, xp, yp) >= slope(x, y, xn, yn)
                return splitat(i, coords)
            else
                return (reverse ∘ splitat)(i, coords)
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
    [x_circ y_circ]
end

"""
Discretises a geometry consisting of x and y coordinates into panels by projecting the x-coordinate of a circle onto the geometry.
"""
function cosine_foil(airfoil :: Foil, n :: Integer = 40)
    upper, lower = split_foil(airfoil.foil)
    upper_cos, lower_cos = cosine_interp.([reverse(upper, dims=1), lower], n)
    Panel([reverse(upper_cos, dims=1); lower_cos])
end

#-----------------WING---------------------#

"""
Definition for a wing consisting of airfoils, span lengths, dihedrals, and sweep angles.
"""
struct HalfWing <: Aircraft
    foils :: Array{Foil} # Airfoil profiles
    chords :: Array{<: Real} # Chord lengths (m)
    spans :: Array{<: Real}  # Leading-edge to leading-edge distance between foils (m)
    dihedrals :: Array{<: Real} # Dihedral angles (deg)
    sweeps :: Array{<: Real} # Leading-edge sweep angles (deg)
    twists :: Array{<: Real} # Twist angles (deg)
    HalfWing(foils, chords, spans, dihedrals, sweeps, twists) = new(foils, chords, spans, 
    deg2rad.(dihedrals), deg2rad.(sweeps), deg2rad.(twists)) # Convert to radians
end

aspect_ratio(span, area) = span^2 / area
mean_geometric_chord(span, area) = area / span
taper_ratio(root_chord, tip_chord) = tip_chord / root_chord
area(span, chord) = span * chord
mean_aerodynamic_chord(chord, taper_ratio) = (2/3) * chord * (1 + taper_ratio + taper_ratio^2)/(1 + taper_ratio)
quarter_chord(chord) = 0.25 * chord

# Difference operators on lists?
fwdsum(xs) = xs[2:end] .+ xs[1:end-1]
fwddiff(xs) = xs[2:end] .- xs[1:end-1]
cendiff(xs) = xs[3:end] .- 2 * xs[2:end-1] .+ xs[1:end-2] 


"""
Computes the span of a half-wing.
"""
span(wing :: HalfWing) = (sum ∘ map)(*, wing.spans, cos.(wing.dihedrals))

"""
Computes the projected area of a half-wing.
"""
function projected_area(wing :: HalfWing)
    mean_chords =  fwdsum(wing.chords) / 2 # Mean chord lengths of sections.
    mean_spans = fwdsum(wing.spans) / 2 # Mean span lengths of sections
    sum(mean_chords .* mean_spans)
end

"""
Computes the mean aerodynamic chord of a half-wing.
"""
function mean_aerodynamic_chord(wing :: HalfWing)
    mean_chords = fwdsum(wing.chords) / 2
    quarter_chords = quarter_chord.(abs.(fwddiff(wing.chords)))
    taper_ratios = map(/, wing.chords[2:end], wing.chords)
    areas = mean_chords .* wing.spans
    macs = mean_aerodynamic_chord(wing.chords)
    mac = (sum ∘ map)(*, macs, areas)/sum(areas)
end

struct Wing <: Aircraft
    left :: HalfWing
    Right :: HalfWing
end

span(wing :: Wing) = span(wing.left) + span(wing.right
)
projected_area(wing :: Wing) = projected_area(wing.left) + projected_area(wing.right)
mean_aerodynamic_chord(wing :: Wing) = (mean_aerodynamic_chord(wing.left) + mean_aerodynamic_chord(wing.right)) / 2
aspect_ratio(wing :: Wing) = aspect_ratio(span(wing), projected_area(wing))

# Kleisli-ish composition to transport global coordinates?
WingPos = Pair(Wing, Point3D)

end