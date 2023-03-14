## Foil type
#==========================================================================================#

abstract type AbstractFoil end

"""
    Foil(x, y, name = "")
    Foil(coordinates, name = "")

Structure consisting of foil coordinates in 2 dimensions with an optional name. 

The coordinates should be provided in counter-clockwise format, viz. from the trailing edge of the upper surface to the trailing edge of the lower surface.
"""
struct Foil{T <: Real} <: AbstractFoil
    x      :: Vector{T}
    y      :: Vector{T}
    name   :: String
end

function Foil(xs, ys, name = "")
    @assert length(xs) == length(ys) "Lengths of xs and ys must match!"

    T = promote_type(eltype(xs), eltype(ys))
    Foil{T}(xs, ys, name) 
end

Foil(coords :: Vector{<: FieldVector{2,<: Real}}, name = "") = Foil(getindex.(coords, 1), getindex.(coords, 2), name)
Foil((name, coords)) = Foil(coords, name)

function Foil(coords :: AbstractMatrix{<: Real}, name = "") 
    @assert size(coords, 2) == 2 "The array must have only two columns for coordinates!"

    @views Foil(coords[:,1], coords[:,2], name)
end

function Base.show(io :: IO, foil :: Foil)
    print(io, foil.name, " ", typeof(foil), " with ", length(foil.x), " points.")
end

## Foil processing
#==========================================================================================#

"""
    coordinates(foil :: Foil)

Generate the array of `Foil` coordinates. 
"""
coordinates(foil :: Foil) = [ foil.x foil.y ]

"""
    arc_length(foil :: Foil)

Compute the arc-length of a `Foil`.
"""
arc_length(foil :: Foil) = let c = coordinates(foil); @views norm(c[2:end,:] .- c[1:end-1,:]) end
chord_length(foil :: Foil) = maximum(foil.x)

"""
    scale(foil :: Foil, scale)

Scale the coordinates of a `Foil` to a scaling value.
"""
scale(foil :: Foil, scale) = Foil(scale .* coordinates(foil), foil.name)

translate(foil :: Foil; vector) = Foil(foil.x .+ vector[1], foil.y .+ vector[2])

@views function rotate(foil :: Foil; angle :: Real, center = zeros(2))
    T = promote_type(eltype(angle), eltype(center))
    trans  = @views [ foil.x .- center[1] foil.y .- center[2] ]     # Translate
    rotate = trans * RotMatrix{2,T}(deg2rad(angle))'    # Rotate
    Foil(rotate[:,1] .+ center[1], rotate[:,2] .+ center[2], foil.name) # Inverse translate
end

@views function interpolate(foil :: Foil, xs)
    upper, lower = split_surface(foil)

    y_u = linear_interpolation(upper[:,1], upper[:,2]).(xs)
    y_l = linear_interpolation(lower[:,1], lower[:,2]).(xs)

    Foil([ xs[end:-1:2]; xs ], [ y_u[end:-1:2]; y_l ], foil.name)
end

reflect(foil :: Foil) = setproperties(foil, y = -foil.y, name = "Inverted " * foil.name)

affine(foil :: Foil; angle, vector) = translate(rotate(foil; angle = angle); vector = vector)

"""
    camber_thickness(foil :: Foil, num :: Integer)

Compute the camber-thickness distribution of a `Foil` with cosine spacing.
"""
camber_thickness(foil :: Foil, num = 40) = coordinates_to_camber_thickness(foil, num + 1)

leading_edge_index(foil :: Foil) = argmin(coordinates(foil)[:,1])

"""
    upper_surface(foil :: Foil)

Get the upper surface coordinates of a `Foil` from leading to trailing edge.
"""
@views upper_surface(foil :: Foil) = reverse(coordinates(foil)[1:leading_edge_index(foil),:], dims = 1)

"""
    lower_surface(foil :: Foil)

Get the lower surface coordinates of a `Foil` from leading to trailing edge.
"""
@views lower_surface(foil :: Foil) = coordinates(foil)[leading_edge_index(foil):end,:]

"""
    split_surface(foil :: Foil)

Split the `Foil` coordinates into upper and lower surfaces.
"""
split_surface(foil :: Foil) = upper_surface(foil), lower_surface(foil)

"""
    cosine_interpolation(foil :: Foil, n :: Integer = 40)

Interpolate a `Foil` profile's coordinates to a cosine by projecting the x-coordinates of a circle onto the geometry with ``2n`` points.
"""
function cosine_interpolation(foil :: Foil, n :: Integer = 40)
    x_min, x_max = extrema(foil.x)
    x_circ = cosine_spacing((x_min + x_max) / 2, x_max - x_min, n)

    interpolate(foil, x_circ)
end

function camber_line(foil :: Foil, n = 60)
    upper, lower = split_surface(foil)
    xs  = LinRange(minimum(foil.x), maximum(foil.x), n + 1)
    y_u = @views linear_interpolation(upper[:,1], upper[:,2]).(xs)
    y_l = @views linear_interpolation(lower[:,1], lower[:,2]).(xs)

    [ xs (y_u + y_l) / 2 ]
end


"""
    make_panels(foil :: Foil)
    make_panels(foil :: Foil, n :: Integer)

Generate a vector of `Panel2D`s from a `Foil`, additionally with cosine interpolation using (approximately) ``n`` points if provided.
"""
make_panels(foil :: Foil) = @views Panel2D.(foil.x[2:end], foil.y[2:end], foil.x[1:end-1], foil.y[1:end-1])[end:-1:1]

make_panels(foil :: Foil, n :: Integer) = make_panels(cosine_interpolation(foil, n ÷ 2))

function camber(foil :: Foil, x_by_c)
    upper, lower = split_surface(foil)
    y_u = @views linear_interpolation(upper[:,1], upper[:,2])(x_by_c * chord_length(foil))
    y_l = @views linear_interpolation(lower[:,1], lower[:,2])(x_by_c * chord_length(foil))

    (y_u + y_l) / 2
end

function control_surface(foil :: Foil, δ, xc_hinge)
    y_hinge  = camber(foil, xc_hinge)
    rot_foil = rotate(foil; angle = -δ, center = [xc_hinge, y_hinge])
    coords   = coordinates(foil)
    @views coords[foil.x .>= xc_hinge,1] = rot_foil.x[foil.x .>= xc_hinge]
    @views coords[foil.x .>= xc_hinge,2] = rot_foil.y[foil.x .>= xc_hinge]
    Foil(coords, foil.name * " Deflected $(δ)° at $xc_hinge (x/c)")
end

"""
    control_surface(foil :: Foil; angle, hinge)

Modify a `Foil` to mimic a control surface by specifying a deflection angle (in degrees, clockwise-positive convention) and a normalized hinge ``x``-coordinate ``∈ [0,1]`` in terms of the chord length.
"""
control_surface(foil :: Foil; angle, hinge) = control_surface(foil, angle, hinge)

maximum_thickness_to_chord(foil :: Foil, n = 40) = maximum_thickness_to_chord(coordinates_to_camber_thickness(foil, n))

## Camber-thickness representation
#==========================================================================================#

"""
    coordinates_to_camber_thickness(coords, n = 40)

Convert 2-dimensional coordinates to its camber-thickness representation after cosine interpolation with ``2n`` points.
"""
function coordinates_to_camber_thickness(foil :: Foil, n = 40)
    # Cosine interpolation and splitting
    upper, lower = split_surface(cosine_interpolation(foil, n))

    camber       = @views (upper[:,2] + lower[:,2]) / 2
    thickness    = @views  upper[:,2] - lower[:,2]

    @views [ upper[:,1] camber thickness ]
end

"""
    camber_thickness_to_coordinates(xs, camber, thickness)

Convert the camber-thickness representation to 2-dimensional coordinates given the ``x``-locations and their corresponding camber and thickness values.
"""
camber_thickness_to_coordinates(xs, camber, thickness) = 
    @views [ [xs camber + thickness / 2][end:-1:2,:];
              xs camber - thickness / 2             ]

camber_thickness_to_coordinates(coords) = @views camber_thickness_to_coordinates(coords[:,1], coords[:,2], coords[:,3])

"""
    camber_coordinates(coords :: Array{2, <: Real})

Generate the camber coordinates on the ``x``-``z`` plane at ``y = 0``.
"""
camber_coordinates(coords) = @views [ coords[:,1] zero(coords[:,1]) coords[:,2] ]

"""
    thickness_coordinates(coords :: Array{2, <: Real})

Generate the thickness coordinates on the ``x``-``z`` plane at ``y = 0``.
"""
thickness_coordinates(coords) = @views [ coords[:,1] zero(coords[:,1]) coords[:,3] ]

function maximum_thickness_to_chord(coords)
    xs, thiccs    = @views coords[:,1], coords[:,3]
    max_thick_arg = argmax(thiccs)
    chord         = @views maximum(coords[:,1])
    @views xs[max_thick_arg] / chord, thiccs[max_thick_arg] / chord
end

"""
    read_foil(path :: String; header = true; name = "")

Generate a `Foil` from a file consisting of 2D coordinates with named arguments to skip the header (first line of the file) or assign a name.

By default, the header is assumed to exist and should contain the airfoil name, which is assigned to the name of the `Foil`.
"""
function read_foil(path :: String; name = "") 
    coords, foil_name = readdlm(path, header = true)
    if name != ""
        return Foil(Float16.(coords[:,1:2]), name)
    else
        return Foil(Float16.(coords[:,1:2]), strip(join(foil_name, " ")))
    end
end