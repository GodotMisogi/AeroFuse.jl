using CoordinateTransformations, Rotations, StaticArrays
using AeroMDAO


## Define type
struct WingCoordinates{T <: Real}
    wing        :: Union{HalfWing{T}, Wing{T}}
    position    :: SVector{3,T}
    orientation :: AngleAxis{T}
end

WingCoordinates(wing; position = zeros(3), angle :: T = 0., axis = [1.,0.,0.]) where T <: Real = WingCoordinates(wing, position, AngleAxis{T}(angle, axis...))

## Methods
component(comp :: WingCoordinates)   = comp.wing
position(comp :: WingCoordinates)    = comp.position
orientation(comp :: WingCoordinates) = comp.orientation

affine_transformation(wing :: WingCoordinates) = Translation(position(wing)) ∘ LinearMap(orientation(wing))

## Tests
wing = Wing(foils     = Foil.(fill(naca4((0,0,1,2)), 2)),
            chords    = [0.7, 0.42],
            twists    = [0.0, 0.0],
            spans     = [1.25],
            dihedrals = [0.],
            sweep_LEs = [6.39])

wing_vlm = WingCoordinates(wing,
                           position = SVector(4., 0, 0.),
                           angle    = -1.,
                           axis     = SVector(0., 1., 0.))

## Mean aerodynamic center
AeroMDAO.AircraftGeometry.mean_aerodynamic_center(wing :: WingCoordinates) = affine_transformation(wing)(AeroMDAO.mean_aerodynamic_center(component(wing)))

mean_aerodynamic_center(wing_vlm)

## Chordline coordinates
chord_coordinates(wing :: WingCoordinates, span_num, chord_num; spacings = ["sine"]) = affine_transformation(wing).(AeroMDAO.chord_coordinates(component(wing), span_num, chord_num; spacings = spacings))

chord_coordinates(wing_vlm, [6], 6)

## Horseshoe mesh
AeroMDAO.AircraftGeometry.mesh_horseshoes(wing :: WingCoordinates, span_num, chord_num; spacings = ["sine"]) = make_panels(affine_transformation(wing).(chord_coordinates(wing, span_num, chord_num; spacings = spacings)))

mesh_horseshoes(wing_vlm, [6], 6)

## Camberline coordinates
AeroMDAO.AircraftGeometry.camber_coordinates(wing :: WingCoordinates, span_num, chord_num; spacings = ["sine"]) = affine_transformation(wing).(AeroMDAO.chord_coordinates(component(wing), span_num, chord_num; spacings = spacings))

camber_coordinates(wing_vlm, [6], 6)

## Camber mesh
AeroMDAO.AircraftGeometry.mesh_cambers(wing :: WingCoordinates, span_num, chord_num; spacings = ["sine"]) = make_panels((camber_coordinates(wing, span_num, chord_num; spacings = spacings)))

mesh_cambers(wing_vlm, [6], 6)