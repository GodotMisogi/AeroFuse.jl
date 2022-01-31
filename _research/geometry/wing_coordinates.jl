using CoordinateTransformations, Rotations, StaticArrays
using AeroMDAO

# An internal definition of coordinates in structs is restrictive and not extensible.
# Consider:
# aff, g :: [Loc] -> [Loc]
# f1, f2 :: HW -> [Loc]
# such that f1 = g ∘ aff ∘ f2, with aff only being applied when required.
#
# Need to think about this monadically. A monad m is defined as a monoid in a category of endofunctors in terms of the monoidal maps join and return in Haskell:
# join :: m (m a) -> m a
# return :: a -> m a
# This corresponds to composition of morphisms in the Kleisli category, implemented as the bind operator in Haskell.
# >>= :: m a -> (a -> m b) -> m b
# xs >>= f = join (fmap f xs)


## Define type
struct WingCoordinates{T <: Real} <: AbstractAircraft
    wing        :: Union{HalfWing{T}, Wing{T}}
    position    :: SVector{3,T}
    orientation :: AngleAxis{T}
end

WingCoordinates(wing; position = zeros(3), angle :: T = 0., axis = SVector(1., 0., 0.)) where T <: Real = WingCoordinates(wing, position, AngleAxis{T}(angle, axis...))

## Methods
component(comp :: WingCoordinates)   = comp.wing
position(comp :: WingCoordinates)    = comp.position
orientation(comp :: WingCoordinates) = comp.orientation

affine_transformation(wing :: WingCoordinates) = Translation(position(wing)) ∘ LinearMap(orientation(wing))

## Tests
wing = Wing(foils     = fill(naca4((0,0,1,2)), 2),
            chords    = [0.7, 0.42],
            twists    = [0.0, 0.0],
            spans     = [1.25],
            dihedrals = [0.],
            LE_sweeps = [6.39])

wing_vlm = WingCoordinates(wing,
                           position = SVector(4., 0, 0.),
                           angle    = -1.,
                           axis     = SVector(0., 1., 0.))

## Mean aerodynamic center
AeroMDAO.AircraftGeometry.mean_aerodynamic_center(wing :: WingCoordinates) = affine_transformation(wing)(AeroMDAO.mean_aerodynamic_center(component(wing)))

mean_aerodynamic_center(wing_vlm)

## Chordline coordinates
chord_coordinates(wing :: WingCoordinates, span_num, chord_num; spacings = [Cosine()]) = affine_transformation(wing).(AeroMDAO.chord_coordinates(component(wing), span_num, chord_num; spacings = spacings))

chord_coordinates(wing_vlm, [6], 6)

## Horseshoe mesh
AeroMDAO.AircraftGeometry.mesh_horseshoes(wing :: WingCoordinates, span_num, chord_num; spacings = [Cosine()]) = make_panels(affine_transformation(wing).(chord_coordinates(wing, span_num, chord_num; spacings = spacings)))

mesh_horseshoes(wing_vlm, [6], 6)

## Camberline coordinates
AeroMDAO.AircraftGeometry.camber_coordinates(wing :: WingCoordinates, span_num, chord_num; spacings = [Cosine()]) = affine_transformation(wing).(AeroMDAO.chord_coordinates(component(wing), span_num, chord_num; spacings = spacings))

camber_coordinates(wing_vlm, [6], 6)

## Camber mesh
AeroMDAO.AircraftGeometry.mesh_cambers(wing :: WingCoordinates, span_num, chord_num; spacings = [Cosine()]) = make_panels((camber_coordinates(wing, span_num, chord_num; spacings = spacings)))

mesh_cambers(wing_vlm, [6], 6)