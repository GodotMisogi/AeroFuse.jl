## Wing analysis case
using AeroMDAO
import LinearAlgebra: norm

## Foils
naca_2412 = naca4((2,4,1,2))

##
foils = fill(naca_2412, 2)

##
@code_warntype maximum_thickness_to_chord.(foils)

##
@code_warntype upper_surface.(foils)

## Wing
wing =  Wing(foils     = foils,
             chords    = [1.0, 0.6],
             twists    = [2.0, 0.0],
             spans     = [4.0],
             dihedrals = [5.],
             sweeps      = [5.]);

##
x_w, y_w, z_w = wing_mac = mean_aerodynamic_center(wing)
S, b, c = projected_area(wing), span(wing), mean_aerodynamic_chord(wing);

## Chord processing
mean_chords = (AeroMDAO.forward_sum ∘ chords)(wing.right) / 2

## Wetted areas
S_wets  = @. mean_chords * wing.right.spans / cos(wing.right.dihedrals)

# Form factors
M       = 10. / 330.

##
@code_warntype AeroMDAO.camber_thickness(wing.left, 50)

##
@code_warntype AeroMDAO.max_thickness_to_chord_ratio_sweeps(wing.right, 60)