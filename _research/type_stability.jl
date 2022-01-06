## Wing analysis case
using AeroMDAO
import LinearAlgebra: norm

## Surfaces



## Foils
naca_2412 = Foil(naca4(2,4,1,2))

##
foil = fill(naca_2412, 2)

## Wing
@code_warntype Wing(foils     = foil,
            chords    = [1.0, 0.6],
            twists    = [2.0, 0.0],
            spans     = [4.0],
            dihedrals = [5.],
            LE_sweeps = [5.]);

##
x_w, y_w, z_w = wing_mac = mean_aerodynamic_center(wing)
S, b, c = projected_area(wing), span(wing), mean_aerodynamic_chord(wing);

##
# Chord processing
mean_chords = (fwdsum ∘ chords)(wing.right) / 2

# Wetted areas
S_wets  = @. mean_chords * wing.right.spans / cos(wing.right.dihedrals)

# Form factors
M       = 10. / 330.
@code_warntype AeroMDAO.max_thickness_to_chord_ratio_sweeps(wing.right, 60)