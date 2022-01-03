##
using AeroMDAO
using BenchmarkTools

## Aircraft definitions
#==========================================================================================#

# Wing
wing = Wing(foils     = Foil.(fill(naca4(2,4,1,2), 2)),
            chords    = [1.0, 0.6],
            twists    = [2.0, 0.0],
            spans     = [4.0],
            dihedrals = [5.],
            LE_sweeps = [5.]);

x_w, y_w, z_w = wing_mac = mean_aerodynamic_center(wing)
S, b, c = projected_area(wing), span(wing), mean_aerodynamic_chord(wing);

# Horizontal tail
htail = Wing(foils     = Foil.(fill(naca4(0,0,1,2), 2)),
             chords    = [0.7, 0.42],
             twists    = [0.0, 0.0],
             spans     = [1.25],
             dihedrals = [0.],
             LE_sweeps = [6.39],
             position  = [4., 0, 0],
             angle     = -2.,
             axis      = [0., 1., 0.])

# Vertical tail
vtail = HalfWing(foils     = Foil.(fill(naca4(0,0,0,9), 2)),
                 chords    = [0.7, 0.42],
                 twists    = [0.0, 0.0],
                 spans     = [1.0],
                 dihedrals = [0.],
                 LE_sweeps = [7.97],
                 position  = [4., 0, 0],
                 angle     = 90.,
                 axis      = [1., 0., 0.])

## WingMesh type
wing_mesh  = WingMesh(wing, [10], 10, 
                      span_spacing = Cosine()
                     )
htail_mesh = WingMesh(htail, [6], 6, 
                      span_spacing = Cosine()
                     )
vtail_mesh = WingMesh(vtail, [6], 5, 
                      span_spacing = Cosine()
                     )

aircraft = ComponentArray(
                          wing  = make_horseshoes(wing_mesh),
                          htail = make_horseshoes(htail_mesh),
                          vtail = make_horseshoes(vtail_mesh)
                         );

ρ       = 1.225
x_ref   = [0.5, 0., 0.]
S, b, c = 9.0, 10.0, 0.9
αs      = -5:5

refs    = References(
                        speed    = 1.0, 
                        density  = ρ,
                        area     = S,
                        span     = b,
                        chord    = c,
                        location = x_ref,
                    )

## Angle of attack sweep
#==========================================================================================#

using Distributed

## Sequential
println("AeroMDAO Aircraft Sequential -")
@benchmark systems = map(α -> solve_case($aircraft, Freestream(alpha = α), $refs), αs)

## Parallel
println("AeroMDAO Aircraft Parallel -")
@benchmark systems = pmap(α -> solve_case($aircraft, Freestream(alpha = α), $refs), αs)