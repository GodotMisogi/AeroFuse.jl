##
using AeroMDAO
using BenchmarkTools

## Aircraft definitions
#==========================================================================================#

# Wing
wing = Wing(foils     = fill(naca4(2,4,1,2), 2),
            chords    = [1.0, 0.6],
            twists    = [2.0, 0.0],
            spans     = [4.0],
            dihedrals = [5.],
            sweeps      = [5.]);

# Horizontal tail
htail = Wing(foils     = fill(naca4(0,0,1,2), 2),
             chords    = [0.7, 0.42],
             twists    = [0.0, 0.0],
             spans     = [1.25],
             dihedrals = [0.],
             sweeps      = [6.39],
             position  = [4., 0, 0],
             angle     = -2.,
             axis      = [0., 1., 0.])

# Vertical tail
vtail = HalfWing(foils     = fill(naca4(0,0,0,9), 2),
                 chords    = [0.7, 0.42],
                 twists    = [0.0, 0.0],
                 spans     = [1.0],
                 dihedrals = [0.],
                 sweeps      = [7.97],
                 position  = [4., 0, 0],
                 angle     = 90.,
                 axis      = [1., 0., 0.])

## WingMesh type
wing_mesh  = WingMesh(wing, [40], 20, 
                      span_spacing = Cosine()
                     )
htail_mesh = WingMesh(htail, [20], 10, 
                      span_spacing = Cosine()
                     )
vtail_mesh = WingMesh(vtail, [20], 10, 
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

## Sequential
println("AeroMDAO Aircraft Sequential -")
@benchmark systems = map(α -> solve_case($aircraft, Freestream(alpha = α), $refs), $αs)

## Thread-based parallelism
using Folds

println("AeroMDAO Aircraft Thread-Parallel -")
@benchmark systems = Folds.map(α -> solve_case($aircraft, Freestream(alpha = α), $refs), $αs)

## Distributed-computing parallelism
using Distributed
# addprocs(3)
println("AeroMDAO Aircraft Proc-Parallel -")
@benchmark systems = pmap(α -> solve_case($(aircraft[:]), Freestream(alpha = α), $refs), $αs)