plot_panel(panel :: Panel3D) = let ps = panel_coordinates(panel); Tuple.([ ps; [ps[1]] ]) end

"""
    plot_panels(panels :: Array{Panel3D})

Get the coordinates of an array of `Panel3D` for plotting.
"""
plot_panels(panels) = plot_panel.(vec(panels))

# foil_coords = [ [ [coord[1]; 0; coord[2]] .* chord .+ loc for coord in foil.coordinates ] for (chord, foil, loc) in zip(wing.right.chords[end:-1:1], wing.right.foils[end:-1:1], wing_coords) ]

function plot_planform(mesh :: Matrix{SVector{3,T}}) where T <: Real
    coords  =   [ 
                    mesh[1,1:end-1]; 
                    mesh[1:end-1,end]; 
                    mesh[end,end:-1:2]; 
                    mesh[end:-1:1,1] 
                ]
                    
    permutedims(combinedimsview(coords))
end

"""
    plot_planform(wing :: AbstractWing)

Get the planform coordinates of an `AbstractWing` for plotting.
"""
plot_planform(wing :: AbstractWing) = plot_planform(coordinates(wing))

plot_streamlines(system :: VortexLatticeSystem, points, length, num_steps) = Tuple.(streamlines(system, points, length, num_steps))

"""
    plot_surface(wing :: AbstractWing)
    plot_surface(wing :: AbstractWing, 
                 span_num :: Vector{Integer}, 
                 chord_num :: Integer)

Get the surface coordinates of an `AbstractWing` for plotting, optionally with specified spanwise and chordwise discretizations.
"""
plot_surface(wing :: AbstractWing) = plot_panels(mesh_wing(wing, [length(spans(wing))], minimum(length, foil.x for foil in foils(wing))))

plot_surface(wing :: AbstractWing, span_num, chord_num) = plot_panels(mesh_wing(wing, span_num, chord_num))

circle3D(r, n) = let arcs = 0:2π/n:2π; [ zeros(length(arcs)) r * cos.(arcs) r * sin.(arcs) ] end

# Fuselage
function plot_surface(fuse :: Fuselage, n_secs = 5)
    n_pts = 20
    xys = cosine_interpolation(fuse, n_secs)
    
    circs = combinedimsview(permutedims.(combinedimsview.([ eachrow(circ) .+ Ref([x; 0.; 0.] + fuse.position) for (x, circ) in zip(xys[:,1] .* fuse.length, circle3D.(xys[:,2], n_pts)) ])))

    plot_circs = [ reduce(vcat, splitdimsview(circs[n:(n+1),:,k:(k+1)])) for n in axes(circs, 1)[1:end-1] for k in axes(circs, 3)[1:end-1] ]

    return plot_circs
end