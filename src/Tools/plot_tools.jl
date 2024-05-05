function plot_panel(panel::AbstractPanel3D)
    p = panel_coordinates(panel)
    combinedimsview([p; [p[1]]], (1))
end
"""
    plot_panels(panels :: Array{Panel3D})

Get the coordinates of an array of `Panel3D` for plotting.
"""
plot_panels(panels) = plot_panel.(panels)

# foil_coords = [ [ [coord[1]; 0; coord[2]] .* chord .+ loc for coord in foil.coordinates ] for (chord, foil, loc) in zip(wing.right.chords[end:-1:1], wing.right.foils[end:-1:1], wing_coords) ]

function plot_planform(mesh)
    coords = [
        mesh[1, 1:end-1];
        mesh[1:end-1, end];
        mesh[end, end:-1:2];
        mesh[end:-1:1, 1]
    ]

    combinedimsview(coords, (1))
end

"""
    plot_planform(wing :: AbstractWing)

Get the planform coordinates of an `AbstractWing` for plotting.
"""
plot_planform(wing::Wing) = plot_planform(coordinates(wing))

plot_streamlines(system::VortexLatticeSystem, points, length, num_steps) =
    combinedimsview(streamlines(system, points, length, num_steps), (1, 3))

"""
    plot_surface(wing :: AbstractWing)
    plot_surface(wing :: AbstractWing, 
                 span_num :: Vector{Integer}, 
                 chord_num :: Integer)

Get the surface coordinates of an `AbstractWing` for plotting, optionally with specified spanwise and chordwise discretizations.
"""
plot_surface(wing::AbstractWing) = plot_panels(
    mesh_wing(wing, [length(spans(wing))], minimum(length, foil.x for foil in foils(wing))),
)

plot_surface(wing::AbstractWing, span_num, chord_num) =
    plot_panels(mesh_wing(wing, span_num, chord_num))

circle3D(r, n) =
    let arcs = 0:2π/n:2π
        [zeros(length(arcs)) r * cos.(arcs) r * sin.(arcs)]
    end

# Fuselage
function plot_surface(fuse::Fuselage, n_secs = 5)
    n_pts = 20
    xys = cosine_interpolation(fuse, n_secs)

    circs = combinedimsview(
        permutedims.(
            combinedimsview.(
                [
                eachrow(circ) .+ Ref([x; 0.0; 0.0] + fuse.position)
                for (x, circ) in zip(
                    xys[:, 1] .* fuse.length,
                    circle3D.(xys[:, 2], n_pts),
                )
            ]
            )
        ),
    )

    # plot_circs = [ reduce(vcat, splitdimsview(circs[n:(n+1),:,k:(k+1)])) for n in axes(circs, 1)[1:end-1] for k in axes(circs, 3)[1:end-1] ]

    return circs
end

## Plots.jl recipes
#============================================#

@recipe function foil_plot(foil::Foil; camber = false, thickness = false)
    xlabel --> "x"
    ylabel --> "y"
    # aspect_ratio --> 1
    label --> foil.name

    if camber
        xcam = camber_line(foil)
        @series begin
            ls := :dash
            label := foil.name * " Camber"
            @views xcam[:, 1], xcam[:, 2]
        end
    end
    if thickness
        xthicc = thickness_line(foil)
        @series begin
            ls := :dot
            label := foil.name * " Thickness"
            @views xthicc[:, 1], xthicc[:, 2]
        end
    end

    foil.x, foil.y
end

@recipe function wing_plot(
    wing::AbstractWing,
    w = 0.25;
    mac = true,
    foils = true,
    n_span = 2,
    n_chord = 30,
)
    wing_plan = plot_planform(wing)

    # aspect_ratio --> true
    # zlim --> span(wing) .* (-0.5, 0.5)

    # Plotting planform coordinates
    @series begin
        @views wing_plan[:, 1], wing_plan[:, 2], wing_plan[:, 3]
    end

    # Plotting airfoil coordinates
    if foils
        foil_coords = combinedimsview(# Converting to 3D array
            surface_coordinates(wing,
                n_span * ones(
                    Int,
                    ifelse(wing.symmetry, 2length(wing.spans), length(wing.spans)),
                ), # Number of spanwise sections
                n_chord, # Number of chordwise points
                span_spacing = Uniform(), # Spanwise distribution
            ), (1, 3),
        )
        for i in axes(foil_coords, 3)
            @series begin
                # seriestype := :scatter
                primary := false
                @views foil_coords[:, 1, i], foil_coords[:, 2, i], foil_coords[:, 3, i]
            end
        end
    end

    # Plot mean aerodynamic center
    if mac
        wing_mac = mean_aerodynamic_center(wing, w)
        @series begin
            seriestype := :scatter
            primary := false
            @views [wing_mac[1]], [wing_mac[2]], [wing_mac[3]]
        end
    end
end

@recipe function mesh_plot(panels::Matrix{<:Panel3D})
    panels = plot_panels(panels)
    label --> ""
    for coords in panels
        @series begin
            seriestype := :path
            primary := false
            # linecolor := :lightgray
            # fillcolor := :lightgray
            @views coords[:, 1], coords[:, 2], coords[:, 3]
        end
    end
end

@recipe function wing_mesh_plot(wing::WingMesh, w = 0.25; mac = true)
    wing_plan = plot_planform(wing.surface)
    wing_pans = plot_panels(camber_panels(wing))
    # aspect_ratio --> true
    # zlim --> span(wing.surface) .* (-0.5, 0.5)

    for coords in wing_pans
        @series begin
            seriestype := :path
            primary := false
            linecolor := :lightgray
            # fillcolor := :lightgray
            @views coords[:, 1], coords[:, 2], coords[:, 3]
        end
    end

    if mac
        wing_mac = mean_aerodynamic_center(wing.surface, w)
        @series begin
            seriestype := :scatter
            primary := false
            # label --> "Mean Aerodynamic Chord"
            @views [wing_mac[1]], [wing_mac[2]], [wing_mac[3]]
        end
    end

    # @series begin
    #     # primary := false
    #     # linewidth := 2
    #     @views wing_plan[:,1], wing_plan[:,2], wing_plan[:,3]
    # end
end

@recipe function fuselage_plot(fuse::Fuselage, n = 20)
    fuse_pans = plot_surface(fuse, n)

    # for coords in fuse_pans
    #     @series begin
    #         seriestype := :path
    #         primary := false
    #         # linecolor := :lightgray
    #         # fillcolor := :lightgray
    #         @views coords[:,1], coords[:,2], coords[:,3]
    #     end
    # end
    for k in axes(fuse_pans, 3)[1:end-1]
        for n in axes(fuse_pans, 1)[1:end-1]
            @series begin
                seriestype := :path
                primary := false

                coo = fuse_pans[n:(n+1), :, (k+1):-1:k]
                coords = @views [coo[:, :, 1]; coo[:, :, 2]]
                @views coords[:, 1], coords[:, 2], coords[:, 3]
            end
        end
    end
end

@recipe function fuselage_plot(fuse::HyperEllipseFuselage; n_secs = 10, n_circ = 20)
    ts = LinRange(0, 1, n_secs)
    fuse_coo = coordinates(fuse, ts, n_circ) # Get coordinates
    fuse_ske = [fuse_coo[end÷4, :, 1:3]; fuse_coo[end÷2, :, 1:3]; fuse_coo[3*end÷4, :, 1:3]]
    fuse_plan = [reshape(fuse_coo, n_circ * n_secs * 3, 3); fuse_ske]
    @series begin
        # primary := false
        # linewidth := 2
        @views fuse_plan[:, 1], fuse_plan[:, 2], fuse_plan[:, 3]
    end
end

@recipe function streamline_plot(
    system::VortexLatticeSystem,
    seed,
    distance,
    num_stream_points = 100,
)
    streams = plot_streamlines(system, seed, distance, num_stream_points)

    for i in axes(streams, 3)
        @series begin
            seriestype := :line
            primary := false
            # linecolor := :green
            @views streams[:, 1, i], streams[:, 2, i], streams[:, 3, i]
        end
    end
end

@recipe function streamline_plot(
    system::VortexLatticeSystem,
    wing::AbstractWing;
    dist = 5 * mean_aerodynamic_chord(wing),
    num = 100,
    span = 20,
    linecolor = :green,
)
    # ys          = LinRange(-span(wing) / 2, span(wing) / 2, span_points)
    # init        = SVector.(0., ys, -0.5) 
    dx, dy, dz = 0, 0, 1e-3
    init = chop_leading_edge(wing, span)
    seed = init .+ Ref([dx, dy, dz])

    streams = plot_streamlines(system, seed, dist, num)

    for i in axes(streams, 3)
        @series begin
            seriestype := :line
            primary := false
            linecolor := linecolor
            @views streams[:, 1, i], streams[:, 2, i], streams[:, 3, i]
        end
    end
end

## Boxes
#===================================#

# Box type
struct Box{T <: Real}
    L::T
    w::T
    h::T
    position::SVector{3, T}
end

function Box(L, w, h, pos)
    T = promote_type(eltype(L), eltype(w), eltype(h), eltype(pos))

    return Box{T}(L, w, h, pos)
end

function coordinates(b::Box)
    L, w, h = b.L, b.w, b.h
    pts = [
        SVector(-L / 2, -w / 2, -h / 2),
        SVector(-L / 2, w / 2, -h / 2),
        SVector(-L / 2, -w / 2, h / 2),
        SVector(-L / 2, w / 2, h / 2),
        SVector(L / 2, -w / 2, h / 2),
        SVector(L / 2, -w / 2, -h / 2),
        SVector(L / 2, w / 2, -h / 2),
        SVector(L / 2, w / 2, h / 2),
    ]

    return Ref(b.position) .+ pts
end

# Plotting connections
const box_connections = [
    (1, 2, 3),
    (4, 2, 3),
    (4, 7, 8),
    (7, 5, 6),
    (2, 4, 7),
    (1, 6, 2),
    (2, 7, 6),
    (7, 8, 5),
    (4, 8, 5),
    (4, 5, 3),
    (1, 6, 3),
    (6, 3, 5),
];

# @recipe function mesh3d(b :: Box)
# pts = coordinates(b)

# @series begin
# Tuple.(pts)
# end
# end
