## 3D cases
#==========================================================================================#

"""
    solve_case(
        components :: Vector{Horseshoe}, 
        fs :: Freestream, 
        refs :: References;
        name = :aircraft :: Symbol, 
        print = false :: Boolean,
        print_components = false :: Boolean
    )

Perform a vortex lattice analysis given a vector of `Horseshoe`s, a `Freestream` condition, and `References` values.
"""
function solve_case(components :: DenseArray, freestream :: Freestream, refs :: References; name = :aircraft, print = false, print_components = false)
    system = VortexLatticeSystem(components, freestream, refs)

    # Printing if needed
    if print_components
        print_coefficients(system, name, components = true)
    elseif print
        print_coefficients(system, name)
    end

    system
end

solve_case(meshes, freestream :: Freestream, refs :: References; name = :aircraft, print = false, print_components = false) = solve_case(ComponentVector(meshes), freestream, refs; name = name, print = print, print_components = print_components)

# Method forwarding
#
# macro forward(ex, fs)
#     @capture(ex, T_.field_) || error("Syntax: @forward T.x f, g, h")
#     T = esc(T)
#     fs = isexpr(fs, :tuple) ? map(esc, fs.args) : [esc(fs)]
#     :($([:($f(x::$T, args...; kwargs...) =
#            (Base.@_inline_meta; $f(x.$field, args...; kwargs...)))
#          for f in fs]...);
#       nothing)
# end

# struct MeshVortexLatticeSystem{M <: WingMesh, VLM <: AbstractVortexLatticeSystem} <: AbstractVortexLatticeSystem
#     meshes :: Vector{M}
#     system :: VLM 
# end

# @forward MeshVortexLatticeSystem.system surface_coefficients, surface_forces, surface_forces, surface_moments, surface_velocities, surface_dynamics, nearfield_coefficients, farfield_coefficients, farfield_forces, nearfield, farfield

## Placeholder for functions I'm not sure where to put
#==========================================================================================#

spanwise_loading(CXs, weighted_areas) = vec(sum(CXs, dims = 1)) .* weighted_areas

## Span-loading
"""
    spanwise_loading(panels, CFs, S)
    spanwise_loading(mesh :: WingMesh, CFs, S)

Compute the spanwise loading of the forces given panels, associated force coefficients, and the reference area.

Alternatively, provide a `WingMesh` type which gets the panels for you.
"""
function spanwise_loading(panels, CFs, S)
    # Get y-coordinates of spanwise strips
    ys = @views vec(mean(x -> midpoint(x)[2], panels, dims = 1))

    # Compute weighted areas for spanwise strips
    area_scale = S ./ vec(sum(panel_area, panels, dims = 1))
    # chords     = vec(sum(average_chord, panels, dims = 1))

    # Compute spanwise coefficients
    span_CFs = permutedims(combinedimsview(spanwise_loading(CFs, area_scale)))

    [ ys span_CFs ]
end

spanwise_loading(wing :: WingMesh, CFs, S) = spanwise_loading(chord_panels(wing), CFs, S)

"""
    spanwise_loading(wing_mesh :: WingMesh, ref :: References, CFs, Γ)

Obtain the spanwise aerodynamic loads `(CDi, CY, CL, CL_norm)` for a given `WingMesh`, reference values `ref`, surface coefficients `CFs`, and circulations ``Γ``. 
    
The chordwise normalized lift coefficient loading `CL_norm` is calculated via the formula ``C_L = 2Γ / ρVc`` based on the reference speed ``V`` and chord ``c`` from `ref`.

Note that the surface coefficients `CFs` should be computed in wind axes to match `CL_norm`.
"""
function spanwise_loading(wing_mesh :: WingMesh, ref :: References, CFs, Γs)
    span_loads = spanwise_loading(wing_mesh, CFs, ref.area) # 
    CL_loads = vec(sum(Γs, dims = 1)) / (0.5 * ref.speed * ref.chord) # Normalized CL loading

    [ span_loads CL_loads ]
end

## Mesh connectivities
@views triangle_connectivities(inds) = [ 
        vec(inds[1:end-1,1:end-1]) vec(inds[1:end-1,2:end]) vec(inds[2:end,2:end])    ;
        vec(inds[2:end,2:end])     vec(inds[2:end,1:end-1]) vec(inds[1:end-1,1:end-1])
    ]

## Extrapolating surface values to neighbouring points
function extrapolate_point_mesh(mesh, weight = 0.75)
    m, n   = size(mesh)
    points = zeros(eltype(mesh), m + 1, n + 1)

    # The quantities are measured at the bound leg of a Horseshoe by default (0.25×)
    @. points[1:end-1,1:end-1] += weight       / 2 * mesh
    @. points[1:end-1,2:end]   += weight       / 2 * mesh
    @. points[2:end,1:end-1]   += (1 - weight) / 2 * mesh
    @. points[2:end,2:end]     += (1 - weight) / 2 * mesh

    points
end