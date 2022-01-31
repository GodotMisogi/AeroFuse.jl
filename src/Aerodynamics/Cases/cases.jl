## 3D cases
#==========================================================================================#

"""
    solve_case(components :: Vector{Horseshoe}, fs :: Freestream, refs :: References;
               name = :aircraft :: Symbol, 
               print = false :: Boolean,
               print_components = false :: Boolean)

Perform a vortex lattice analysis given a vector of `Horseshoe`s, a `Freestream` condition, and `Reference` values.
"""
function solve_case(components, freestream :: Freestream, refs :: References; name = :aircraft, print = false, print_components = false, finite_core = false)
    system = solve_system(components, freestream, refs, finite_core)

    # Printing if needed
    if print_components
        print_coefficients(system, name, components = true)
    elseif print
        print_coefficients(system, name)
    end

    system
end

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
    chords     = vec(sum(average_chord, panels, dims = 1))

    # Compute spanwise coefficients
    CFs      = combinedimsview(CFs)
    span_CXs = @views spanwise_loading(CFs[1,:,:], area_scale)
    span_CYs = @views spanwise_loading(CFs[2,:,:], area_scale)
    span_CZs = @views spanwise_loading(CFs[3,:,:], area_scale)
    span_Cls = @views spanwise_loading(CFs[4,:,:], area_scale) ./ chords
    span_Cms = @views spanwise_loading(CFs[5,:,:], area_scale) ./ chords
    span_Cns = @views spanwise_loading(CFs[6,:,:], area_scale) ./ chords

    [ ys span_CXs span_CYs span_CZs span_Cls span_Cms span_Cns ]
end

spanwise_loading(wing :: WingMesh, CFs, S) = spanwise_loading(chord_panels(wing), CFs, S)

## Mesh connectivities
triangle_connectivities(inds) = @views [ 
                                         vec(inds[1:end-1,1:end-1]) vec(inds[1:end-1,2:end]) vec(inds[2:end,2:end])    ;
                                         vec(inds[2:end,2:end])     vec(inds[2:end,1:end-1]) vec(inds[1:end-1,1:end-1])
                                       ]

## Extrapolating surface values to neighbouring points
function extrapolate_point_mesh(mesh, weight = 0.75)
    m, n   = size(mesh)
    points = zeros(eltype(mesh), m + 1, n + 1)

    # The quantities are measured at the bound leg of a Horseshoe by default (0.25×)
    @views @. points[1:end-1,1:end-1] += weight       / 2 * mesh
    @views @. points[1:end-1,2:end]   += weight       / 2 * mesh
    @views @. points[2:end,1:end-1]   += (1 - weight) / 2 * mesh
    @views @. points[2:end,2:end]     += (1 - weight) / 2 * mesh

    points
end