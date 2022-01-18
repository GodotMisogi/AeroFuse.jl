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

## Span-loading
function span_loads(panels, CFs_wind, S)
    CFs  = combinedimsview(CFs_wind)
    CDis = @views CFs[1,:,:]
    CYs  = @views CFs[2,:,:]
    CLs  = @views CFs[3,:,:]

    area_scale = S ./ vec(sum(panel_area, panels, dims = 1))
    span_CDis  = vec(sum(CDis, dims = 1)) .* area_scale
    span_CYs   = vec(sum(CYs,  dims = 1)) .* area_scale
    span_CLs   = vec(sum(CLs,  dims = 1)) .* area_scale

    ys = @views vec(mean(x -> midpoint(x)[2], panels, dims = 1))

    [ ys span_CDis span_CYs span_CLs ]
end

span_loads(wing :: WingMesh, CFs_wind, S) = span_loads(chord_panels(wing), CFs_wind, S)

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
    @views points[1:end-1,1:end-1] += weight       / 2 * mesh
    @views points[1:end-1,2:end]   += weight       / 2 * mesh
    @views points[2:end,1:end-1]   += (1 - weight) / 2 * mesh
    @views points[2:end,2:end]     += (1 - weight) / 2 * mesh

    points
end