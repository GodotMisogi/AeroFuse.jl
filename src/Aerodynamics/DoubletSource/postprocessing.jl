function surface_velocities(prob :: DoubletSourceSystem)
    # Panel properties
    ps   = prob.surface_panels
    Δrs  = @views @. distance(ps[2:end], ps[1:end-1])
    αs   = @views panel_tangent.(ps[2:end])
    
    @views surface_velocities(prob.singularities[1:end-1], Δrs, αs, velocity(prob.freestream), false)
end

function surface_coefficients(prob :: DoubletSourceSystem)
    # Panel properties
    ps   = prob.surface_panels
    Δrs  = @views @. distance(ps[2:end], ps[1:end-1])
    xs   = @views combinedimsview(panel_points(ps)[2:end-1])[1,:]
    θs   = @views panel_angle.(ps[2:end])

    # Inviscid edge velocities
    u_es = @views surface_velocities(prob)

    # Aerodynamic coefficients
    cps  = @. pressure_coefficient(prob.freestream.magnitude, u_es)
    cls  = @. -cps * Δrs * cos(θs)
    cms  = @. -cls * xs * cos(prob.freestream.angle)

    cls, cms, cps
end

@views function surface_velocities(prob :: DoubletSourceSystem3D)
    # Due to influence
    vs = prob.velocity_influence_matrix * prob.singularities .+ Ref(prob.Umag * velocity(prob.freestream))
    npancd, npansp = size(prob.surface_panels)
    npanf = npancd * npansp
    vs = permutedims(reshape(vs, npansp, npancd))

    # Due to principal value
    φs = permutedims(reshape(prob.singularities[1:npanf], npansp, npancd))
    rs = collocation_point.(prob.surface_panels)
    ns = normal_vector.(prob.surface_panels)

    rso = [permutedims(rs[1,:]); rs[1:end-2,:]; permutedims(rs[end-1,:])]
    rno = [permutedims(rs[2,:]); rs[3:end,:]  ; permutedims(rs[end,:])  ]
    rea = [rs[:,1]               rs[:,1:end-2]  rs[:,end-1] ]
    rwe = [rs[:,2]               rs[:,3:end]    rs[:,end]   ]

    φso = [permutedims(φs[1,:]); φs[1:end-2,:]; permutedims(φs[end-1,:])]
    φno = [permutedims(φs[2,:]); φs[3:end,:]  ; permutedims(φs[end,:])  ]
    φea = [φs[:,1]               φs[:,1:end-2]  φs[:,end-1] ]
    φwe = [φs[:,2]               φs[:,3:end]    φs[:,end]   ]

    ∇Γ = ((rwe .- rno) .* (φwe + φno) + (rso .- rwe) .* (φso + φwe) + (rea .- rso) .* (φea + φso) + (rno .- rea) .* (φno + φea)) .× ns ./ norm.((rno-rso) .× (rwe-rea))
    ∇Γ[1,:] ./= 2
    ∇Γ[end,:] ./= 2
    vs -= ∇Γ/2

    return vs
end

function surface_coefficients(prob :: DoubletSourceSystem3D, wing :: Wing)
    # Panel properties
    ps = prob.surface_panels
    # rs = midpoint(ps)
    ns = normal_vector.(ps)
    As = panel_area.(ps) .* ns
    A = projected_area(wing)

    # Inviscid edge velocities
    vs = surface_velocities(prob)

    # Aerodynamic coefficients
    CP = pressure_coefficient.(prob.Umag, vs)
    CF = -sum(CP .* As) / A
    # Cm  = -sum(CP .* (rq .× As)) / A / cbar

    cos, _, sin = velocity(prob.freestream)
    CL = CF[3] * cos - CF[1] * sin
    CD = CF[3] * sin + CF[1] * cos

    CL, CD, CP
end

"""
    plot_scalar_field(ps, scalarfield)

Return necessary elements for plotting a scalar field over the aircraft/wing surface.
"""
@views function plot_scalar_field(ps, scalarfield)
    X = zeros(3, length(ps) * 2)
    for i ∈ eachindex(ps)
        p = 2 * i
        X[:,p-1] .= xs(ps[i])[[1,2,3]]
        X[:, p ] .= xs(ps[i])[[3,4,1]]
    end

    Y = zeros(3, length(ps) * 2)
    for i ∈ eachindex(ps)
        p = 2 * i
        Y[:,p-1] .= ys(ps[i])[[1,2,3]]
        Y[:, p ] .= ys(ps[i])[[3,4,1]]
    end

    Z = zeros(3, length(ps) * 2)
    for i ∈ eachindex(ps)
        p = 2 * i
        Z[:,p-1] .= zs(ps[i])[[1,2,3]]
        Z[:, p ] .= zs(ps[i])[[3,4,1]]
    end

    C = zeros(3, length(scalarfield) * 2)
    for i ∈ eachindex(ps)
        p = 2 * i
        C[:,p-1] .= scalarfield[i]
        C[:, p ] .= scalarfield[i]
    end
    xyz = reshape([X[:] Y[:] Z[:]]', :)
    pts = connect(xyz, Point{3})
    cnt = connect(1:length(X), TriangleFace)
    
    return pts, cnt, C
end

"""
    offbody_velocity(pt :: Point3D, prob :: DoubletSourceSystem3D)

Calculate offbody velocity at a given point `pt` due to panel influence of `prob`
"""
function offbody_velocity(pt, prob)
    ps = prob.surface_panels
    φs = prob.singularities

    vel = zeros(3)
    for i ∈ eachindex(ps)
        vel += quadrilateral_doublet_velocity(ps[i], pt) * φs[i]
    end

    return vel
end


lift_coefficient(prob :: DoubletSourceSystem) = 2 * last(prob.singularities) / prob.freestream.magnitude
