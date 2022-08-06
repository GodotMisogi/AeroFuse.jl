## 3D Panels
#==========================================================================================#
abstract type AbstractPanel3D <: AbstractPanel end

"""
    Panel3D(p1, p2, p3, p4)

Four Cartesian coordinates `p1, p2, p3, p4` representing corners of a panel in 3 dimensions. The following commutative diagram (math joke) depicts the order:

```
z → y
↓
x
        p1 —→— p4
        |       |
        ↓       ↓
        |       |
        p2 —→— p3
```
"""
struct Panel3D{T <: Real} <: AbstractPanel3D
    p1 :: Point3D{T}
    p2 :: Point3D{T}
    p3 :: Point3D{T}
    p4 :: Point3D{T}
end

struct WakePanel3D{T <: Real} <: AbstractPanel3D
    p1 :: Point3D{T}
    p2 :: Point3D{T}
    p3 :: Point3D{T}
    p4 :: Point3D{T}
end

# function Panel3D(p1 :: Point3D{T}, p2 :: Point3D{T}, p3 :: Point3D{T}, p4 :: Point3D{T}) where T <: Real
# 	Panel3D{T}(p1, p2, p3, p4)
# end

# function WakePanel3D(p1 :: Point3D{T}, p2 :: Point3D{T}, p3 :: Point3D{T}, p4 :: Point3D{T}) where T <: Real
# 	WakePanel3D{T}(p1, p2, p3, p4)
# end


Panel3D(p1, p2, p3, p4) = let T = promote_type(eltype(p1), eltype(p2), eltype(p3), eltype(p4)); Panel3D{T}(p1, p2, p3, p4) end

Panel3D((p1, p2, p3, p4)) = Panel3D(p1, p2, p3, p4)

collocation_point(panel :: AbstractPanel3D, a = 0.5) = (p1(panel) + p2(panel) + p3(panel) + p4(panel)) / 4

## Some algebraic operations on Panel3D
+(p :: AbstractPanel3D, v) = Panel3D(p.p1 + v, p.p2 + v, p.p3 + v, p.p4 + v)
-(p :: AbstractPanel3D, v) = Panel3D(p.p1 - v, p.p2 - v, p.p3 - v, p.p4 - v)
×(p :: AbstractPanel3D, v) = Panel3D(p.p1 × v, p.p2 × v, p.p3 × v, p.p4 × v)
(trans::AffineMap)(p::AbstractPanel3D) = Panel3D(trans(p.p1), trans(p.p2), trans(p.p3), trans(p.p4))

Base.length(:: Panel3D) = 1

average_chord(panel :: AbstractPanel3D) = (p2(panel) - p1(panel) + p3(panel) - p4(panel)) / 2
average_width(panel :: AbstractPanel3D) = (p4(panel) - p1(panel) + p3(panel) - p2(panel)) / 2


"""
    panel_coordinates(panel :: Panel3D)

Compute the coordinates of a `Panel3D`.
"""
panel_coordinates(panel :: Panel3D) = SVector( p1(panel), p2(panel), p3(panel), p4(panel) )


"""
    make_panels(xyzs)

Convert an array of coordinates corresponding to a wing, ordered from root to tip and leading-edge to trailing-edge, into panels.
"""
make_panels(xyzs :: AbstractArray{<: SVector{3,<: Real}}) = @views Panel3D.(xyzs[1:end-1,1:end-1], xyzs[2:end,1:end-1], xyzs[2:end,2:end], xyzs[1:end-1,2:end])


"""
    midpoint(panel :: AbstractPanel3D)

Compute the midpoint of an `AbstractPanel3D`.
"""
midpoint(panel :: AbstractPanel3D) = (p1(panel) + p2(panel) + p3(panel) + p4(panel)) / 4


"""
    panel_normal(panel :: AbstractPanel3D)

Compute the unit normal vector of an `AbstractPanel3D` normalisation.
"""
panel_normal(panel :: AbstractPanel3D) = normalize(×(p4(panel) - p2(panel), p3(panel) - p1(panel)))


"""
    transform_normal(panel :: Panel3D, h_l, g_l)

Transform the normal vector ``n̂₀`` of a `Panel3D` about the hinge axis ``ĥₗ`` by the control gain ``gₗ``.

The transformation is the following: ``n̂ₗ = gₗ ĥₗ × n̂₀`` (FVA, Drela, 6.36).
"""
transform_normal(panel :: Panel3D, h_l, g_l) = g_l * cross(h_l, panel_normal(panel))


"""
    panel_area(panel :: AbstractPanel3D)

Compute the area of a planar quadrilateral 3D panel.
"""
function panel_area(panel :: AbstractPanel3D)
    coord = panel_coordinates(panel)

	d_ij(i :: Int64, j :: Int64, local_coordinates) = norm(local_coordinates[j] - local_coordinates[i])
    d12 = d_ij(1, 2, coord)
    d23 = d_ij(2, 3, coord)
    d34 = d_ij(3, 4, coord)
    d41 = d_ij(4, 1, coord)
    d13 = d_ij(1, 3, coord)

    s1 = 0.5 * (d12 + d23 + d13)
    s2 = 0.5 * (d34 + d41 + d13)

    A1 = sqrt(s1 * (s1 - d12) * (s1 - d23) * (s1 - d13))
    A2 = sqrt(s2 * (s2 - d34) * (s2 - d41) * (s2 - d13))

    return A1 + A2
end


"""
    wetted_area(panels :: Array{Panel3D})

Compute the total wetted area by summing the areas of an array of `Panel3D`.
"""
wetted_area(panels) = sum(panel -> panel_area(panel), panels)


"""
    reflect_xz(panel :: Panel3D)

Reflect a Panel3D with respect to the ``x``-``z`` plane of its reference coordinate system.
"""
reflect_xz(panel :: Panel3D) = Panel3D((reflect_xz ∘ p1)(panel), (reflect_xz ∘ p2)(panel), (reflect_xz ∘ p3)(panel), (reflect_xz ∘ p4)(panel))

# Compute local axis coordinates
function local_coordinate_system(stream, normie) 
    s_hat = normalize(stream)
    n_hat = normalize(normie)
    c_hat = s_hat × n_hat

    reduce(hcat, @SMatrix [ c_hat; s_hat; n_hat ])
end

# Compute local axis coordinates
local_coordinate_system(panel :: Panel3D) = local_coordinate_system((panel.p4 - panel.p1 + panel.p3 - panel.p2) / 2, panel_normal(panel))

"""
    transform_panel(panel :: AbstractPanel3D, point :: Point3D) -> panel :: AbstractPanel3D, point :: Point3D

Transform point and panel from GCS into panel LCS.
"""
function transform_panel(panel :: AbstractPanel3D, point :: Point3D)
    # Translate p1 to xy plane
    tr1 = Translation(0., 0., -panel.p1.z)
    local_panel = tr1(panel)
    local_point = tr1(point)

    # Find normal of panel
    n = panel_normal(panel)

    # Find rotation axis
    d = n × SVector(0.,0.,1.)

    # No rotation if already aligned
    if norm(d) <= 1e-7
        return local_panel, local_point
    end

    # Find rotation angle
    θ = acos(n.z / norm(n))
    tr2 = recenter(LinearMap(AngleAxis(θ, d.x, d.y, d.z)), local_panel.p1)

    local_panel = tr2(local_panel)
    local_point = tr2(local_point)

    # tr = tr1 ∘ tr2
    # invtr = inv(tr)

    return local_panel, local_point
end


"""
	get_transformation(panel :: AbstractPanel3D) -> tr :: AffineMap

Calculate required transformation from GCS to panel LCS.
"""
function get_transformation(panel :: AbstractPanel3D)
    # Translate p1 to xy plane
	p1 = panel.p1
    tr1 = Translation(0., 0., -p1.z)

    # Find normal of panel
    n = panel_normal(panel)

    # Find rotation axis
    d = panel_normal(panel) × SVector(0.,0.,1.)

    # No rotation if already aligned
    if norm(d) <= 1e-7
        return local_panel, local_point
    end

    # Find rotation angle
    θ = acos(n.z / norm(n))
    tr2 = recenter(LinearMap(AngleAxis(θ, d.x, d.y, d.z)), tr1(p1))

    return tr1 ∘ tr2
end


"""
wake_panel(panels :: AbstractPanel3D, bound, α)

Calculate required transformation from GCS to panel LCS.
"""
function wake_panel(panels :: AbstractArray{<: AbstractPanel3D}, bound, α, β)
    pt1 = 0.5 * ( p1(first(panels)) + p2(last(panels)) )
	pt4 = 0.5 * ( p4(first(panels)) + p3(last(panels)) )
	dx, dy, dz = velocity(Freestream(α, β, zeros(3))) * bound
	pt2 = pt1 + Point3D(dx, dy, dz)
	pt3 = pt4 + Point3D(dx, dy, dz)

	return WakePanel3D(pt1, pt2, pt3, pt4)
end