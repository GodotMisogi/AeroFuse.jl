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


Panel3D(p1, p2, p3, p4) = let T = promote_type(eltype(p1), eltype(p2), eltype(p3), eltype(p4)); Panel3D{T}(p1, p2, p3, p4) end

WakePanel3D(p1, p2, p3, p4) = let T = promote_type(eltype(p1), eltype(p2), eltype(p3), eltype(p4)); WakePanel3D{T}(p1, p2, p3, p4) end

Panel3D((p1, p2, p3, p4)) = Panel3D(p1, p2, p3, p4)

collocation_point(panel :: AbstractPanel3D, a = 0.5) = (p1(panel) + p2(panel) + p3(panel) + p4(panel)) / 4

## Some algebraic operations on Panel3D
+(p :: AbstractPanel3D, v) = Panel3D(p.p1 + v, p.p2 + v, p.p3 + v, p.p4 + v)
-(p :: AbstractPanel3D, v) = Panel3D(p.p1 - v, p.p2 - v, p.p3 - v, p.p4 - v)
×(p :: AbstractPanel3D, v) = Panel3D(p.p1 × v, p.p2 × v, p.p3 × v, p.p4 × v)

(T :: AffineMap)(p :: AbstractPanel3D) = typeof(p)(T(p.p1), T(p.p2), T(p.p3), T(p.p4))
(T :: LinearMap)(p :: AbstractPanel3D) = typeof(p)(T(p.p1), T(p.p2), T(p.p3), T(p.p4))

Base.length(:: Panel3D) = 1

# Transformation laws with vectors
+(p, v) = Panel3D(p.p1 + v, p.p2 + v, p.p3 + v, p.p4 + v)
-(p, v) = Panel3D(p.p1 - v, p.p2 - v, p.p3 - v, p.p4 - v)
×(p, v) = Panel3D(p.p1 × v, p.p2 × v, p.p3 × v, p.p4 × v)

average_chord(panel :: Panel3D) = (p2(panel) - p1(panel) + p3(panel) - p4(panel)) / 2
average_width(panel :: Panel3D) = (p4(panel) - p1(panel) + p3(panel) - p2(panel)) / 2

"""
    panel_coordinates(panel :: Panel3D)

Compute the coordinates of a `Panel3D`.
"""
panel_coordinates(panel :: AbstractPanel3D) = [ p1(panel), p2(panel), p3(panel), p4(panel) ]

"""
    make_panels(xyzs)

Convert an array of coordinates corresponding to a wing, ordered from root to tip and leading-edge to trailing-edge, into panels.
"""
make_panels(xyzs) = @views Panel3D.(xyzs[1:end-1,1:end-1], xyzs[2:end,1:end-1], xyzs[2:end,2:end], xyzs[1:end-1,2:end])

"""
    transform(panel :: Panel3D, rotation, translation)

Perform an affine transformation on the coordinates of a `Panel3D` given a rotation matrix and translation vector.
"""
function transform(panel :: Panel3D, rotation, translation) 
    T = LinearMap(rotation) ∘ Translation(translation)
    Panel3D(T(panel.p1), T(panel.p2), T(panel.p3), T(panel.p4))
end

"""
    midpoint(panel :: AbstractPanel3D)

Compute the midpoint of an `AbstractPanel3D`.
"""
midpoint(panel :: AbstractPanel3D) = (p1(panel) + p2(panel) + p3(panel) + p4(panel)) / 4

"""
    normal_vector(panel :: Panel3D)

Compute the normal vector of an `AbstractPanel3D`.
"""
normal_vector(panel :: AbstractPanel3D) = let p31 = panel.p3 - panel.p1, p42 = panel.p4 - panel.p2; normalize(p31 × p42) end

"""
    transform_normal(panel :: Panel3D, h_l, g_l)

Transform the normal vector ``n̂₀`` of a `Panel3D` about the hinge axis ``ĥₗ`` by the control gain ``gₗ``.

The transformation is the following: ``n̂ₗ = gₗ ĥₗ × n̂₀`` (Flight Vehicle Aerodynamics, M. Drela, Eq. 6.36).
"""
transform_normal(panel :: Panel3D, h_l, g_l) = g_l * cross(h_l, normal_vector(panel))

"""
    panel_area(panel :: AbstractPanel3D)

Compute the area of a planar quadrilateral 3D panel.
"""
panel_area(panel :: Panel3D) = let p31 = panel.p3 - panel.p1, p42 = panel.p4 - panel.p2; 1/2 * norm(p31 × p42) end

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

# Determine panel's local axis system
function local_coordinate_system(stream, normie) 
    s_hat = normalize(stream)
    n_hat = normalize(normie)
    l_hat = n_hat × s_hat

    return [ s_hat l_hat n_hat ]
end

# Compute local axis coordinates
local_coordinate_system(panel :: AbstractPanel3D) = local_coordinate_system((panel.p2 + panel.p3) / 2 - midpoint(panel), normal_vector(panel))

function transform_panel(panel :: AbstractPanel3D, point :: Point3D)
    T = get_transformation(panel)
    return T(panel), T(point) 
end

"""
    get_transformation(panel :: AbstractPanel3D)

Generate the mapping to transform a point from global coordinates to an `AbstractPanel3D`'s local coordinate system.
"""
get_transformation(p, P = I(3)) = LinearMap(P * local_coordinate_system(p)') ∘ Translation(-midpoint(p))


"""
    wake_panel(panels :: AbstractPanel3D, bound, α)

Calculate required transformation from GCS to panel LCS.
"""
function wake_panel(panels :: DenseArray{<: AbstractPanel3D}, bound, U)
    pt1 = 0.5 * ( p1(first(panels)) + p2(last(panels)) )
    pt4 = 0.5 * ( p4(first(panels)) + p3(last(panels)) )
    dr  = U * bound
    pt2 = pt1 + Point3D(dr...)
    pt3 = pt4 + Point3D(dr...)

    return WakePanel3D(pt1, pt2, pt3, pt4)
end

function panel_area_exact(panel :: AbstractPanel3D)
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