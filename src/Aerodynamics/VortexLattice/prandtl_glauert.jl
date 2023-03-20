## Prandtl-Glauert transformation
#======================================================#

prandtl_glauert_scale_coordinates(x, y, z, β) = SVector(x, β * y, β * z)
prandtl_glauert_scale_coordinates(r, β) = @views prandtl_glauert_scale_coordinates(r[1], r[2], r[3], β)

prandtl_glauert_scale_normal(nx, ny, nz, β) = SVector(β * nx, ny, nz)
prandtl_glauert_scale_normal(n, β) = @views prandtl_glauert_scale_normal(n[1], n[2], n[3], β)

prandtl_glauert_scale_coordinates(horseshoe :: Horseshoe, β) = setproperties(horseshoe,
    r1 = prandtl_glauert_scale_coordinates(horseshoe.r1, β),
    r2 = prandtl_glauert_scale_coordinates(horseshoe.r2, β),
    rc = prandtl_glauert_scale_coordinates(horseshoe.rc, β),
    normal = prandtl_glauert_scale_normal(horseshoe.normal, β),
)

prandtl_glauert_scale_coordinates(ring :: VortexRing, β) = setproperties(ring,
    r1 = prandtl_glauert_scale_coordinates(ring.r1, β),
    r2 = prandtl_glauert_scale_coordinates(ring.r2, β),
    r3 = prandtl_glauert_scale_coordinates(ring.r3, β),
    r4 = prandtl_glauert_scale_coordinates(ring.r4, β),
    rc = prandtl_glauert_scale_coordinates(ring.rc, β),
    normal = prandtl_glauert_scale_normal(ring.normal, β),
)

prandtl_glauert_inverse_scale_velocity(vx, vy, vz, β) = SVector(vx / β^2, vy / β, vz / β)
prandtl_glauert_inverse_scale_velocity(v, β) = @views prandtl_glauert_inverse_scale_velocity(v[1], v[2], v[3], β)