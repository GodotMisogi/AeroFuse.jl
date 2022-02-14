## Prandtl-Glauert transformation
#======================================================#

prandtl_glauert_scale_coordinates(x, y, z, β) = SVector(x, β * y, β * z)
prandtl_glauert_scale_coordinates(r, β) = @views prandtl_glauert_scale_coordinates(r[1], r[2], r[3], β)

prandtl_glauert_scale_normal(nx, ny, nz, β) = SVector(β * nx, ny, nz)
prandtl_glauert_scale_normal(n, β) = @views prandtl_glauert_scale_normal(n[1], n[2], n[3], β)

prandtl_glauert_scale_coordinates(horseshoe :: Horseshoe, β) = 
    setproperties(horseshoe,
                  r1                = prandtl_glauert_scale_coordinates(horseshoe.r1, β), 
                  r2                = prandtl_glauert_scale_coordinates(horseshoe.r2, β), 
                  collocation_point = prandtl_glauert_scale_coordinates(horseshoe.collocation_point, β),
                  normal            = prandtl_glauert_scale_normal(horseshoe.normal, β),
                )

prandtl_glauert_inverse_scale_velocities(vx, vy, vz, β) = SVector(vx / β^2, vy / β, vz / β)
prandtl_glauert_inverse_scale_velocities(v, β) = @views prandtl_glauert_inverse_scale_velocities(v[1], v[2], v[3], β)