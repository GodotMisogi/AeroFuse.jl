# Induced velocities
induced_velocity(r, horseshoes, Γs, U_hat) = sum(x -> velocity(r, x[1], x[2], U_hat), zip(horseshoes, Γs))

velocity(hs, Γs, r, n, U_hat, Ω_hat) = dot(induced_velocity(r, hs, Γs, -U_hat) - (U_hat + Ω_hat × r), n)

# Residual computation
aerodynamic_residual!(R_A, hs, rs, nms, Γs, U_hat, Ω_hat) = R_A .= velocity.(Ref(hs), Ref(Γs), rs, nms, Ref(U_hat), Ref(Ω_hat))