## Linearized stability analysis
#==========================================================================================#

struct AircraftData{T}
    dvs :: Matrix{T}
    refs :: References{T}
    m :: T
    I :: SMatrix{3,3,T}
end

"""
    longitudinal_stability_derivatives(dvs, U, m, Iyy, Q, S, c)

Compute the stability derivatives for the forces and moments in the longitudinal plane ``[X,Z,M]_{u,w,q}``.

The inputs are force and moment coefficient stability derivatives matrix `dvs`, the freestream speed `U`, the mass `m` and longitudinal moment of inertia ``I_{yy}``, the dynamic pressure Q, reference area ``S`` and chord length ``c``.
"""
function longitudinal_stability_derivatives(dvs, U, m, Iyy, Q, S, c, a = 330.)
    QS  = Q * S
    c1  = QS / (m * U) 
    c2  = QS / (Iyy * U)

    X_u = c1 * dvs.CX_Ma / a
    Z_u = c1 * dvs.CZ_Ma / a
    M_u = c2 * dvs.Cm_Ma / a * c

    X_w = c1 * dvs.CX_al
    Z_w = c1 * dvs.CZ_al
    M_w = c2 * dvs.Cm_al * c

    X_q = c1 / 2 * dvs.CX_qb
    Z_q = c1 / 2 * dvs.CZ_qb
    M_q = c2 / 2 * dvs.Cm_qb * c

    X_u, Z_u, M_u, X_w, Z_w, M_w, X_q, Z_q, M_q
end

longitudinal_stability_matrix(X_u, Z_u, M_u, X_w, Z_w, M_w, X_q, Z_q, M_q, U₀, g) = 
    [ X_u X_w     0      -g
      Z_u Z_w (U₀ + Z_q)  0
      M_u M_w    M_q      0
       0   0      1       0 ]

"""
    lateral_stability_derivatives(dvs, U, m, Iyy, Q, S, c)

Compute the stability derivatives for the forces and moments in the lateral plane ``[Y,L,N]_{v,p,r}``.

The inputs are force and moment coefficient stability derivatives matrix `dvs`, the freestream speed `U`, the mass `m` and lateral moment sof inertia ``I_{xx}, I_{zz}``, the dynamic pressure Q, reference area ``S`` and span length ``b``.
"""
function lateral_stability_derivatives(dvs, U, m, Ixx, Izz, Q, S, b)
    QS  = Q * S
    c1  = QS / (m * U) 
    c2  = QS / (Ixx * U)
    c3  = QS / (Izz * U)

    Y_v = c1 * dvs.CY_be
    L_v = c2 * dvs.Cl_be * b
    N_v = c3 * dvs.Cn_be * b

    Y_p = c1 / 2 * dvs.CY_pb * b
    L_p = c2 / 2 * dvs.Cl_pb * b^2
    N_p = c3 / 2 * dvs.Cn_pb * b^2

    Y_r = c1 / 2 * dvs.CY_rb * b
    L_r = c2 / 2 * dvs.Cl_rb * b^2
    N_r = c3 / 2 * dvs.Cn_rb * b^2

    Y_v, L_v, N_v, Y_p, L_p, N_p, Y_r, L_r, N_r
end

lateral_stability_matrix(Y_v, L_v, N_v, Y_p, L_p, N_p, Y_r, L_r, N_r, U₀, θ₀, g) =
    [ Y_v Y_p (Y_r - U₀) (g * cos(θ₀))
      L_v L_p     L_r         0      
      N_v N_p     N_r         0      
       0   1       0          0       ]