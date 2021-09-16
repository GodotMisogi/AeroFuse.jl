## Aerodynamic models

# Tangler-Ostowari model
function post_stall_CL(α, AR, CL_s, α_s)
    # Processing
    sα, cα     = sincos(α)
    sα_s, cα_s = sincos(α_s)

    # Parameters
    C1 = 1.1 + 0.018AR
    A1 = C1 / 2
    A2 = (CL_s - C1 * sα_s * cα_s) * sα_s / cα_s^2

    # Lift coefficient
    CL = A1 * 2sα * cα + A2 * cα^2 / sα
end

function post_stall_CD(α, AR, CD_s, α_s, foil_tc)
    # Processing
    sα, cα     = sincos(α)
    sα_s, cα_s = sincos(α_s)

    # Parameters
    CD_max = (1 + 0.065AR) / (0.9 + foil_tc)
    B1 = CD_max
    B2 = (CD_s - CD_max * sα_s) / cα_s

    # Drag coefficient
    CD = B1 * sα + B2 * cα
end

lifting_line_slope(CLα, AR, e) = CLα / (1 + CLα/(π * e * AR))
induced_drag(CL, AR, e) = CL^2 / (π * e * AR)

biplane_induced_drag(L, q, b, σ) = (1 + σ)/(π * q) * L^2 / b^2

quartic_CD(α) = 0.008 + 1.107α^2 + 1.792α^4


struct Airfoil
    cl_func
    cd_func
    cm_func
end

cl(af :: Airfoil, args...) = af.cl_func(args...)
cd(af :: Airfoil, args...) = af.cd_func(args...)
cm(af :: Airfoil, args...) = af.cm_func(args...)


# Post-stall modified functions
# global_cl(α, Cl_func :: F, Cl_s, α_s, AR) where F = ifelse(α < α_s, Cl_func, post_stall_CL(α, AR, Cl_s, α_s))
# global_cd(α, Cd_func :: F, Cd_s, α_s, AR, foil_tc) where F = ifelse(α < α_s, Cd_func, post_stall_CD(α, AR, Cd_s, α_s, foil_tc))

# Basic lift/drag coefficient relations
linear_cl(α, Cl_max, Cl_α) = Cl_max - Cl_α * α
parabolic_cd(α, Cd_max, Cd_0, k) = k * (α - Cd_0)^2 + Cd_0 # ???

global_cl(α, Cl_α, Cl_s, α_s, AR) = ifelse(α < α_s, linear_cl(α, Cl_s, Cl_α), post_stall_CL(α, AR, Cl_s, α_s))
global_cd(α, Cd_0, Cd_s, α_s, AR, foil_tc, k) = ifelse(α < α_s, parabolic_cd(α, Cd_s, Cd_0, k), post_stall_CD(α, AR, Cd_s, α_s, foil_tc))