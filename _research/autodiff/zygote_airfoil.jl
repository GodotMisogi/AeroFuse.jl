using Zygote
using StaticArrays
using AeroMDAO

##
@Zygote.adjoint (T :: Type{<:SVector})(xs :: Number ...) = T(xs...), dv -> (nothing, dv...)
@Zygote.adjoint (T :: Type{<:SVector})(x :: AbstractVector) = T(x), dv -> (nothing, dv)

## Custom adjoints
@Zygote.adjoint AeroMDAO.PanelGeometry.p1(p :: Panel2D) = p.p1, x̄ -> (Panel2D(x̄, SVector(0., 0.)),)
@Zygote.adjoint AeroMDAO.PanelGeometry.p2(p :: Panel2D) = p.p2, ȳ -> (Panel2D(SVector(0., 0.), ȳ),)
@Zygote.adjoint Panel2D(a, b) = Panel2D(a, b), p̄ -> (p̄.p1, p̄.p2)

## 2D doublet-source panel method
function naca_zygote(digits)
    airfoil = (Foil ∘ naca4)(digits)
    uniform = Uniform2D(1.0, 5.0)
    solve_case(airfoil, uniform, num_panels = 60)
end

##
digits = (2,4,1,2)
@time naca_zygote(digits)

##
gradient(naca_zygote, digits)