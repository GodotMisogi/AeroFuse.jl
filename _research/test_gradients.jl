## 
using Revise

##
using AeroMDAO
using BenchmarkTools
using ProfileView
using StaticArrays
import Base.Iterators
import Base: getindex

# getindex(::Nothing, :: Integer) = nothing

## 2D doublet-source panel method
alpha_u = [0.2, 0.3, 0.2, 0.15, 0.2]
alpha_l = [-0.2, -0.1, -0.1, -0.001]
dzs = (1e-4, 1e-4)
airfoil = kulfan_CST(alpha_u, alpha_l, (0., 0.), 0., 80)
panels = make_2Dpanels(airfoil);


uniform = Uniform2D(1.0, 5.0)

##
function optimize_CST(α)
    cl = solve_case(panels, Uniform2D(1.0, α))
end

optimize_CST(5.0)


## Gradient package
using Zygote

@Zygote.adjoint (T::Type{<:SVector})(xs::Number...) = T(xs...), dv -> (nothing, dv...)
@Zygote.adjoint (T::Type{<:SVector})(x::AbstractVector) = T(x), dv -> (nothing, dv)

@Zygote.adjoint enumerate(xs) = enumerate(xs), diys -> (map(last, diys),)

_ndims(::Base.HasShape{d}) where {d} = d
_ndims(x) = Base.IteratorSize(x) isa Base.HasShape ? _ndims(Base.IteratorSize(x)) : 1

@Zygote.adjoint function Iterators.product(xs...)
                    d = 1
                    Iterators.product(xs...), dy -> ntuple(length(xs)) do n
                        nd = _ndims(xs[n])
                        dims = ntuple(i -> i<d ? i : i+nd, ndims(dy)-nd)
                        d += nd
                        reshape(sum(y->y[n], dy; dims=dims), axes(xs[n]))
                    end
                end
                
@Zygote.adjoint function Iterators.Zip(xs)
                    back(dy::NamedTuple{(:is,)}) = tuple(dy.is)
                    back(dy::AbstractArray) = ntuple(length(xs)) do d
                      dx = map(y->y[d], dy)
                      length(dx) == length(xs[d]) ? dx : vcat(dx, falses(length(xs[d])-length(dx)))
                    end |> tuple
                    Iterators.Zip(xs), back
                  end

@Zygote.adjoint AeroMDAO.Panel2D(p1, p2) = AeroMDAO.Panel2D(p1, p2), p̄ -> (p̄.p1, p̄.p2)
@Zygote.adjoint AeroMDAO.point1(p :: AeroMDAO.Panel2D) = p.p1, p̄1 -> (AeroMDAO.Panel2D(p̄1, SVector(zeros(2)...)),)
@Zygote.adjoint AeroMDAO.point2(p :: AeroMDAO.Panel2D) = p.p2, p̄2 -> (AeroMDAO.Panel2D(SVector(zeros(2)...), p̄2),)

##
gradient(optimize_CST, 5.0)




## 3D vortex-lattice method
num_secs = 3

foil = naca4((4,4,1,2))
foils = [ foil for i ∈ 1:num_secs ]
airfoils = Foil.(foils)
wing_twists = [2., 0., -2.]
wing_spans = [0.5, 0.5]
wing_dihedrals = [0, 11.3]
wing_sweeps = [1.14, 8]

ρ = 1.225
uniform = Uniform(10.0, 5.0, 0.0)

## Define objective
function test_model(wing_chords)
    wing_right = HalfWing(airfoils, wing_chords, wing_spans, wing_dihedrals, wing_sweeps, wing_twists)
    wing = Wing(wing_right, wing_right)

    ## Assembly
    ref = (0.25 * mean_aerodynamic_chord(wing), 0, 0)
    coeffs = solve_case(wing, uniform, ref, span_num = 10, chord_num = 5, print = false)

    return coeffs
end

lift_coefficient(chords) = test_model(chords)[1]
drag_coefficient(chords) = test_model(chords)[2]
lift_to_drag_ratio(chords) = lift_coefficient(chords)/drag_coefficient(chords)

##
lift_coefficient([0.16, 0.12, 0.08])

##
lift_coefficient([0.2, 0.1, 0.05])

##
lift_coefficient'([0.16, 0.12, 0.08])
