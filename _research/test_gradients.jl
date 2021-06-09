## 
using Revise

##
using AeroMDAO
using StaticArrays
import Base.Iterators
import Base: getindex

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
@Zygote.adjoint AeroMDAO.p1(p :: AeroMDAO.Panel2D) = p.p1, p̄1 -> (AeroMDAO.Panel2D(p̄1, SVector(zeros(2)...)),)
@Zygote.adjoint AeroMDAO.p2(p :: AeroMDAO.Panel2D) = p.p2, p̄2 -> (AeroMDAO.Panel2D(SVector(zeros(2)...), p̄2),)

##
function optimize_CST(alpha_u, alpha_l)
  airfoil = (Foil ∘ kulfan_CST)(alpha_u, alpha_l, (0., 0.), 0., 80)
  uniform = Uniform2D(1.0, 5.0)
  solve_case(airfoil, uniform, num_panels = 60)
end

##
alpha_u = [0.2, 0.3, 0.2, 0.15, 0.2]
alpha_l = [-0.2, -0.1, -0.1, -0.001]
gradient(optimize_CST, alpha_u, alpha_l)
