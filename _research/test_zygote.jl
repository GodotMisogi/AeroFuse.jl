##
using Zygote
import Base: +, -, zero
import Base.Iterators
using StaticArrays
using Revise
using AeroMDAO

##
@Zygote.adjoint (T :: Type{<:SVector})(xs :: Number ...) = T(xs...), dv -> (nothing, dv...)
@Zygote.adjoint (T :: Type{<:SVector})(x :: AbstractVector) = T(x), dv -> (nothing, dv)

# @Zygote.adjoint enumerate(xs) = enumerate(xs), diys -> (map(last, diys),)

_ndims(::Base.HasShape{d}) where {d} = d
_ndims(x) = Base.IteratorSize(x) isa Base.HasShape ? _ndims(Base.IteratorSize(x)) : 1

@Zygote.adjoint function Iterators.product(xs...)
    back(::AbstractArray{Nothing}) = nothing
    function back(dy::AbstractArray)
        d = 1
        ntuple(length(xs)) do n
        first(dy)[n] === nothing && return nothing
        nd = _ndims(xs[n])
        dims = ntuple(i -> i<d ? i : i+nd, ndims(dy)-nd)
        d += nd
        return reshape(sum(Zygote.StaticGetter{n}(), dy; dims=dims), axes(xs[n]))
        end
    end
    Iterators.product(xs...), back
end

@Zygote.adjoint function Iterators.Zip(xs)
    back(dy::NamedTuple{(:is,)}) = tuple(dy.is)
    back(dy::AbstractArray) = ntuple(length(xs)) do d
      dx = map(y->y[d], dy)
      length(dx) == length(xs[d]) ? dx : vcat(dx, falses(length(xs[d])-length(dx)))
    end |> tuple
    Iterators.Zip(xs), back
end



## 2D doublet-source panel method
alpha_u = [0.2, 0.3, 0.2, 0.15, 0.2]
alpha_l = [-0.2, -0.1, -0.1, -0.001]
dzs = (1e-4, 1e-4)
airfoil = (Foil ∘ kulfan_CST)(alpha_u, alpha_l, (0., 0.), 0., 80)

## Objective function
function test_zygote(alpha_u, alpha_l)
    airfoil = (Foil ∘ kulfan_CST)(alpha_u, alpha_l, (0., 0.), 0., 60)
    uniform = Uniform2D(1.0, 5.0)
    solve_case(airfoil, uniform, 60)
end

##
@time test_zygote(alpha_u, alpha_l)

##
gradient(test_zygote, alpha_u, alpha_l)


##
struct Point
    x :: Real
    y :: Real
end
    
width(p::Point) = p.x
height(p::Point) = p.y

a::Point + b::Point = Point(width(a) + width(b), height(a) + height(b))
a::Point - b::Point = Point(width(a) - width(b), height(a) - height(b))
dist(p::Point) = sqrt(width(p)^2 + height(p)^2)
zero(:: Point) = Point(0, 0)
                

@Zygote.adjoint width(p::Point) = p.x, x̄ -> (Point(x̄, 0),)
@Zygote.adjoint height(p::Point) = p.y, ȳ -> (Point(0, ȳ),)
@Zygote.adjoint Point(a, b) = Point(a, b), p̄ -> (p̄.x, p̄.y)


##
xs = Point.((1:4), (2:5))
not_perimeter(p1 :: Point, p2 :: Point) = ifelse(p1 === p2, 0.5, 1 / width(p1)  + 1 / height(p2))

np(xs) = sum([ not_perimeter(p1, p2) for p1 in xs, p2 in xs ])

##
not_perimeter(xs[1], xs[2])

##
np(xs)

##
test(n) = sum([ i+j for i in 1:n^2, j in 1:n ])
test(5)

gradient(test, 5)
##
gradient(np, xs)
# grad_pts = gradient(xs -> sum(map(pts -> not_perimeter(pts...), product(xs, xs))), xs)

