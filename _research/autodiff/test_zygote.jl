##
using Zygote
import Base: +, -, zero
import Base.Iterators
using StaticArrays

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

