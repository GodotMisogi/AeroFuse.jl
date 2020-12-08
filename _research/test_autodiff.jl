##
import Base: +, -, zero
import Base.Iterators
using StaticArrays
import Base.Iterators.product

struct Point
  x :: Real
  y :: Real
end

width(p::Point) = p.x
height(p::Point) = p.y

dist(p::Point) = sqrt(width(p)^2 + height(p)^2)

##
not_perimeter(p1 :: Point, p2 :: Point) = p1 === p2 ? 0.5 : 1 / width(p1)  + 1 / height(p2)

function sumthing(x1, x2)
    pts = Point.(x1, x2)
    prod_pts = product(pts,pts)
    (sum ∘ map)(pts -> not_perimeter(pts...), prod_pts)
end

##
distmat(ps :: AbstractVector{Point}) = [ not_perimeter(p_i, p_j) for p_i in ps, p_j in ps ]

function funky(xs :: AbstractVector{<: Real}, ys :: AbstractVector{<: Real}, rhs :: AbstractVector{<: Real})
  ps = Point.(xs, ys)
  inf = distmat(ps) 
  boco = rhs

  φs = inf \ boco
end

##
using ReverseDiff: GradientTape, GradientConfig, gradient, gradient!, compile, DiffResults

## Test 1
#========================================#
inputs = ([1., 3., 1.], [2., 4., 5.], [1., 3., 5.])

const test_tape = GradientTape(sum ∘ funky, inputs)
const compiled_test_tape = compile(test_tape)

results = similar.(inputs)
all_results = map(DiffResults.GradientResult, results)
cfg = GradientConfig(inputs)

gradient!(results, compiled_test_tape, inputs)

## Test 2
#========================================#
inputs = (1:5, 5:9)

const f_tape = GradientTape(sumthing, inputs)
const compiled_f_tape = compile(f_tape)

results = (similar(1:5), similar(5:9))
all_results = map(DiffResults.GradientResult, results)
cfg = GradientConfig(inputs)

gradient!(results, compiled_f_tape, inputs)