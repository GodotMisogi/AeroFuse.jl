module MathTools

using StaticArrays
using Base.Iterators
using Base: product
using Interpolations

struct Point2D{T <: Real} <: FieldVector{2, T} 
    x :: T
    y :: T
end

StaticArrays.similar_type(::Type{<:Point2D}, ::Type{T}, s::Size{(2,)}) where {T} = Point2D{T}

x(p :: Point2D) = p.x
y(p :: Point2D) = p.y

struct Point3D{T <: Real} <: FieldVector{3, T} 
    x :: T
    y :: T
    z :: T
end

StaticArrays.similar_type(::Type{<:Point3D}, ::Type{T}, s::Size{(3,)}) where {T} = Point3D{T}

x(p :: Point3D) = p.x
y(p :: Point3D) = p.y
z(p :: Point3D) = p.z

# Copying NumPy's linspace function
linspace(min, max, step) = min:(max - min)/step:max
columns(M) = tuple((view(M, :, i) for i in 1:size(M, 2))...)

## Haskell Master Race
#===========================================================================#

span(pred, iter) = (takewhile(pred, iter), dropwhile(pred, iter))
splitat(n, xs) = @views (xs[1:n,:], xs[n+1:end,:])  

lisa(pred, iter) = span(!pred, iter)

Base.Iterators.partition(pred, xs, f, g) = f.(filter(pred, xs)), g.(filter(!pred, xs))
Base.Iterators.partition(pred, xs, f = identity) = partition(pred, xs, f, f)

## Renaming math operations
#===========================================================================#

"""
"Lenses" to access subfields on lists of objects.
"""
|>(obj, fields :: Array{Symbol}) = foldl(getproperty, fields, init = obj)
|>(list_objs :: Array{T}, fields :: Array{Symbol}) where T <: Any = list_objs .|> [fields]

field << obj = getfield(obj, field)

# Convert homogeneous struct entries to lists
structtolist(x) = [ name << x for name ∈ (fieldnames ∘ typeof)(x) ]

## Renaming math operations
#===========================================================================#

⊗(A, B)    = kron(A, B)

×(xs, ys)   = product(xs, ys)
dot(V₁, V₂) = sum(V₁ .* V₂)
# ×(xs, ys) = (collect ∘ zip)(xs' ⊗ (ones ∘ length)(ys), (ones ∘ length)(xs)' ⊗ ys)

## Basic transformations/geometry
#===========================================================================#

# Transforms (x, y) to the coordinate system with (x_s, y_s) as origin oriented at α_s.
affine_2D(x, y, x_s, y_s, α_s) = rotation(x - x_s, y - y_s, α_s)
inverse_rotation(x, y, angle) = SVector(x * cos(angle) - y * sin(angle), x * sin(angle) + y * cos(angle))
rotation(x, y, angle) = SVector(x * cos(angle) + y * sin(angle), -x * sin(angle) + y * cos(angle))

# Matrix versions
rotation(θ)         = [ cos(θ) sin(θ) ;
                       -sin(θ) cos(θ) ]
inverse_rotation(θ) = rotation(-θ)

# Cartesian-polar coordinates
magnitude(U, W) = sqrt(U^2 + W^2)
angle(U, W) = atan(W, U)


slope(x1, y1, x2, y2) = (y2 - y1) / (x2 - x1)

## Array conversions
#===========================================================================#

tuparray(xs)  = tuple.(eachcol(xs)...)
vectarray(xs) = SVector.(eachcol(xs)...)

extend_yz(coords) = @views [ coords[:,1] zero(coords[:,1]) coords[:,2] ]

reflect_mapper(f, xs) = @views [ f(xs[:,end:-1:1]) xs ]

## Difference operations
#===========================================================================#

forward_difference_matrix(n) = [ I zeros(n) ] - [ zeros(n) I ]

forward_sum(xs) = @views @. xs[2:end] + xs[1:end-1]
forward_difference(xs) = @views @. xs[2:end] - xs[1:end-1]
forward_division(xs) = @views @. xs[2:end] / xs[1:end-1]
ord2diff(xs) = @views @. xs[3:end] - 2 * xs[2:end-1] + xs[1:end-2] 

adj3(xs) = @views zip(xs[1:end-2,:], xs[2:end-1,:], xs[3:end,:])

# Central differencing schema for pairs except at endpoints
function midpair_map(f :: H, xs; dims :: Int64) where {H}
    if dims == 1
        @views  [permutedims(f.(xs[1,:], xs[2,:]));     f.(xs[1:end-2,:], xs[3:end,:]);     permutedims(f.(xs[end-1,:], xs[end,:])) ]
    elseif dims == 2
        @views  [f.(xs[:,1], xs[:,2])       f.(xs[:,1:end-2], xs[:,3:end])      f.(xs[:,end-1], xs[:,end])  ]
    else
        ArgumentError("Array with order > 2 is currently not supported!")
    end
end

# stencil(xs, n) = [ xs[n+1:end] xs[1:length(xs) - n] ]
# parts(xs) = let adj = stencil(xs, 1); adj[1,:], adj[end,:] end

# Lazy? versions
# stencil(xs, n) = zip(xs[n+1:end], xs[1:length(xs) - n])
# parts(xs) = (first ∘ stencil)(xs, 1), (last ∘ stencil)(xs, 1)

# function midgrad(xs) 
#     first_two_pairs, last_two_pairs = permutedims.(parts(xs))
#     central_diff_pairs = stencil(xs, 2)
    
#     [first_two_pairs; central_diff_pairs; last_two_pairs]
# end


## Spacing formulas
#===========================================================================#

uniform_spacing(x1, x2, n) = range(x1, x2, length = n)
linear_spacing(x_center, len, n :: Integer) = @. x_center + len * 0:1/(n-1):1
cosine_spacing(x_center, diameter, n :: Integer = 40) = @. x_center + (diameter / 2) * cos(-π:π/(n-1):0)

function sine_spacing(x1, x2, n :: Integer = 40)
    d = x2 - x1
    if n < 0
       @. x2 - d * sin(π/2 * (1. - (1:1/(n+1):0)))[end:-1:1]
    else
       @. x1 + d * sin(π/2 * (0:1/(n-1):1))
    end
end

function cosine_interp(coords, n :: Integer = 40)
    xs = @views coords[:,1]
    ys = @views coords[:,2]

    d = maximum(xs) - minimum(xs)
    x_center = (maximum(xs) + minimum(xs)) / 2
    x_circ = cosine_spacing(x_center, d, n)
    
    itp_circ = LinearInterpolation(xs, ys)
    y_circ = itp_circ(x_circ)

    [ x_circ y_circ ]
end

## Iterator methods
#===========================================================================#

# Need to improve this via recursion; the for loop seems really unnecessary
function accumap(f, n, xs)
    data = [ xs ]
    for i = 1:n
        ys = map(f, xs)
        data = [ reduce(vcat, data); ys ]
        xs = ys
    end
    return reduce(hcat, data)
end

## Helper functions
#===========================================================================#

"""
    weighted(x1, x2, μ)

Compute the weighted value between two values ``x_1`` and ``x_2`` with weight ``\\mu \\in [0,1]``.
"""
weighted(x1, x2, μ) = (1 - μ) * x1 + μ * x2

"""
Compute the weighted average (μ) of two vectors. 
"""
weighted_vector(x1, x2, μ) = weighted.(x1, x2, μ)

## Macros
#===========================================================================#

# macro getter(obj)
#     return :((function name(x)
#                 getfield(name, x) 
#             end))
# end

reshape_array(arr, inds, sizes) = @views [ reshape(arr[i1+1:i2], size) for (i1, i2, size) in zip(inds[1:end-1], inds[2:end], sizes) ]

end