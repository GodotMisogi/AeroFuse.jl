module MathTools

using Base.Iterators

# Copying NumPy's linspace function
linspace(min, max, step) = min:(max - min)/step:max

#-------------HASKELL MASTER RACE--------------#

span(pred, iter) = (takewhile(pred, iter), dropwhile(pred, iter))
splitat(n, xs) = (xs[1:n,:], xs[n+1:end,:])  
lisa(pred, iter) = span(!pred, iter)

# Zieg Heil!


#--------------------HACKS----------------------#

"""
"Lenses" to access subfields on lists of objects.
"""
|>(obj, fields :: Array{Symbol}) = foldl(getproperty, fields, init = obj)
|>(list_objs :: Array{T}, fields :: Array{Symbol}) where {T <: Any} = list_objs .|> [fields]

field << obj = getfield(obj, field)

# Convert homogeneous struct entries to lists
structtolist(x) = [ name << x for name ∈ (fieldnames ∘ typeof)(x) ]

#--------------------------Convenient math------------------------#

⊗(A, B) = kron(A, B)
dot(V₁,V₂) = sum(V₁ .* V₂)
×(xs, ys) = (collect ∘ zip)(xs' ⊗ (ones ∘ length)(ys), (ones ∘ length)(xs)' ⊗ ys)

# Transforms (x, y) to the coordinate system with (x_s, y_s) as origin oriented at α_s.
affine_2D(x, y, x_s, y_s, α_s) = rotation(x - x_s, y - y_s, α_s)
inverse_rotation(x, y, angle) = (x * cos(angle) - y * sin(angle), x * sin(angle) + y * cos(angle))
rotation(x, y, angle) = (x * cos(angle) + y * sin(angle), -x * sin(angle) + y * cos(angle))

slope(x1, y1, x2, y2) = (y2 - y1)/(x2 - x1)

#---------------------Improving readiblity with arrays------------------------#

stencil(xs, n) = [ xs[n+1:end] xs[1:length(xs) - n] ]
parts(xs) = let adj = stencil(xs, 1); adj[1,:], adj[end,:] end

# Central differencing schema for pairs except at the trailing edge
function midgrad(xs) 
    first_two_pairs, last_two_pairs = permutedims.(parts(xs))
    central_diff_pairs = stencil(xs, 2)
    
    [first_two_pairs; central_diff_pairs; last_two_pairs]
end

# Difference operators
fwdsum(xs) = xs[2:end] .+ xs[1:end-1]
fwddiff(xs) = xs[2:end] .- xs[1:end-1]
ord2diff(xs) = xs[3:end] .- 2 * xs[2:end-1] .+ xs[1:end-2] 


adj3(xs) = [ xs[1:end-2,:] xs[2:end-1,:] xs[3:end,:] ]


end