abstract type AbstractLine end

"""
    Line(r1 :: SVector{3,<: Real}, r2 :: SVector{3,<: Real})

Define a line segment based on two 3-dimensional vectors.
"""
struct Line{T <: Real} <: AbstractLine
    r1 :: SVector{3,T}
    r2 :: SVector{3,T}
end

Line(r1, r2) = let T = promote_type(eltype(r1), eltype(r2)); Line{T}(r1, r2) end
Line((r1, r2)) = Line(r1, r2)

r1(line :: Line) = line.r1
r2(line :: Line) = line.r2
vector(line :: Line) = r2(line) - r1(line)
center(line :: Line) = (r1(line) + r2(line)) / 2

points(lines) = [ r1.(lines); [(r2 ∘ last)(lines)] ]

transform(line :: Line, rotation, translation) = let trans = Translation(translation) ∘ LinearMap(rotation); Line((trans ∘ r1)(line), (trans ∘ r2)(line)) end

r1(r, line :: Line) = r - r1(line)
r2(r, line :: Line) = r - r2(line)