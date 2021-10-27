using StaticArrays

struct BSpline{T <: Real, S <: Integer}
    t :: Vector{T}  # Knot vector
    k :: S          # Spline multiplicity
end

function BSpline(t :: AbstractVector{T}, k :: S) where T <: Real where S <: Integer 
    sort!(t)
    clamp = ones(k - 1)
    t_clamped = [ first(t) * clamp; t; last(t) * clamp ]
    BSpline{T,S}(t_clamped, k)
end

"""
    knot_basis(t, T, i, k)

Evaluates a basis element ``N_{i,k}`` at ``t`` with set of internal knots ``\\mathbf T``.
"""
function knot_basis(x, T, i, k) 
    if k == 1
        ifelse(T[i] <= x < T[i+1], 1., 0.)
    elseif i + k >= length(T) || i < 1 
        0.
    elseif k < 1
        @error "How did you know I can't count?"
    else
        w1 = spline_weight(x, T, i, k - 1) 
        w2 = spline_weight(x, T, i + 1, k - 1)
        w1 * knot_basis(x, T, i, k - 1) +  (1 - w2) * knot_basis(x, T, i + 1, k - 1)
    end
end

spline_weight(x, T, i, k) = ifelse(T[i+k] != T[i], (x - T[i]) / (T[i+k] - T[i]), 0.)

function knot_derivative(t, T, i, k)
    d1 = (T[i+k] - T[i])
    d2 = (T[i+k+1] - T[i+1])

    a = k * ifelse(d1 == 0, 0., d1)
    b = k * ifelse(d2 == 0, 0., d2)

    a * knot_basis(t, T, i, k - 1) + b * knot_basis(t, T, i + 1, k - 1)
end

function bspline(t, Ps, Ts, k = 2) 
    basis = knot_basis.(t, Ref(Ts), 1:length(Ps), k)
    sum(p_i .* N_ik for (p_i, N_ik) in zip(Ps, basis))
end

Spline(Ts, Ps, m) = x -> bspline(x, Ps, Ts, m)

## Tests
using Interpolations, Plots
x = sort(2*rand(10))
y = x.^2 .+ rand.()

itp = interpolate(x,y,FritschCarlsonMonotonicInterpolation())
xq = minimum(x):0.01:maximum(x)

plot(x,y, seriestype=:scatter, label="Data", legend=:topleft, framestyle=:box)
plot!(xq, itp.(xq), label="Interpolation")