# Matrix setup
#==========================================================================================#

influence_coefficient(source :: AbstractPotentialFlowElement, target :: AbstractPotentialFlowElement) = dot(velocity(control_point(target), source, one(eltype(control_point(target)))), normal_vector(target))

influence_coefficient(source :: AbstractPotentialFlowElement, target :: AbstractPotentialFlowElement, V_hat) = dot(velocity(control_point(target), source, 1., V_hat), normal_vector(target))


"""
    influence_matrix(elements)
    influence_matrix(elements, u_hat)

Assemble the Aerodynamic Influence Coefficient (AIC) matrix given an array of `AbstractPotentialFlowElement` and the freestream direction ``û``.
"""
influence_matrix(elements) = [ influence_coefficient(source, target) for target ∈ elements, source ∈ elements ]
influence_matrix(elements, V_hat) = [ influence_coefficient(source, target, V_hat) for target ∈ elements, source ∈ elements ]

function influence_matrix!(A, elements)
    size(A) == (length(elements), length(elements)) ||
        throw(DimensionMismatch("Influence matrix must have one row and column per element."))
    for j in axes(A, 2), i in axes(A, 1)
        A[i,j] = influence_coefficient(elements[j], elements[i])
    end
    return A
end

"""
    boundary_condition(elements, U, Ω)

Assemble the boundary condition vector given an array of `AbstractPotentialFlowElement`, the freestream velocity ``U``, and a quasi-steady rotation vector  ``Ω``.
"""
boundary_condition(elements, U, Ω) = map(hs -> dot(U + Ω × control_point(hs), normal_vector(hs)), elements)

# Added velocity for slipstream model
boundary_condition(elements, U, Ups, Ω) = map((hs, Up) -> dot(U + Ω × control_point(hs) + Up, normal_vector(hs)), elements, Ups)

"""
    apply_bc_rows!(AIC, boco, elements, U, Ups, Ω)

After the generic `velocity·normal` AIC and boundary-condition vector are assembled, give each
collocation element a chance to rewrite its own row and right-hand side to impose a
non-standard boundary condition. Dispatched per element via [`apply_bc_row!`]; standard Neumann
elements leave their row untouched, so this is a no-op unless the system carries an element
type that overrides it (e.g. the slender-body cylinder condition of a `FuselageLine`).
"""
function apply_bc_rows!(AIC, boco, elements, U, Ups, Ω)
    for i in eachindex(elements)
        apply_bc_row!(AIC, boco, i, elements[i], elements, U, Ups, Ω)
    end
    return AIC, boco
end

# Default: the standard Neumann row already assembled stands unchanged.
apply_bc_row!(AIC, boco, i, source :: AbstractPotentialFlowElement, elements, U, Ups, Ω) = nothing

"""
    has_bc_override(element :: AbstractPotentialFlowElement)

Whether `element` overrides [`apply_bc_row!`] with a non-standard boundary condition. Element
types defining such an override must also return `true` here, so that device solvers know which
rows to rewrite without downloading the full matrix.
"""
has_bc_override(::AbstractPotentialFlowElement) = false

## Backend dispatch
#==========================================================================================#

"""
    component_blocks(elements)

Split a named `ComponentArray` of potential-flow elements into `(range, block)` pairs, where
`range` indexes the block in the flattened array and `block` is a concretely typed vector of its
elements. A plain array of elements is returned as a single block.
"""
function component_blocks(elements :: ComponentArray)
    offset = 0
    blocks = map(collect(keys(elements))) do key
        block = [ el for el in vec(getproperty(elements, key)) ] # Narrow to the concrete element type
        range = offset + 1 : offset + length(block)
        offset += length(block)
        (range, block)
    end
    offset == length(elements) || throw(DimensionMismatch("Component blocks do not tile the element array."))
    return filter(b -> !isempty(b[2]), blocks)
end

component_blocks(elements :: AbstractArray{<: AbstractPotentialFlowElement}) = [ (1:length(elements), [ el for el in vec(elements) ]) ]

"""
    device_solve_linear(backend, elements, U, Ups, Ω)

Assemble and solve the influence system of [`solve_linear`](@ref) on a compute `backend`.
`backend = nothing` runs the serial host implementation; `KernelAbstractions.Backend`
methods are provided by the `KernelAbstractions` package extension.
"""
device_solve_linear(::Nothing, elements, U, Ups, Ω) = solve_linear(elements, U, Ups, Ω)

"""
    device_induced_sum(backend, f, points, elements, strengths, V_hat)

Evaluate ``∑_j f(r_i, e_j, s_j, V̂)`` at every point ``r_i`` of `points` over `elements` with
`strengths`, returning a vector of velocities in the layout of `points`. `f` is an element
velocity kernel such as `velocity` or `trailing_velocity`. `backend = nothing` runs on the host;
`KernelAbstractions.Backend` methods are provided by the package extension.
"""
device_induced_sum(::Nothing, f :: F, points, elements, strengths, V_hat) where F =
    map(r -> sum(x -> f(r, x[1], x[2], V_hat), zip(elements, strengths)), points)

# Matrix-free setup for nonlinear analyses
#==========================================================================================#

induced_velocity(r, elements, strengths, U_hat) = sum(x -> velocity(r, x[1], x[2], U_hat), zip(elements, strengths))

induced_trailing_velocity(r, elements, strengths, U_hat) = sum(x -> trailing_velocity(r, x[1], x[2], U_hat), zip(elements, strengths))

# In-place versions
@views function induced_velocity!(vel, r, elements, strengths, U_hat)
    for i in eachindex(elements)
        vel += velocity(r, elements[i], strengths[i], U_hat)
    end

    return vel
end

@views function induced_trailing_velocity!(vel, r, elements, strengths, U_hat)
    for i in eachindex(elements)
        vel += trailing_velocity(r, elements[i], strengths[i], U_hat)
    end

    return vel
end

function induced_velocity(r, hs, strengths, U, Ω)
    vel = zero(r)
    induced_velocity!(vel, r, hs, strengths, -normalize(U)) - (U + Ω × r)
    # @timeit "Induced Velocity" induced_velocity(r, hs, strengths, -normalize(U)) - (U + Ω × r)
end

function induced_trailing_velocity(r, elements, strengths, U, Ω) 
    vel = zero(r)
    induced_trailing_velocity!(vel, r, elements, strengths, -normalize(U)) - (U + Ω × r)
    # induced_trailing_velocity(r, elements, strengths, -normalize(U)) - (U + Ω × r)
end

# Matrix-free residuals evaluate the standard normal-velocity condition only;
# prescribed fields and element-specific boundary-row overrides are not included.
# Residual computations
residual(r, n, hs, strengths, U, Ω) = dot(induced_velocity(r, hs, strengths, U, Ω), n)

solve_nonlinear(elements, strengths, U_hat, Ω_hat) = map(hs -> residual(control_point(hs), normal_vector(hs), elements, strengths, U_hat, Ω_hat), elements)

solve_nonlinear!(R, elements, strengths, U_hat, Ω_hat) = map!(hs -> residual(control_point(hs), normal_vector(hs), elements, strengths, U_hat, Ω_hat), R, elements)
