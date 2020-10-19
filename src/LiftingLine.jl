module LiftingLine

include("AeroMDAO.jl")
include("MathTools.jl")
include("PanelMethods.jl")

using .AeroMDAO: HalfWing, Wing
using .MathTools: <<
using LinearAlgebra

abstract type Laplace end

struct Uniform3D <: Laplace
    mag :: Float64
    alpha :: Float64 
    beta :: Float64
end

velocity(uni :: Uniform3D) = let ang = deg2rad(uni.alpha), beta = deg2rad(uni.beta); 
    uni.mag .* (cos(alpha) * cos(beta), cos(alpha) * sin(beta), sin(alpha)) end

struct Doublet2D <: Laplace
    str :: Float64
    x0 :: Float64
    y0 :: Float64 
end

potential(dub :: Doublet2D, x, y) = -dub.str / (2π) * (y - dub.y0) / ((x - dub.x0)^2 + (y - dub.y0)^2)

abstract type Panel end

struct Panel3D <: Panel
    p1 :: Tuple{Float64, Float64, Float64}
    p2 :: Tuple{Float64, Float64 ,Float64}
    p3 :: Tuple{Float64, Float64, Float64}
    p4 :: Tuple{Float64, Float64 ,Float64}
    # function Panel3D(p1, p2, p3, p4)  # Constructor needed for checking orientation
    #     if 
end

panel_dist(panel_1 :: Panel, panel_2 :: Panel) = norm(collocation_point(panel_2) .- collocation_point(panel_1))
collocation_point(panel :: Panel3D) = panel.p1 .+ panel.p2 .+ panel.p3 .+ panel.p4
panel_normal(panel :: Panel3D) = cross(panel.p2 .- panel.p1, panel.p3 .- panel.p2)

struct Line
    r1 :: Array{Float64, 2}
    r2 :: Array{Float64, 2}
end

velocity(line :: Line, r, Γ, ε = 1e-6) = let r1 = r .- line.r1, r2 = r .- line.r2, r1_x_r2 = r1 × r2;
    (r1 || r2 || r1_x_r2 ) < ε ? [0,0,0] : Γ/(4π) * r1_x_r2 / norm(r1_x_r2) * r1 .- r2 * (r1 / norm(r1) .- r2 / norm(r2)) end

mutable struct Horseshoe
    vortex_lines :: Array{Lines}
    strength :: Float64
    function Horseshoe(vortex_line :: Line, angle :: Float64, strength = 1, bound = 1e5)
        r1, r2 = vortex_line.r1, vortex_line.r2
        tline_1 = Line([ bound, r1[2], r1[3] * (sin ∘ deg2rad)(angle) ], r1)
        tline_2 = Line(r2, [ bound, r2[2], r2[3] * (sin ∘ deg2rad)(angle) ])

        new([vortex_line, tline_1, tline_2], strength)
    end
end

velocity(horseshoe :: Horseshoe) = sum(velocity.(horseshoe.vortex_lines, horseshoe.strength))
downwash_velocity(horseshoe :: Horseshoe) = let vels = velocity.(horseshoe.vortex_lines, horseshoe.strength); vels[1] + vels[3] end

#---------------------------------Matrix setup--------------------------------------#

influence_coefficient(panel_1 :: Panel3D, panel_2 :: Panel3D) = dot(velocity(panel_1, 1, collocation_point(panel_2), panel_normal(panel_2))

influence_matrix(panels :: Array{<: Panel}) = [ influence_coefficient(panel_j, panel_i) for panel_j in panels, for panel_i in panels ]

boundary_condition(panels :: Array{<: Panel}, uniform :: Uniform3D) = - [ dot(velocity(uniform), panel_normal(panel)) for panel in panels ]

solve_case(panels :: Array{<: Panel}, uniform :: Uniform3D) = influence_matrix(panels) \ boundary_condition(panels, uniform)

#-------------------------Force evaluations------------------------------------#

lift_coefficient(Γ, Δy, speed) = 2Γ * Δy / speed
# induced_drag_coefficient()


end