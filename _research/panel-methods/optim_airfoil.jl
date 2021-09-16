##
using Revise
using AeroMDAO
using Optim
using Plots
plotlyjs()

## Surface optimization
#============================================#

# Objective function
function optimize_CST(αs, α, n_upper :: Integer, le = 0.)
    airfoil = (Foil ∘ kulfan_CST)(αs[1:n_upper], αs[n_upper+1:end], (0., 0.), le, 80)
    uniform = Uniform2D(1.0, α)
    cl = solve_case(airfoil, uniform, 60)
end

## Test
dzs = (0., 0.)
α_u0 = [0.1, 0.1, 0.1, 0.1, 0.1]
α_l0 = [-0.2, -0.1, -0.1, -0.1]
α_ul0 = [α_u0; α_l0]

n_upper = length(α_u0)  # Number of upper surface variables
le = 0.0                # Leading edge modification
α = 0.                  # Angle of attack
cl0 = 1.2               # Target lift coefficient
CST_cl = @time optimize_CST(α_ul0, α, n_upper, le)
println("CST Cl: $CST_cl")

## Plot
CST_test = kulfan_CST(α_u0, α_l0, dzs, le, 80)
plot(CST_test[:,1], CST_test[:,2],
     label = "Initial Surface", aspectratio = 1)

## Optimization
l_bound = [ fill(1e-12, length(α_u0)); fill(-Inf, length(α_l0)) ]
u_bound = [ fill(Inf, length(α_u0)); fill(5e-1, length(α_l0)) ]

resi_CST = optimize(x -> abs(optimize_CST(x, α, n_upper, le) - cl0),
                    l_bound, u_bound,               # Bounds
                    α_ul0,                          # Initial value
                    Fminbox(GradientDescent()),
                    autodiff = :forward,
                    Optim.Options(
                                #   extended_trace = true,
                                  show_trace = true
                                 )
                    )


## Plot
CST_opt = kulfan_CST(resi_CST.minimizer[1:n_upper], resi_CST.minimizer[n_upper+1:end], dzs, le, 80)

fig = plot(aspectratio = 1, dpi = 300)
plot!(CST_test[:,1], CST_test[:,2],
        label = "Initial Surface")
plot!(CST_opt[:,1], CST_opt[:,2],
        label = "Optimal Surface")
# savefig(fig, "_research/tmp/CST_surface.pdf")

## Optimal test
CST_cl = optimize_CST(resi_CST.minimizer, α, n_upper)
println("CST Cl: $CST_cl")

## Camber-thickness optimization
#============================================#

# Camber optimization
function optimize_camber_CST(αs, α, num_cam)
    airfoil = (Foil ∘ camber_CST)(αs[1:num_cam], αs[num_cam+1:end], (0., 0.), 0., 80)
    uniform = Uniform2D(1.0, α)
    cl = solve_case(airfoil, uniform, 80)
end

## Test
α_c0 = zeros(6)
α_t0 = fill(0.4, 6)
α_ct0 = [α_c0; α_t0]

num_cam = length(α_c0)  # Number of camber variables
α = 0.                  # Angle of attack
cl0 = 1.2               # Target lift coefficient
camber_cl = @time optimize_camber_CST(α_ct0, α, num_cam)
println("Camber Cl: $camber_cl")

## Plot
camber_test = camber_CST(α_ct0[1:num_cam], α_ct0[num_cam+1:end], (0., 0.), 0., 80)
camthick_test = foil_camthick(camber_test)
plot(camber_test[:,1], camber_test[:,2],
     label = "Initial Surface", aspectratio = 1)
# plot(camthick_test[:,1], camthick_test[:,2],
#                           label = "Initial Camber", aspectratio = 1)
# plot!(camthick_test[:,1], camthick_test[:,3],
#                           label = "Initial Thickness", aspectratio = 1)

## Optimization
l_bound = [ fill(-1e-12, length(α_c0)); [0.1, 0.12, 0.12, 0.05, 0.05, 0.05] ]
u_bound = fill(Inf, length(α_ct0))

resi_cam = optimize(x -> abs(optimize_camber_CST(x, α, num_cam) - cl0),
                    l_bound, u_bound,               # Bounds
                    α_ct0,                          # Initial value
                    Fminbox(GradientDescent()),
                    autodiff = :forward,
                    Optim.Options(
                                #   extended_trace = true,
                                  show_trace = true
                                 )
                    )

## Plot
camber_opt = camber_CST(resi_cam.minimizer[1:num_cam], resi_cam.minimizer[num_cam+1:end], (0., 0.), 0., 80)

surf_fig = plot(aspectratio = 1, dpi = 300)
plot!(camber_test[:,1], camber_test[:,2],
      label = "Initial Surface", aspectratio = 1)
plot!(camber_opt[:,1], camber_opt[:,2],
      label = "Optimal Surface for α = $α, Cl = $cl0")
savefig(surf_fig, "_research/tmp/surface.pdf")

##
camthick_opt = foil_camthick(camber_opt)
camthick_fig = plot(aspectratio = 1, dpi = 300)
plot!(camthick_test[:,1], camthick_test[:,2],
      label = "Initial Camber", aspectratio = 1)
plot!(camthick_test[:,1], camthick_test[:,3],
      label = "Initial Thickness", aspectratio = 1)
plot!(camthick_opt[:,1], camthick_opt[:,2],
      label = "Optimal Camber for α = $α, Cl = $cl0", aspectratio = 1)
plot!(camthick_opt[:,1], camthick_opt[:,3],
      label = "Optimal Thickness for α = $α, Cl = $cl0", aspectratio = 1)
savefig(camthick_fig, "_research/tmp/camthick.pdf")

## Optimal test
camber_cl = optimize_camber_CST(resi_cam.minimizer, α, num_cam)
println("Camber Cl: $camber_cl")