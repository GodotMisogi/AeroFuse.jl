## Mesh convergence tests
using AeroMDAO
using Plots
using DataFrames

## Wing
wing = Wing(foils     = fill(naca4((0,0,1,2)), 2),
            chords    = [1.0, 0.6],
            twists    = [0.0, 0.0],
            spans     = [5.0],
            dihedrals = [5.],
            sweeps      = [5.]);
print_info(wing, "Wing")
S, b, c = projected_area(wing), span(wing), mean_aerodynamic_chord(wing);

## Assembly
ρ       = 1.225
ref     = [0.25, 0., 0.]
V, α, β = 1.0, 5.0, 5.0
Ω       = [0.0, 0.0, 0.0]
fs      = Freestream(V, α, β, Ω)

## Panelling error closure
panel_case(x, y) = solve_case(wing, fs;
                              rho_ref   = ρ,
                              r_ref     = ref,
                              area_ref  = S,
                              span_ref  = b,
                              chord_ref = c,
                              span_num  = y,
                              chord_num = x,
                              viscous   = false,
                              x_tr      = [0.3, 0.3],
                              spacing   = Uniform()
                             );

## Distributions
xs = 10:2:50
ys = xs .* 2# fill(maximum(xs) * 4, length(xs))

num = length(xs)

coeff_data = zeros(num, 9 + 2)
for (i, (x, y)) in enumerate(zip(xs, ys))
    @time nf_coeffs, ff_coeffs, CFs, CMs, horseshoe_panels, camber_panels, horseshoes, Γs =
    panel_case(x, y)

    println("Distribution: nx = $x, ny = $y")
    print_coefficients(nf_coeffs, ff_coeffs)

    coeff_data[i,1]     = x
    coeff_data[i,2]     = y
    coeff_data[i,3:5]   = ff_coeffs
    coeff_data[i,6:end] = nf_coeffs
end

##
data = DataFrame(coeff_data, :auto)
rename!(data, [:chord_num, :span_num, :CD_ff, :CY_ff, :CL_ff, :CD_nf, :CY_nf, :CL_nf, :Cl, :Cm, :Cn])

##
errors = @. abs((data[3:end,3:end] - data[2:end-1,3:end]) / data[2:end-1,3:end])

##
plot(xlabel = "Chordwise Panels", ylabel = "Error", yscale = :log10)

# CDs Errors
plot!(data[2:end-1,"chord_num"], errors[!,"CD_nf"], label = "CD Nearfield")
plot!(data[2:end-1,"chord_num"], errors[!,"CD_ff"], label = "CD Farfield")

# CLs Errors
plot!(data[2:end-1,"chord_num"], errors[!,"CL_nf"], label = "CL Nearfield")
plot!(data[2:end-1,"chord_num"], errors[!,"CL_ff"], label = "CL Farfield")

# CYs Errors
plot!(data[2:end-1,"chord_num"], errors[!,"CY_nf"], label = "CY Nearfield")
plot!(data[2:end-1,"chord_num"], errors[!,"CY_ff"], label = "CY Farfield")

##
CL = plot(data[2:end,"chord_num"], data[2:end,"CL_nf"], marker = :dot, label = "Nearfield", ylabel = "CL")
plot!(data[2:end,"chord_num"], data[2:end,"CL_ff"], marker = :dot, label = "Farfield")

CD = plot(data[2:end,"chord_num"], data[2:end,"CD_nf"], marker = :dot, label = "Nearfield", ylabel = "CD")
plot!(data[2:end,"chord_num"], data[2:end,"CD_ff"], marker = :dot, label = "Farfield")

CY = plot(data[2:end,"chord_num"], data[2:end,"CY_nf"], marker = :dot, label = "Nearfield", ylabel = "CY")
plot!(data[2:end,"chord_num"], data[2:end,"CY_ff"], marker = :dot, label = "Farfield")

Cl = plot(data[2:end,"chord_num"], data[2:end,"Cl"], marker = :dot, label = "Nearfield", ylabel = "Cl")
Cm = plot(data[2:end,"chord_num"], data[2:end,"Cm"], marker = :dot, label = "Nearfield", ylabel = "Cm")
Cn = plot(data[2:end,"chord_num"], data[2:end,"Cn"], marker = :dot, label = "Nearfield", ylabel = "Cn")


plot(CL, CD, CY, Cl, Cm, Cn, layout = (2,3), size = (1280, 720))
