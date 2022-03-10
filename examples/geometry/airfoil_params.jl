##
using AeroMDAO

## Test function
function doublet_source_case(foil, uniform)
    @time system = solve_case(foil,
                              uniform;
                              num_panels = 80)
                            
    cl = lift_coefficient(system)
    cls, cms, _ = surface_coefficients(system)

    println("Cl: $cl")
    println("Σᵢ Clᵢ: $(sum(cls))")
    println("Σᵢ Cmᵢ: $(sum(cms))")
end

uniform = Uniform2D(angle = 5)

## Foil processing
naca_foil = naca4((2,4,1,2)) # NACA 4-digit

doublet_source_case(naca_foil, uniform)

## Cosine spacing
cos_foil = cosine_spacing(naca_foil, 81)

doublet_source_case(cos_foil, uniform)

## Camber-thickness transformations
xcamthick = camber_thickness(cos_foil)
foiler = Foil(camber_thickness_to_coordinates(xcamthick[:,1], xcamthick[:,2], xcamthick[:,3]))

doublet_source_case(foiler, uniform)

## Fitting to CST
num_dv = 8

# Coordinates fitting
up, low  = split_surface(naca_foil)
alpha_u  = coordinates_to_CST(up, num_dv)
alpha_l  = coordinates_to_CST(low, num_dv)
cst_foil = kulfan_CST(alpha_u, alpha_l, (0., 0.), (0., 0.))

doublet_source_case(cst_foil, uniform)

## Camber-thickness fitting
alpha_c, alpha_t = camber_thickness_to_CST(cos_foil, num_dv)
cam_foil = camber_CST(alpha_c, alpha_t, 0., 0)

doublet_source_case(cam_foil, uniform)

## Control surface deflection
con_foil = control_surface(cos_foil; angle = 8., hinge = 0.75)

doublet_source_case(con_foil, uniform)

## Plotting library
using Plots
gr(dpi = 300);

## Cosine and camber-thickness
plot(aspect_ratio = 1)
plot!(naca_foil.x, naca_foil.y, label = "Original")
plot!(cos_foil.x, cos_foil.y, label = "Cosine")
plot!(xcamthick[:,1], xcamthick[:,2], ls = :dot, label = "Camber")
plot!(xcamthick[:,1], xcamthick[:,3], ls = :dot, label = "Thickness")
plot!(foiler.x, foiler.y, label = "Inverse Camber-Thickness")
plot!(cst_foil.x, cst_foil.y, label = "CST Coordinates Fit")
plot!(cam_foil.x, cam_foil.y, label = "CST Camber-Thickness Fit")
plot!(con_foil.x, con_foil.y, ls = :dash, label = "Control Surface Deflection")