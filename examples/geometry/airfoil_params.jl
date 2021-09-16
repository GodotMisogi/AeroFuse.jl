##
using AeroMDAO

## Test function
function doublet_source_case(coords, uniform)
    @time cl, cls, cms, cps, panels = solve_case(Foil(coords),
                                             uniform;
                                             num_panels = 80)

    println("Cl: $cl")
    println("Σᵢ Clᵢ: $(sum(cls))")
    println("Σᵢ Cmᵢ: $(sum(cms))")
end

uniform = Uniform2D(1., 5.)

## Foil processing
coords = naca4(2,4,1,2) # NACA 4-digit

doublet_source_case(coords, uniform)

## Cosine spacing
cos_foil = cosine_foil(coords, 60)

doublet_source_case(cos_foil, uniform)

## Camber-thickness transformations
xcamthick = foil_camthick(cos_foil)
foiler = camthick_foil(xcamthick[:,1], xcamthick[:,2], xcamthick[:,3])

doublet_source_case(foiler, uniform)

## Fitting to CST
num_dv = 8

# Coordinates fitting
up, low  = split_foil(coords)
alpha_u  = coords_to_CST(up, num_dv)
alpha_l  = coords_to_CST(low, num_dv)
cst_foil = kulfan_CST(alpha_u, alpha_l, (0., 0.), 0.0)

doublet_source_case(cst_foil, uniform)

## Camber-thickness fitting
alphas   = camthick_to_CST(cos_foil, num_dv)
cam_foil = camber_CST(alphas..., (0., 0.), 0)

doublet_source_case(cam_foil, uniform)

## Plotting library
using Plots
plotly(dpi = 300);

## Cosine and camber-thickness
plot(aspect_ratio = 1)
plot!(coords[:,1], coords[:,2],       label = "Original")
plot!(cos_foil[:,1], cos_foil[:,2],   label = "Cosine")
plot!(xcamthick[:,1], xcamthick[:,2], label = "Camber")
plot!(xcamthick[:,1], xcamthick[:,3], label = "Thickness")
plot!(cst_foil[:,1], cst_foil[:,2],   label = "CST Coordinates Fit")
plot!(cam_foil[:,1], cam_foil[:,2],   label = "CST Camber-Thickness Fit")