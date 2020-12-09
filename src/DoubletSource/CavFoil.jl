module CavFoil

include("DoubletSource.jl")

using .DoubletSource

function cavity_split_foil(foil :: Foil, cavity_start, cavity_end)
    upper, lower = split_foil(foil.coords)

    xs = [ first(coords) for coords in upper ]
    c = (abs ∘ map)(-, extrema(xs))
    cavity_le = cavity_start
    l = (cavity_end - cavity_start) / 100 / c

    leading_cav_upper = filter(x -> first(x) <= cavity_start * c / 100, upper)
    wetted_cav_upper = filter(x -> first(x) >= cavity_start * c / 100 && first(x) <= cavity_end * c / 100, upper)
    cav_trailing_upper = filter(x -> first(x) >= cavity_end * c / 100, upper)
    leading_cav_lower = filter(x -> first(x) <= cavity_start * c / 100, lower)
    wetted_cav_lower = filter(x -> first(x) >= cavity_start * c / 100 && first(x) <= cavity_end * c / 100, lower)
    cav_trailing_lower = filter(x -> first(x) >= cavity_end * c / 100, lower)

    return (leading_cav_upper, wetted_cav_upper, cav_trailing_upper, leading_cav_lower, wetted_cav_lower, cav_trailing_lower)
end

function cavity_length_foil(foil :: Foil, cavity_start, cavity_end, N_panels = 40, N_cavpanels = 10)

    leading_cav_upper, wetted_cav_upper, cav_trailing_upper, leading_cav_lower, wetted_cav_lower, cav_trailing_lower = cavity_split_foil(foil, cavity_start, cavity_end)

    wet_foil = cosine_foil(wetted_cav_lower, N_cavpanels)

    if cavity_le == 0
        # Upper: Leading (cavity) to trailing
        wu = cosine_interp(wetted_cav_upper, Int(trunc(N_cavpanels)))
        tu = cosine_interp(cav_trailing_upper, Int(trunc((N_panels - N_cavpanels)/3)))
        # Lower: Leading (cavity) to trailing
        wl = cosine_interp(wetted_cav_lower, Int(trunc((N_panels - N_cavpanels)/3)))
        tl = cosine_interp(cav_trailing_lower, Int(trunc((N_panels - N_cavpanels)/3)))
        # Concatenate coordinates and reorder
        wet_upper = [wu; tu]
        wet_lower = [wl; tl]
    else
        # Upper: Leading to trailing
        lu = cosine_interp(leading_cav_upper, Int(trunc((N_panels - N_cavpanels)/4)))
        wu = cosine_interp(wetted_cav_upper, Int(trunc(N_cavpanels)))
        tu = cosine_interp(cav_trailing_upper, Int(trunc((N_panels - N_cavpanels)/4)))
        # Lower: Leading to trailing
        ll = cosine_interp(leading_cav_lower, Int(trunc((N_panels - N_cavpanels)/4)))
        wl = cosine_interp(wetted_cav_lower, Int(trunc(N_cavpanels)))
        tl = cosine_interp(cav_trailing_lower, Int(trunc((N_panels - N_cavpanels)/4)))
        # Concatenate coordinates and reorder
        wet_upper = [ tu[1:end-1]; wu; lu[2:end] ]
        wet_lower = [ ll[1:end-1]; wl; tl[2:end] ]
    end
end

# Cavity termination models 
cavity_pressure(sf, sl, A = 0.5, ν = 1.0, λ = 0.1) = sf < (1 - λ) * sl ? 0 : A * ((sf - (1 - λ) * sl) / (sl - (1 - λ) * sl))^ν

function cavity_influence_matrix(panels :: AbstractVector{Panel2D}, cavity_panels :: AbstractVector{Panel2D})

    wetted_panels, wetted_cavpanels = cavity_split_panels(panels)

    wetdub_wetdub = doublet_matrix(wetted_panels, wetted_panels)
    cavdub_wetdub = doublet_matrix(cavity_panels, wetted_panels)
    wetdub_cavsrc = doublet_matrix(wetted_panels, cavity_panels)
    cavdub_cavsrc = doublet_matrix(cavity_panels, cavity_panels)
    cavsrc_wetdub = source_matrix(cavity_panels, wetted_panels)
    cavsrc_cavsrc = source_matrix(cavity_panels, cavity_panels)
    kutta         = kutta_condition(wetted_panels)
    wake_vec      = wake_vector(wetted_panels)

    sl = sum(panel_length.(wetted_cavpanels))

    trans = cumsum(panel_length.(cavity_panels)) .* (1 .- cavity_pressure.(cumsum(panel_length.(wetted_cavpanels)), sl))

    cavdub_wetdub_trans = cavdub_wetdub * trans
    cavdub_cavsrc_trans = cavdub_cavsrc * trans

    cavclose = cavity_length ./ (1 - cavity_pressure.(cumsum(panel_length.(wetted_cavpanels)), sl))

    AIC = [ wetdub_wetdub    cavsrc_wetdub    zeros(num_wet)         wake_vec      ;
            wetdub_cavsrc    cavsrc_cavsrc  cavdub_cavsrc_trans cavdub_wetdub_trans;
           zeros(num_wet)'     cav_close'            0                  0          ;
              kutta'        zeros(num_cav)'          0                  0          ]
                            
    # Unkown doublet strength at leading edge of cavity
    AIC[num_wetted + 1:end - 2, leading_cavpanel] += sum(cavdub_cavsrc, dims = 2) ./ 2
    AIC[num_wetted + 1:end - 2, leading_cavpanel + 1] += sum(cavdub_cavsrc, dims = 2) ./ 2
end

function cavity_boundary_vector(panels :: AbstractVector{Panel2D}, cavity_panels :: AbstractVector{Panel2D}, freestream :: Uniform2D)
    cavdub_wetdub = doublet_matrix(cavity_panels, panels)
    wetdub_cavsrc = doublet_matrix(panels, cavity_panels)
    
    # Cavity doublets boundary condition
    cavdub_bound = [ potential(freestream, collocation_point(cav_panel)...) - potential(freestream, collocation_point(first(cavity_panels))...) for cav_panel in cavity_panels ]

    # Influence of wetted sources on wetted doublets
    wetsrc_wetdub = source_matrix(panels, cavity_panels)

    # Wetted sources boundary condition (∂Φ/∂n|_wet = - U · n)
    wetsrc_bound = -source_strengths(panels, freestream)

    # ----------- Cavity sources ------------#

    # Influence of cavity doublets on cavity sources

    # Influence of wetted sources on cavity sources
    wetsrc_cavsrc = source_matrix(cavity_panels, panel)

    # ------- Cavity closure condition -------#

    cav_closure = -sum(cavclose .* -source_strengths(cavity_panels, freestream))

    RHS = [ cavdub_wetdub * cavdub_bound + wetsrc_wetdub * wetsrc_bound;
            cavdub_cavsrc * cavdub_bound + wetsrc_cavsrc * wetsrc_bound;
                                  cav_closure                          ;
                                       0                               ]

end    


function solve_strengths(panels :: AbstractVector{Panel2D}, cavity_panels :: AbstractVector{Panel2D}, freestream :: Uniform2D)
    AIC = cavity_influence_matrix(panels, cavity_panels)
    RHS = cavity_boundary_vector(panels, cavity_panels, freestream)
    
    # Solve system
    strengths = AIC \ RHS

    φs = strengths[1:num_wetted], strengths[num_wetted+1:end-2]
    # Cavity velocity, cavitation number and wake panel strength
    qc = strengths[end - 1]
    woke_panel.doublet_strength = strengths[end]
    σ = qc^2 / freestream.mag^2 - 1
    println("Cavity Velocity: ", qc, ", Cavitation Number: ", σ)
end

function cavity_height(cavity_panels :: AbstractVector{Panel2D})
    heights = zeros(num_cavity)
    cavity_length = sum([ panel.length for panel in cavity_panels ])
    for n in 1:num_cavity
        sfj += wetted_cavpanels[n].length
        heights[n] = cavity_length * (sum(u .* cavity_panels[n].normal) + cavity_panels[n].source_strength) / (qc * (1 - cavPressure(sfj, sl)))
    end

    for (h1, h2, p1, p2) in zip(heights[1:end - 1], heights[2:end], cavity_panels[1:end - 1], cavity_panels[2:end])
        # Update coordinates
        angle = p1.normal .+ p2.normal
        p1.xe = p2.xs += (h1 + h2) / 2 * angle[1]
        p1.ye = p2.ys += (h1 + h2) / 2 * angle[2]

        # Panel 1 update
        p1.length = mag([p1.xe - p1.xs, p1.ye - p1.ys])
        p1.angle = atan(p1.ye - p1.ys, p1.xe - p1.xs)
        p1.normal = (-sin(p1.angle), cos(p1.angle))
        p1.tangent = (cos(p1.angle), sin(p1.angle))
        p1.xc, p1.yc = (p1.xe + p1.xs) / 2.0, (p1.ye + p1.ys) / 2.0

        # Panel 2 update
        p2.length = mag([p2.xe - p2.xs, p2.ye - p2.ys])
        p2.angle = atan(p2.ye - p2.ys, p2.xe - p2.xs)
        p2.normal = (-sin(p2.angle), cos(p2.angle))
        p2.tangent = (cos(p2.angle), sin(p2.angle))
        p2.xc, p2.yc = (p2.xe + p2.xs) / 2.0, (p2.ye + p2.ys) / 2.0
    end

    return heights
end

end