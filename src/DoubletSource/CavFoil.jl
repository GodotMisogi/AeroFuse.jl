module CavFoil

include("DoubletSource.jl")

using .DoubletSource

function cavity_split_foil(foil :: Foil, cavity_start, cavity_end)
    upper, lower = split_foil(foil.coords)

    xs = [ first(coords) for coords in upper ]
    c = abs(maximum(xs) - minimum(xs))
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

function solve_strengths(panels :: AbstractVector{Panel2D}, cavity_panels :: AbstractVector{Panel2D})

    # Cavity - Transition condition
    trans = zeros(num_cavity)
    for n in 1:num_cavity
        sfj += wetted_cavpanels[n].length
        scj += cavity_panels[n].length
        trans[n] = scj * (1 - cavPressure(sfj, sl))
    end 

    wetdub_wetdub = doublet_matrix(panels, panels)
    cavdub_wetdub = doublet_matrix(cavity_panels, panels)
    wetdub_cavsrc = doublet_matrix(panels, cavity_panels)
    cavdub_cavsrc = doublet_matrix(cavity_panels, cavity_panels)
    cavsrc_wetdub = source_matrix(cavity_panels, panels)
    cavsrc_cavsrc = source_matrix(cavity_panels, cavity_panels)
    kutta = kutta_condition(panels)

    cavdub_wetdub_trans = cavdub_wetdub * trans
    cavdub_cavsrc_trans = cavdub_cavsrc * trans

    # Cavity - Closure condition
    cavclose = zeros(num_cavity)
    for n in 1:num_cavity
        sfj += wetted_cavpanels[n].length
        cavclose[n] = cavity_length / (1 - cavPressure(sfj, sl))
    end

    AIC = [ wetdub_wetdub  cavsrc_wetdub  cavdub_wetdub_trans;
            wetdub_cavsrc  cavsrc_cavsrc  cavdub_cavsrc_trans;
                kutta'       cav_close'             0        ]
                            
    # Unkown doublet strength at leading edge of cavity
    AIC[num_wetted + 1:end - 2, leading_cavpanel] += sum(cavdub_cavsrc, dims = 2) ./ 2
    AIC[num_wetted + 1:end - 2, leading_cavpanel + 1] += sum(cavdub_cavsrc, dims = 2) ./ 2
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