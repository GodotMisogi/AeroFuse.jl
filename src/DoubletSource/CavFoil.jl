module CavFoil

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

end