struct AircraftSurfaces <: Aircraft
    surfs :: Vector{Aircraft}
end

wing_dictionary(c_r = 1.0, λ = 0.0, ι = 0.0, τ_r = 0.0, τ_t = 0.0, b = 1.0) = 
    Dict(
        :root_chord => c_r,
        :taper      => λ,
        :incidence  => ι,
        :root_twist => τ_r,
        :tip_twist  => τ_t,
        :span       => b,
        )

function VanillaAirplane(; surfs = Dict(:wing  => wing_dictionary(1.0, 0.6, 0.0, 0.0, 0.0, 5.0),
                                        :htail => wing_dictionary(0.6, 0.4, 0.0, 0.0, 0.0, 1.2),
                                        :vtail => wing_dictionary(0.6, 0.4, 0.0, 0.0, 0.0, 0.6)))

    wing  = WingSection(root_chord = surfs[:wing][:root_chord],
                        taper      = surfs[:wing][:taper],
                        span       = surfs[:wing][:span],
                        angle      = surfs[:wing][:incidence],
                        axis       = [0., 1., 0.])

    htail = WingSection(root_chord = surfs[:htail][:root_chord],
                        taper      = surfs[:htail][:taper],
                        span       = surfs[:htail][:span],
                        angle      = surfs[:htail][:incidence],
                        axis       = [0., 1., 0.])

    vtail = WingSection(root_chord = surfs[:vtail][:root_chord],
                        taper      = surfs[:vtail][:taper],
                        span       = surfs[:vtail][:span],
                        position   = [1., 0., 0.],
                        angle      = 90.,
                        axis       = [1., 0., 0.])

    AircraftSurfaces([wing, htail, vtail])
end