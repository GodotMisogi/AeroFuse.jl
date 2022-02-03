struct AircraftSurfaces <: AbstractAircraft
    surfs :: Vector{Aircraft}
end

wing_dictionary(c_r = 1.0, λ = 0.0, ι = 0.0, τ_r = 0.0, τ_t = 0.0, b = 1.0, δ = 5.0, Λ = 0.0, pos = zeros(3)) = 
    Dict(
        :root_chord => c_r,
        :taper      => λ,
        :incidence  => ι,
        :root_twist => τ_r,
        :tip_twist  => τ_t,
        :span       => b,
        :dihedral   => δ,
        :sweep      => Λ,
        :position   => pos,
        )

function VanillaAirplane(; surfs = Dict(:wing  => wing_dictionary(1.0, 0.6, 0.0, 2.0, 2.0, 10.0, 11.3, 2.29, [0., 0., 0.]),
                                        :htail => wing_dictionary(0.7, 0.6, 0.0, 0.0, 0.0, 2.5,  0.0, 6.39, [4., 0., 0.]),
                                        :vtail => wing_dictionary(0.7, 0.6, 0.0, 0.0, 0.0,  1.0,  0.0, 7.97, [4., 0., 0.])))

    wing  = WingSection(root_chord = surfs[:wing][:root_chord],
                        taper      = surfs[:wing][:taper],
                        root_twist = surfs[:wing][:root_twist],
                        tip_twist  = surfs[:wing][:tip_twist],
                        span       = surfs[:wing][:span],
                        dihedral   = surfs[:wing][:dihedral],
                        sweep        =  surfs[:wing][:sweep],
                        angle      = surfs[:wing][:incidence],
                        axis       = [0., 1., 0.])

    htail = WingSection(root_chord = surfs[:htail][:root_chord],
                        taper      = surfs[:htail][:taper],
                        root_twist = surfs[:htail][:root_twist],
                        tip_twist  = surfs[:htail][:tip_twist],
                        span       = surfs[:htail][:span],
                        dihedral   = surfs[:htail][:dihedral],
                        sweep        =  surfs[:htail][:sweep],
                        angle      = surfs[:htail][:incidence],
                        position   = surfs[:htail][:position],
                        axis       = [0., 1., 0.])

    vtail = HalfWingSection(root_chord = surfs[:vtail][:root_chord],
                            taper      = surfs[:vtail][:taper],
                            root_twist = surfs[:vtail][:root_twist],
                            tip_twist  = surfs[:vtail][:tip_twist],
                            span       = surfs[:vtail][:span],
                            dihedral   = surfs[:vtail][:dihedral],
                            sweep        =  surfs[:vtail][:sweep],
                            position   = [1., 0., 0.],
                            angle      = 90.,
                            axis       = [1., 0., 0.])

    AircraftSurfaces([wing, htail, vtail])
end