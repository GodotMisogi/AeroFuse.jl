"""
Generates a wake Panel3D aligned with a normalised direction vector.
"""
function wake_panel(panel :: Panel3D, freestream :: Freestream, bound = 1.)
    wake = bound .* (normalize ∘ velocity)(freestream)
    Panel3D(panel.p2, panel.p2 .+ wake, panel.p3 .+ wake, panel.p3)
end

"""
Generates wake Panel3Ds given the Panel3Ds of a wing corresponding to the trailing edge and a Freestream. The wake length and number of wake panels must be provided.
"""
make_wake(last_panels :: Array{Panel3D}, freestream :: Freestream, wake_length, wake_num) = (permutedims ∘ accumap)(panel -> wake_panel(panel, freestream, wake_length / wake_num), wake_num, last_panels)