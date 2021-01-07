function kutta_condition!(AIC)
    AIC[end,1] = 1
    AIC[end,2] = -1
    AIC[end,end-2] = 1
    AIC[end,end-1] = -1
end

function matrix_assembly!(AIC, boco, panels :: AbstractVector{<: Panel2D}, woke_panel :: Panel2D, u)
    for i in eachindex(panels)
        for j in eachindex(panels)
            if i == j 
                AIC[i,j] = 0.5
                boco[i] -= source_influence(panels[j], panels[i]) * dot(u, panel_normal(panels[j]))
            else
                AIC[i,j] = doublet_influence(panels[j], panels[i])
                boco[i] -= source_influence(panels[j], panels[i]) * dot(u, panel_normal(panels[j]))
            end
        end
        AIC[i,end] = doublet_influence(woke_panel, panels[i])
    end

    kutta_condition!(AIC)

    AIC
end