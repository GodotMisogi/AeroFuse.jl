## Printing
#==========================================================================================#

"""
    print_info(wing :: AbstractWing)

Print the relevant geometric characteristics of an `AbstractWing`.
"""
function print_info(wing :: AbstractWing, head = "")
    labels = [ "Aspect Ratio", "Span (m)", "Area (m²)", "Mean Aerodynamic Chord (m)", "Mean Aerodynamic Center (m)" ]
    wing_info = properties(wing)
    data = Any[ labels wing_info ]
    header = [ head, "Value" ]
    h1 = TextHighlighter( (data,i,j) -> (j == 1), foreground = :cyan, bold = true)

    pretty_table(data; column_labels = header, alignment = [:c, :c], highlighters = [h1], 
        table_format = TextTableFormat(; @text__no_vertical_lines), formatters = [fmt__round(8)])
end