using PlotlyJS
using WebIO

include("compute_set_enrichment.jl")
include("make_vector_01.jl")
include("sort_vectors.jl")


function plot_set_enrichment(
    element_values::Vector{Float64},
    elements::Vector{String},
    set_elements::Vector{String};
    height::Real = 500,
    width::Real = 800,
    element_values_line_width::Int64 = 2,
    element_values_line_color::String = "#ffb61e",
    set_elements_marker_size::Int64 = 16,
    set_elements_marker_line_width::Int64 = 2,
    set_elements_marker_color::String = "#000000",
    cumulative_sums_line_width::Real = 2,
    cumulative_sums_line_color::String = "#8db255",
    title_text::String = "Set Enrichment Plot",
    element_value_name = "Element<br>Value",
    yaxis_annotation_font_size = 16,
)
        
    n_element = length(element_values)

    annotation_template = Dict(
        "xref" => "paper",
        "yref" => "paper",
        "x" => -0.12,
        "font_size" => yaxis_annotation_font_size,
        "xanchor" => "center",
        "yanchor" => "middle",
        "showarrow" => false,
    )

    layout = Layout(
        height = height,
        width = width,
        legend_orientation = "h",
        legend_y = -0.2,
        title_text = "<b>$title_text</b>",
        xaxis1_zeroline = false,
        xaxis1_showspikes = true,
        xaxis1_spikemode = "across",
        xaxis1_spikedash = "solid",
        xaxis1_spikethickness = 1,
        xaxis1_automargin = true,
        xaxis1_title_text = "<b>Element Rank (n=$n_element)</b>",
        yaxis1_domain = (0, 0.3,),
        yaxis1_showline = true,
        yaxis1_automargin = true,
        yaxis2_domain = (0.3, 0.4,),
        yaxis2_showticklabels = false,
        yaxis2_showgrid = false,
        yaxis2_automargin = true,
        yaxis3_domain = (0.4, 0.9,),
        yaxis3_showline = true,
        yaxis3_automargin = true,
        margin_l = 130,
        annotations = [
            merge(
                annotation_template,
                Dict("y" => 0.15, "text" => "<b>$element_value_name</b>")
            ),
            merge(
                annotation_template,
                Dict("y" => 0.35, "text" => "<b>Set<br>Member</b>")
            ),
            merge(
                annotation_template,
                Dict("y" => 0.65, "text" => "<b>Set<br>Enrichment</b>")
            ),
        ],
    )

    x = 1:n_element

    element_values, elements = sort_vectors(
        [element_values, elements,],
        reverse = true,
    )

    element_values_trace = scatter(
        name = "Element Value",
        x = x,
        y = element_values,
        text = elements,
        line_width = element_values_line_width,
        line_color = element_values_line_color,
        fill = "tozeroy",
    )
    
    set_elements_01 = make_vector_01(elements, set_elements,)

    set_elements_bit = BitVector(set_elements_01)

    set_elements_trace = scatter(
        name = "Set Element",
        yaxis = "y2",
        x = x[set_elements_bit],
        y = zeros(sum(set_elements_01)),
        text = elements[set_elements_bit],
        mode = "markers",
        marker_symbol = "line-ns-open",
        marker_size = set_elements_marker_size,
        marker_line_width = set_elements_marker_line_width,
        marker_color = set_elements_marker_color,
        hoverinfo = "name+x+text",
    )
    
    set_enrichments, ks, auc = compute_set_enrichment(
        element_values,
        elements,
        set_elements_01;
        compute_cumulative_sums = true,
    )

    set_enrichments_trace = scatter(
        name = "Set Enrichment",
        yaxis = "y3",
        x = x,
        y = set_enrichments,
        text = elements,
        line_width = cumulative_sums_line_width,
        line_color = cumulative_sums_line_color,
        fill = "tozeroy",
    )

    plot(
        [element_values_trace, set_elements_trace, set_enrichments_trace,],
        layout,
    )

end
