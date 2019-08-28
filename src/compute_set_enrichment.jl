using DataFrames
using Distributed

include("make_ins.jl")
include("sum_vector.jl")

function compute_set_enrichment(
    element_values::Vector{Float64},
    elements::Vector{String},
    set_elements::Vector{String};
    sort_element_values::Bool = true,
    element_index::Union{Nothing,Dict{String,Int64,}} = nothing,
    compute_cumulative_sum::Bool = true,
)
    
    if sort_element_values

        element_values_sort_indices = sortperm(element_values)

        element_values = element_values[element_values_sort_indices]

        elements = elements[element_values_sort_indices]

    end

    if element_index === nothing
        
        ins = make_ins(elements, set_elements,)
    
    else
        
        ins = make_ins(element_index, set_elements,)
        
    end

    element_values_abs = abs.(element_values)

    element_values_abs_ins_sum = sum_values(element_values_abs, ins,)

    n_element = length(elements)
    
    d_down = -1 / (n_element - sum_values(ins))
    
    value = 0.0

    if compute_cumulative_sum

        cumulative_sum = Vector{Float64}(undef, n_element,)

    else

        cumulative_sum = nothing

    end
    
    ks = 0.0

    ks_abs = 0.0

    auc = 0.0
    
    @inbounds @fastmath @simd for index = 1:n_element
        
        if ins[index] == 1
            
            value += element_values_abs[index] / element_values_abs_ins_sum
            
        else
            
            value += d_down
            
        end
        
        if compute_cumulative_sum

            cumulative_sum[index] = value

        end

        value_abs = abs(value)
        
        if ks_abs < value_abs

            ks = value

            ks_abs = value_abs

        end

        auc += value
            
    end
    
    return cumulative_sum, ks, auc
    
end


function compute_set_enrichment(
    element_values::Vector{Float64},
    elements::AbstractVector{String},
    set_elements::Dict{String,Vector{String},};
    sort_element_values::Bool = true,
)

    if sort_element_values

        element_values_sort_indices = sortperm(element_values)

        element_values = element_values[element_values_sort_indices]

        elements = elements[element_values_sort_indices]

    end

    if length(set_elements) < 10
        
        element_index = nothing
        
    else
        
        element_index = Dict(element => index for (element, index,) in zip(
            elements,
            1:length(elements),
        ))
        
    end

    set_enrichment = Dict{
        String,
        Tuple{Union{Vector{Nothing,Float64},},Float64,Float64,},
    }()

    for (set, set_elements_,) in set_elements

        set_enrichment[set] = compute_set_enrichment(
            element_values,
            elements,
            set_elements_;
            sort_element_values = false,
            element_index = element_index,
            compute_cumulative_sum = false,
        )

    end

    return set_enrichment
    
end


function compute_set_enrichment(
    element_x_sample::DataFrame,
    set_elements::Dict{String,Vector{String},},
    statistic::String,
)
    
    elements = element_x_sample[:, Symbol("Gene")]

    set_x_sample = DataFrame(Symbol("Gene Set") => sort(collect(keys(set_elements))))

    if statistic == "ks"

        set_enrichment_score_index = 2

    elseif statistic == "auc"

        set_enrichment_score_index = 3

    elseif statistic == "js"

        set_enrichment_score_index = 4

    end

    for sample in names(element_x_sample)[2:end]

        has_element_value = findall(!ismissing, element_x_sample[:, sample],)

        set_enrichment = compute_set_enrichment(
            Float64.(element_x_sample[has_element_value, sample]),
            elements[has_element_value],
            set_elements,
        )

        set_x_sample

        set_x_sample[!, sample] = collect(set_enrichment[set][set_enrichment_score_index] for set in set_x_sample[:, Symbol("Gene Set")])

    end

    return set_x_sample
    
end
