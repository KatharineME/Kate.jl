function make_ins(elements::Vector{String}, elements_to_check::Vector{String},)
    
    n_element = length(elements)
    
    ins = Vector{Int64}(undef, n_element,)
    
    element_to_check_nothing = Dict(element => nothing for element in elements_to_check)
    
    @inbounds @fastmath @simd for index = 1:n_element

        if haskey(element_to_check_nothing, elements[index],)
            
            ins[index] = 1

        else
            
            ins[index] = 0

        end

    end
    
    return ins

end


function make_ins(
    element_index::Dict{String,Int64,},
    elements_to_check::Vector{String},
)
    
    ins = fill(0, length(element_index),)
    
    @inbounds @fastmath @simd for element in elements_to_check

        index = get(element_index, element, nothing,)

        if index !== nothing

            ins[index] = 1

        end
        
    end
    
    return ins

end
