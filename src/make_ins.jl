function make_ins(
    elements::Array{String, 1},
    elements_to_check::Array{String, 1},
)
    
    n_element = length(elements)
    
    ins = Array{Int64, 1}(
        undef,
        n_element,
    )
    
    element_to_check_nothing = Dict(element=>nothing for element in elements_to_check)
    
    @inbounds @fastmath @simd for index in 1:n_element

        if haskey(
            element_to_check_nothing,
            elements[index],
        )
            
            ins[index] = 1

        else
            
            ins[index] = 0

        end

    end
    
    ins

end


function make_ins(
    element_index::Dict{String, Int64},
    elements_to_check::Array{String, 1},
)
    
    ins = fill(
        0,
        length(element_index),
    )
    
    @inbounds @fastmath @simd for element in elements_to_check

        index = get(
            element_index,
            element,
            nothing,
        )

        # TODO: Use the best practice to check for nothing
        if index !== nothing

            ins[index] = 1

        end
        
    end
    
    ins

end
