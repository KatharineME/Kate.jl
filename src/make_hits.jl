function make_hits(
    elements::Array{String, 1},
    elements_to_find::Array{String, 1},
)
    
    n_element = length(elements)
    
    hits = Array{Int64, 1}(
        undef,
        n_element,
    )
    
    element_to_find_nothing = Dict(element=>nothing for element in elements_to_find)
    
    @inbounds @fastmath @simd for index in 1:n_element

        if haskey(
            element_to_find_nothing,
            elements[index],
        )
            
            hits[index] = 1

        else
            
            hits[index] = 0

        end

    end
    
    hits

end


function make_hits(
    element_index::Dict{String, Int64},
    elements_to_find::Array{String, 1},
)
    
    hits = fill(
        0,
        length(element_index),
    )
    
    @inbounds @fastmath @simd for element in elements_to_find

        index = get(
            element_index,
            element,
            nothing,
        )

        if index !== nothing

            hits[index] = 1

        end
        
    end
    
    hits

end
