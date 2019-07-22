function make_hits(
    elements::Array{String, 1},
    elements_to_find::Array{String, 1},
)
    
    n = length(elements)
    
    hits = Array{Int64, 1}(
        undef,
        n,
    )
    
    elements_to_find_ = Dict(e=>nothing for e in elements_to_find)
    
    @inbounds @fastmath @simd for i in 1:n

        if haskey(
            elements_to_find_,
            elements[i],
        )
            
            hit = 1

        else
            
            hit = 0

        end
        
        hits[i] = hit

    end
    
    hits

end
