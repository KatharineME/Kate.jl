function sum_hits(hits::Array{Int64, 1})
    
    sum_ = 0
    
    @inbounds @fastmath @simd for index in 1:length(hits)
        
         sum_ += hits[index]
        
    end
    
    sum_
    
end 
