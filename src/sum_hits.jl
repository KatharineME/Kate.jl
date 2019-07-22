function sum_hits(hits::Array{Int64, 1})
    
    sum_ = 0
    
    @inbounds @fastmath @simd for i in 1:length(hits)
        
         sum_ += hits[i]
        
    end
    
    sum_
    
end 
