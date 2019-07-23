function sum_hit_scores(
    scores::Array{Float64, 1},
    hits::Array{Int64, 1},
)
    
    sum_ = 0.0
    
    @inbounds @fastmath @simd for i in 1:length(scores)
        
        if hits[i] == 1
        
            sum_ += scores[i]
            
        end
        
    end
    
    sum_
    
end 
