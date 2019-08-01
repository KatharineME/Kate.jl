function sum_values(values::Array{
    Int64,
    1,
})::Int64
    
    sum_ = 0
    
    @inbounds @fastmath @simd for index in 1:length(values)
        
        sum_ += values[index]
        
    end
    
    sum_
    
end


function sum_values(
    values::Array{
        Float64,
        1,
    },
    ins::Array{
        Int64,
        1,
    },
)::Float64
    
    sum_ = 0.0
    
    @inbounds @fastmath @simd for index in 1:length(values)
        
        if ins[index] == 1
        
            sum_ += values[index]
            
        end
        
    end
    
    sum_
    
end
