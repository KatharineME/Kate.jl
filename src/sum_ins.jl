function sum_ins(ins::Array{
    Int64,
    1
})
    
    sum_ = 0
    
    @inbounds @fastmath @simd for index in 1:length(ins)
        
        sum_ += ins[index]
        
    end
    
    sum_
    
end
