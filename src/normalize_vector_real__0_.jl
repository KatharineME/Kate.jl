using Statistics


function normalize_vector_real__0_(vector::Vector{T} where T <: Real)
    
    vector = Vector{Float64}(vector)
    
    is_not_nan = .!isnan.(vector)
    
    vector_not_nan = vector[is_not_nan]
    
    vector[is_not_nan] .= (vector_not_nan .- mean(vector_not_nan)) /
                          std(vector_not_nan)
    
    vector
    
end
