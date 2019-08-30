function sum_vector_number(
    vector_number::Vector{<:Number},
    vector_01::Vector{Int64},
)

    sum_ = eltype(vector_number)(0)

    @inbounds @fastmath @simd for index = 1:length(vector_number)

        if vector_01[index] == 1

            sum_ += vector_number[index]

        end

    end

    return sum_

end
