function sum_vector(vector::Vector{Float64}, vector_01::Vector{Int64},)

    sum_ = 0.0

    @inbounds @fastmath @simd for index = 1:length(vector)

        if _01[index] == 1

            sum_ += vector[index]

        end

    end

    return sum_

end
