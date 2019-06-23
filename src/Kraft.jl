module Kraft

println(@__FILE__)

for name in readdir(@__DIR__)

    if name != splitdir(@__FILE__)[end] && endswith(name, ".jl")

        include(name)

    end

end

end
