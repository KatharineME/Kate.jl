module analyze_sequence

for name::String in readdir(@__DIR__)
    if endswith(name, ".jl") && name != splitdir(@__FILE__)[end]

        include(name)

    end

end

end
