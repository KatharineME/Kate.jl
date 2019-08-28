function sort_elements(element_values::Vector, elements::Vector{String},)

    sort_indices = sortperm(element_values)

    return element_values[sort_indices], elements[sort_indices]

end
