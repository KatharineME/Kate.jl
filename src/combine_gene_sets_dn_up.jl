using DataFrames

function combine_gene_sets_dn_up(gene_set_x_element)

    gene_sets = gene_set_x_element[:, Symbol("Gene Set")]

    gene_set_x_element_ = DataFrame(Dict(n => Array{t, 1}() for (n, t) in zip(
        names(gene_set_x_element),
        eltypes(gene_set_x_element),
    )))

    for gene_set in gene_sets
        
        if endswith(
            gene_set,
            "_DN",
        ) || endswith(
            gene_set,
            "_UP",
        )
            
            gene_set_ = gene_set[1:end - 3]
            
            if gene_set_ in gene_set_x_element_[:, Symbol("Gene Set")]
                
                continue
                
            end
            
            dn_values = gene_set_x_element[gene_sets .== "$(gene_set_)_DN", 2:end]
            
            up_values = gene_set_x_element[gene_sets .== "$(gene_set_)_UP", 2:end]
            
            if size(
                dn_values,
                1,
            ) == 0 || size(
                up_values,
                1,
            ) == 0
                
                continue
                
            end
            
            combined_values = up_values .- dn_values
            
            insertcols!(
                combined_values,
                1,
                Symbol("Gene Set") => gene_set_,
            )
            
            push!(
                gene_set_x_element_,
                combined_values[1, :],
            )
            
        end
            
    end

    sort(
        vcat(
            gene_set_x_element,
            gene_set_x_element_,
        ),
        Symbol("Gene Set"),
    )

end
