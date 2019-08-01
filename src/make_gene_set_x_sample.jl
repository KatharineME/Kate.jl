using CSV
using DataFrames

include("compute_gene_set_enrichment.jl")


function make_gene_set_x_sample(
    gene_x_sample::DataFrame,
    gene_set_genes::Dict{
        String,
        Array{
            String,
            1,
        },
    };
    gene_set_x_sample_file_path::Union{
        String,
        Nothing,
    } = nothing,
)
    
    genes = gene_x_sample[:, Symbol("Gene")]

    gene_set_x_sample = DataFrame(Symbol("Gene Set") => collect(keys(gene_set_genes)))

    for sample in names(gene_x_sample)[2:end]

        println(sample)

        has_gene_value = findall(
            !ismissing,
            gene_x_sample[:, sample],
        )

        gene_set_enrichment = compute_gene_set_enrichment(
            Float64.(gene_x_sample[has_gene_value, sample]),
            genes[has_gene_value],
            gene_set_genes,
        )
        
        gene_set_x_sample[!, sample] = collect(enrichment[4] for enrichment in values(gene_set_enrichment))

    end

    if gene_set_x_sample_file_path !== nothing
        
        CSV.write(
            gene_set_x_sample_file_path,
            gene_set_x_sample;
            delim = '\t',
        )
        
    end
    
    gene_set_x_sample
    
end
