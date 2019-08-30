using CSV

include("combine_gene_sets_dn_up.jl")
include("compute_set_enrichment.jl")
include("read_gmt.jl")

function gsea(
    gene_x_sample_tsv_file_path::String,
    gmt_file_paths::Vector{String},
    output_directory_path::String;
    sample_normalization_method::Union{Nothing,String} = "-0-",
    gene_sets::Union{Nothing,Vector{String}} = nothing,
    n_required_gene_set_element::Int64 = 1,
    weight::Union{Nothing,Real} = nothing,
    statistic::String = "ks",
)

    gene_x_sample = CSV.read(gene_x_sample_tsv_file_path)

    genes = gene_x_sample[!, 1]
    
    gene_set_genes = read_gmt(gmt_file_paths)

    if gene_sets !== nothing

        Dict(gene_set => gene_set_genes[gene_set] for gene_set in gene_sets if haskey(
            gene_set_genes,
            gene_set
        ))

    end

    gene_set_genes = Dict(gene_set => genes_ for (
        gene_set,
        genes_
    ) = gene_set_genes if n_required_gene_set_element < length(intersect(
        genes_,
        genes
    )))

    sample_normalization_method, nothing, "log", "rank", "-0-", "0-1", weight

    gene_set_x_sample = combine_gene_sets_dn_up(compute_set_enrichment(
        gene_x_sample,
        gene_set_genes,
        statistic,
    ))
    
    mkpath(output_directory_path)
    
    CSV.write(
        joinpath(output_directory_path, "gene_set_x_sample.tsv",),
        gene_set_x_sample;
        delim = '\t',
    )
    
    return nothing

end
