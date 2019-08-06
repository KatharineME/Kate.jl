using CSV

include("compute_gene_set_enrichment.jl")
include("read_gmt.jl")

function gsea(
    gene_x_sample_tsv_file_path::String,
    gmt_file_paths::Array{
        String,
        1,
    },
    output_directory_path::String,
)

    gene_x_sample = CSV.read(gene_x_sample_tsv_file_path)
    
    gene_set_genes = read_gmt(gmt_file_paths)
    
    gene_set_x_sample = compute_gene_set_enrichment(
        gene_x_sample,
        gene_set_genes,
    )
    
    mkpath(output_directory_path)
    
    CSV.write(
        joinpath(
            output_directory_path,
            "gene_set_x_sample.tsv",
        ),
        gene_set_x_sample;
        delim='\t',
    )
    
    nothing

end
