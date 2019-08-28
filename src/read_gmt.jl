function read_gmt(gmt_file_path::String)

    gene_set_genes = Dict{String,Vector{String},}()

    for line in readlines(gmt_file_path)

        line_split = String.(split(line, '\t',))

        gene_set_genes[line_split[1]] = line_split[3:end]

    end

    return gene_set_genes
    
end


function read_gmt(gmt_file_paths::Vector{String})

    gene_set_genes = Dict{String,Vector{String},}()

    for gmt_file_path in gmt_file_paths

        gene_set_genes = merge(gene_set_genes, read_gmt(gmt_file_path),)

    end

    return gene_set_genes

end
