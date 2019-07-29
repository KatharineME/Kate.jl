function read_gmt(gmt_file_path::String)
    
    gene_set_genes = Dict{String, Array{String, 1}}()
    
    for line in readlines(gmt_file_path)
        
        line_split = String.(split(
            line,
            '\t',
        ))
        
        gene_set_genes[line_split[1]] = line_split[3:end]
        
    end
    
    gene_set_genes
    
end


function read_gmt(gmt_file_paths::Array{String, 1})
    
    gene_set_genes = Dict{String, Array{String, 1}}()
    
    for gmt_file_path in gmt_file_paths
        
        gene_set_genes = merge(
            gene_set_genes,
            read_gmt(gmt_file_path),
        )
        
    end
    
    gene_set_genes
    
end
