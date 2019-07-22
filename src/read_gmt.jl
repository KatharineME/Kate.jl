function read_gmt(gmt_file_path::String)
    
    gene_set_dict = Dict{String, Array{String, 1}}()
    
    for line in readlines(gmt_file_path)
        
        line_split = String.(split(
            line,
            '\t',
        ))
        
        gene_set_dict[line_split[1]] = line_split[3:end]
        
    end
    
    gene_set_dict
    
end