include("make_hits.jl")
include("sum_hit_scores.jl")
include("sum_hits.jl")


function compute_gene_set_enrichment(
    genes::Array{String, 1},
    scores::Array{Float64, 1},
    gene_set_genes::Array{String, 1};
    gene_index::Union{Dict{String, Int64}, Nothing}=nothing,
    sort_scores::Bool=true,
)
    
    n_gene = length(genes)
    
    cumulative_sum = Array{Float64, 1}(
        undef,
        n_gene,
    )
    
    if gene_index === nothing
        
        hits = make_hits(
            genes,
            gene_set_genes,
        )
    
    else
        
        hits = make_hits(
            gene_index,
            gene_set_genes,
        )
        
    end

    if sort_scores

        sort!(
            scores,
            rev=true,
        )

    end

    abs_scores = abs.(scores)
        
    hit_scores_sum = sum_hit_scores(
        abs_scores,
        hits,
    )
    
    d_down = -1 / (n_gene - sum_hits(hits))
    
    value = 0.0
    
    auc = 0.0
    
    min_ = 0.0
    
    max_ = 0.0
    
    @inbounds @fastmath @simd for index in 1:n_gene
        
        if hits[index] == 1
            
            value += abs_scores[index] / hit_scores_sum
            
        else
            
            value += d_down
            
        end
        
        cumulative_sum[index] = value
        
        auc += value
        
        if value < min_
            
            min_ = value
            
        elseif max_ < value
            
            max_ = value
            
        end
            
    end
    
    cumulative_sum, auc, min_, max_
    
end


function compute_gene_set_enrichment(
    genes::Array{String, 1},
    scores::Array{Float64, 1},
    gene_set_dict::Dict{String, Array{String, 1}};
    sort_scores::Bool=true,

)
    
    if length(gene_set_dict) < 5
        
        gene_index = nothing
        
    else
        
        gene_index = Dict(gene=>index for (gene, index) in zip(
            genes,
            1:length(genes),
        ))
        
    end

    if sort_scores

        sort!(
            scores,
            rev=true,
        )

        sort_scores = false

    end

    gene_set_result_dict = Dict{String, Dict{String, Union{Float64, Array{Float64, 1}}}}()

    for (gene_set_name, gene_set_genes) in gene_set_dict
        
        cumulative_sum, auc, min_, max_ = compute_gene_set_enrichment(
            genes,
            scores,
            gene_set_genes;
            gene_index=gene_index,
            sort_scores=false,
        )
        
        gene_set_result_dict[gene_set_name] = Dict(
            "cumulative_sum"=>cumulative_sum,
            "auc"=>auc,
            "min"=>min_,
            "max"=>max_,
        )
        
    end

    gene_set_result_dict
    
end
