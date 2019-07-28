include("make_hits.jl")
include("sum_hit_scores.jl")
include("sum_hits.jl")


function compute_gene_set_enrichment(
    scores::Array{Float64, 1},
    genes::AbstractArray{String, 1},
    gene_set_genes::Array{String, 1};
    sort_scores::Bool=true,
    gene_index::Union{Dict{String, Int64}, Nothing}=nothing,
    compute_cumulative_sum::Bool=true,
)
    
    if sort_scores

        sort_indices = sortperm(
            scores,
            rev=true,
        )

        scores = scores[sort_indices]

        genes = genes[sort_indices]

    end

    abs_scores = abs.(scores)

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

    hit_scores_sum = sum_hit_scores(
        abs_scores,
        hits,
    )

    n_gene = length(genes)
    
    d_down = -1 / (n_gene - sum_hits(hits))
    
    value = 0.0

    if compute_cumulative_sum

        cumulative_sum = Array{Float64, 1}(
            undef,
            n_gene,
        )

    else

        cumulative_sum = nothing

    end
    
    min_ = 0.0
    
    max_ = 0.0

    auc = 0.0
    
    @inbounds @fastmath @simd for index in 1:n_gene
        
        if hits[index] == 1
            
            value += abs_scores[index] / hit_scores_sum
            
        else
            
            value += d_down
            
        end
        
        if compute_cumulative_sum

            cumulative_sum[index] = value

        end
        
        if value < min_
            
            min_ = value
            
        elseif max_ < value
            
            max_ = value
            
        end

        auc += value
            
    end
    
    cumulative_sum, min_, max_, auc
    
end


function compute_gene_set_enrichment(
    scores::Array{Float64, 1},
    genes::AbstractArray{String, 1},
    gene_set_dict::Dict{String, Array{String, 1}};
    sort_scores::Bool=true,
)

    if sort_scores

        sort_indices = sortperm(
            scores,
            rev=true,
        )

        scores = scores[sort_indices]

        genes = genes[sort_indices]

    end

    if length(gene_set_dict) < 5
        
        gene_index = nothing
        
    else
        
        gene_index = Dict(gene=>index for (gene, index) in zip(
            genes,
            1:length(genes),
        ))
        
    end

    gene_set_enrichment_dict = Dict{String, Tuple{Union{Array{Float64, 1}, Nothing}, Float64, Float64, Float64}}()

    for (gene_set_name, gene_set_genes) in gene_set_dict

        gene_set_enrichment_dict[gene_set_name] = compute_gene_set_enrichment(
            scores,
            genes,
            gene_set_genes;
            sort_scores=false,
            gene_index=gene_index,
            compute_cumulative_sum=false,
        )        

    end

    gene_set_enrichment_dict
    
end
