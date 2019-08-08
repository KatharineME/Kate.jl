using DataFrames
using Distributed

include("make_ins.jl")
include("sum_values.jl")


function compute_gene_set_enrichment(
    gene_values::Array{
        Float64,
        1,
    },
    genes::AbstractArray{
        String,
        1,
    },
    gene_set_genes::Array{
        String,
        1,
    };
    sort_gene_values::Bool = true,
    gene_index::Union{
        Dict{
            String,
            Int64,
        },
        Nothing,
    } = nothing,
    compute_cumulative_sum::Bool = true,
)
    
    if sort_gene_values

        gene_values_sort_indices = sortperm(gene_values)

        gene_values = gene_values[gene_values_sort_indices]

        genes = genes[gene_values_sort_indices]

    end

    if gene_index === nothing
        
        ins = make_ins(
            genes,
            gene_set_genes,
        )
    
    else
        
        ins = make_ins(
            gene_index,
            gene_set_genes,
        )
        
    end

    gene_values_abs = abs.(gene_values)

    gene_values_abs_ins_sum = sum_values(
        gene_values_abs,
        ins,
    )

    n_gene = length(genes)
    
    d_down = -1 / (n_gene - sum_values(ins))
    
    value = 0.0

    if compute_cumulative_sum

        cumulative_sum = Array{
            Float64,
            1,
        }(
            undef,
            n_gene,
        )

    else

        cumulative_sum = nothing

    end
    
    ks = 0.0

    ks_abs = 0.0

    auc = 0.0
    
    @inbounds @fastmath @simd for index in 1:n_gene
        
        if ins[index] == 1
            
            value += gene_values_abs[index] / gene_values_abs_ins_sum
            
        else
            
            value += d_down
            
        end
        
        if compute_cumulative_sum

            cumulative_sum[index] = value

        end

        value_abs = abs(value)
        
        if ks_abs < value_abs

            ks = value

            ks_abs = value_abs

        end

        auc += value
            
    end
    
    cumulative_sum, ks, auc
    
end


function compute_gene_set_enrichment(
    gene_values::Array{
        Float64,
        1,
    },
    genes::AbstractArray{
        String,
        1,
    },
    gene_set_genes::Dict{
        String,
        Array{
            String,
            1,
        },
    };
    sort_gene_values::Bool = true,
)

    if sort_gene_values

        gene_values_sort_indices = sortperm(gene_values)

        gene_values = gene_values[gene_values_sort_indices]

        genes = genes[gene_values_sort_indices]

    end

    if length(gene_set_genes) < 10
        
        gene_index = nothing
        
    else
        
        gene_index = Dict(gene => index for (
            gene,
            index,
        ) in zip(
            genes,
            1:length(genes),
        ))
        
    end

    gene_set_enrichment = Dict{
        String,
        Tuple{
            Union{
                Array{
                    Float64,
                    1,
                },
                Nothing,
            },
            Float64,
            Float64,
        },
    }()

    for (
        gene_set,
        gene_set_genes_,
    ) in gene_set_genes

        gene_set_enrichment[gene_set] = compute_gene_set_enrichment(
            gene_values,
            genes,
            gene_set_genes_;
            sort_gene_values = false,
            gene_index = gene_index,
            compute_cumulative_sum = false,
        )

    end

    gene_set_enrichment
    
end


function compute_gene_set_enrichment(
    gene_x_sample::DataFrame,
    gene_set_genes::Dict{
        String,
        Array{
            String,
            1,
        },
    },
)
    
    genes = gene_x_sample[:, Symbol("Gene")]

    gene_set_x_sample = DataFrame(Symbol("Gene Set") => collect(keys(gene_set_genes)))

    for sample in names(gene_x_sample)[2:end]

        has_gene_value = findall(
            !ismissing,
            gene_x_sample[:, sample],
        )

        gene_set_enrichment = compute_gene_set_enrichment(
            Float64.(gene_x_sample[has_gene_value, sample]),
            genes[has_gene_value],
            gene_set_genes,
        )
        
        gene_set_x_sample[!, sample] = collect(enrichment[2] for enrichment in values(gene_set_enrichment))

    end

    gene_set_x_sample
    
end
