include("print_and_run_cmd.jl")


function count_transcript(
    _1_fq_gz::String,
    _2_fq_gz::String,
    cdna_fa_gz::String,
    output_dir::String,
    n_job::Int,
)

    println("Counting transcript ...")

    cdna_fa_gz_kallisto_index::String = "$cdna_fa_gz.kallisto_index"

    if !ispath(cdna_fa_gz_kallisto_index)

        print_and_run_cmd(`kallisto index --index $cdna_fa_gz_kallisto_index $cdna_fa_gz`)

    end

    mkpath(output_dir)

    print_and_run_cmd(`kallisto quant --threads $n_job --index $cdna_fa_gz_kallisto_index --output-dir $output_dir $_1_fq_gz $_2_fq_gz`)

    nothing

end
