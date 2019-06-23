include("print_and_run_cmd.jl")


function count_transcript(
    _1_fq_gz::String,
    _2_fq_gz::String,
    kallisto_index::String,
    output_dir::String,
    n_job::Int,
)

    println("Counting transcript ...")

    print_and_run_cmd(`kallisto quant --threads $n_job --index $kallisto_index --output-dir $output_dir $_1_fq_gz $_2_fq_gz`)

end
