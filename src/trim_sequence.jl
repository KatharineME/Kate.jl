include("print_and_run_cmd.jl")


function trim_sequence(
    _1_fq_gz::String,
    _2_fq_gz::String,
    output_dir::String,
    n_job::Int,
)

    println("Trimming sequence ...")

    print_and_run_cmd(`skewer --threads $n_job -x AGATCGGAAGAGC --compress --output $output_dir $_1_fq_gz $_2_fq_gz`)

end
