include("print_and_run_cmd.jl")


function trim_sequence(
    _1_fq_gz::String,
    _2_fq_gz::String,
    output_prefix::String,
    n_job::Int,
)

    println("Trimming sequence ...")

    output_dir::String = splitdir(output_prefix)[1]

    mkpath(output_dir)

    print_and_run_cmd(`skewer --threads $n_job -x AGATCGGAAGAGC --compress --output $output_prefix $_1_fq_gz $_2_fq_gz`)

end
