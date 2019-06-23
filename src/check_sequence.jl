include("print_and_run_cmd.jl")


function check_sequence(
    fq_gzs::Tuple{String, String},
    output_dir::String,
    n_job::Int,
)

    println("Checking sequence ...")

    mkpath(output_dir)

    n_job = minimum((length(fq_gzs), n_job))

    print_and_run_cmd(`fastqc --threads $n_job --outdir $output_dir $fq_gzs`)

end
