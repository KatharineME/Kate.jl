include("print_and_run_cmd.jl")


function check_sequence(
    fq_gzs::Vector{String},
    output_dir::String,
    n_job::Int,
)

    println("Checking sequence ...")

    print_and_run_cmd(`fastqc --threads $n_job --outdir $output_dir $fq_gzs`)

end
