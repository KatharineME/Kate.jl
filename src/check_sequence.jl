include("print_and_run_cmd.jl")


function check_sequence(
    fq_gzs::Tuple{Vararg{String}},
    output_dir::String,
    n_job::Int,
)

    println("Checking sequence ...")

    mkpath(output_dir)

    print_and_run_cmd(`fastqc --threads $(minimum((length(fq_gzs), n_job))) --outdir $output_dir $fq_gzs`)

    nothing

end
