function check_sequence(
        fq_gzs::Vector{String},
        n_job::Int,
        output_dir::String,
    )

    println("Checking sequence ...")

    kraft.print_and_run_cmd(`fastqc --threads $n_job $fq_gzs --outdir $output_dir`)

end
