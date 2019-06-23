function trim_sequence(
        _1_fq_gz::String,
        _2_fq_gz::String,
        n_job::Int,
        output_dir::String,
    )

    println("Trimming sequence ...")

    Kraft.print_and_run_cmd(`skewer --threads $n_job -x AGATCGGAAGAGC $_1_fq_gz $_2_fq_gz --compress --output $output_dir`)

end
