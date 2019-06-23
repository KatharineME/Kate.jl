function count_transcript(
        _1_fq_gz::String,
        _2_fq_gz::String,
        kallisto_index::String,
        n_job::Int,
        output_dir::String,
    )

  println("Counting transcript ...")

  Kraft.print_and_run_cmd(`kallisto quant --index $kallisto_index --output-dir $output_dir --threads $n_job $_1_fq_gz $_2_fq_gz`)

end
