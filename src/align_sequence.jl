function align_sequence(
        _1_fq_gz::String,
        _2_fq_gz::String,
        sample_name::String,
        dna_fa_gz_mmi::String,
        bam::String,
        n_job::Int,
        output_dir::String,
    )

    println("Aligning sequence ...")

    kraft.print_and_run_cmd(pipeliene(
            `minimap2 -x sr -t $n_job -R "@RG\tID:$sample_name\tSM:$sample_name" -a $dna_fa_gz_mmi $_1_fq_gz $_2_fq_gz`,
            `samtools sort -n --threads $n_job`,
            `samtools fixmate -m --threads $n_job - -`,
            `samtools sort --threads $n_job`,
            "/tmp/no_markdup.bam"
            ))

    kraft.print_and_run_cmd(`samtools markdup -s --threads $n_job /tmp/no_markdup.bam $bam`)

    kraft.print_and_run_cmd(`rm --force /tmp/no_markdup.bam`)

    kraft.print_and_run_cmd(`samtools index -@ $n_job $bam`)

    kraft.print_and_run_cmd(pipeliene(`samtools flagstat --threads $n_job $bam`, "$bam.flagstat"))

end
