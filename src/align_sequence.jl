include("print_and_run_cmd.jl")


function align_sequence(
    _1_fq_gz::String,
    _2_fq_gz::String,
    sample_name::String,
    dna_fa_gz_mmi::String,
    bam::String,
    n_job::Int,
)

    println("Aligning sequence ...")

    bam_prefix::String, bam_extension::String = splitext(bam)

    output_dir::String = splitdir(bam)[1]

    mkpath(output_dir)

    no_markdup_bam::String = "$bam_prefix.no_markdup$bam_extension"

    print_and_run_cmd(pipeline(
        `minimap2 -x sr -t $n_job -R "@RG\tID:$sample_name\tSM:$sample_name" -a $dna_fa_gz_mmi $_1_fq_gz $_2_fq_gz`,
        `samtools sort -n --threads $n_job`,
        `samtools fixmate -m --threads $n_job - -`,
        `samtools sort --threads $n_job`,
        no_markdup_bam,
    ))

    print_and_run_cmd(`samtools markdup --threads $n_job -s $no_markdup_bam $bam`)

    print_and_run_cmd(`rm --force $no_markdup_bam`)

    print_and_run_cmd(`samtools index -@ $n_job $bam`)

    print_and_run_cmd(pipeline(
        `samtools flagstat --threads $n_job $bam`,
        "$bam.flagstat",
    ))

end
