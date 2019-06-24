include("print_and_run_cmd.jl")


function align_sequence(
    _1_fq_gz::String,
    _2_fq_gz::String,
    sample_name::String,
    dna_fa_gz::String,
    bam::String,
    n_job::Int,
    job_gb_memory::Int,
)

    println("Aligning sequence ...")

    dna_fa_gz_mmi::String = "$dna_fa_gz.mmi"

    if !ispath(dna_fa_gz_mmi)

        print_and_run_cmd(`minimap2 -t $n_job -d $dna_fa_gz_mmi $dna_fa_gz`)

    end

    output_dir::String = splitdir(bam)[1]

    mkpath(output_dir)

    print_and_run_cmd(pipeline(
        `minimap2 -x sr -t $n_job -K $(job_gb_memory)G -R "@RG\tID:$sample_name\tSM:$sample_name" -a $dna_fa_gz_mmi $_1_fq_gz $_2_fq_gz`,
        `samtools sort --threads $n_job -m $(job_gb_memory)G -n`,
        `samtools fixmate --threads $n_job -m - -`,
        `samtools sort --threads $n_job -m $(job_gb_memory)G`,
        "$bam.tmp",
    ))

    print_and_run_cmd(`samtools markdup --threads $n_job -s $bam.tmp $bam`)

    rm("$bam.tmp")

    print_and_run_cmd(`samtools index -@ $n_job $bam`)

    print_and_run_cmd(pipeline(
        `samtools flagstat --threads $n_job $bam`,
        "$bam.flagstat",
    ))

    nothing

end
