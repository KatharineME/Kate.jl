using Dates

include("print_and_run_cmd.jl")


function align_sequence(
    fq1::String,
    fq2::String,
    sa::String,
    fa::String,
    bam::String,
    n_job::Int,
    job_gb_memory::Int,
)

    start_time = now()

    println("($start_time) Aligning sequence ...")

    dna_fasta_gz_mmi::String = "$fa.mmi"

    if !ispath(dna_fasta_gz_mmi)

        print_and_run_cmd(`minimap2 -t $n_job -d $dna_fasta_gz_mmi $fa`)

    end

    output_dir::String = splitdir(bam)[1]

    mkpath(output_dir)

    print_and_run_cmd(pipeline(
        `minimap2 -x sr -t $n_job -K $(job_gb_memory)G -R "@RG\tID:$sa\tSM:$sa" -a $dna_fasta_gz_mmi $fq1 $fq2`,
        # `samtools sort --threads $n_job -m $(job_gb_memory)G -n`,
        `samtools sort --threads $n_job -n`,
        `samtools fixmate --threads $n_job -m - -`,
        # `samtools sort --threads $n_job -m $(job_gb_memory)G`,
        `samtools sort --threads $n_job`,
        "$bam.tmp",
    ))

    print_and_run_cmd(`samtools markdup --threads $n_job -s $bam.tmp $bam`)

    rm("$bam.tmp")

    print_and_run_cmd(`samtools index -@ $n_job $bam`)

    print_and_run_cmd(pipeline(`samtools flagstat --threads $n_job $bam`, "$bam.flagstat"))

    end_time = now()

    println("($end_time) Done in $(canonicalize(Dates.CompoundPeriod(end_time - start_time))).")

end
