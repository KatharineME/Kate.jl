using Dates

include("print_and_run_cmd.jl")


function align_sequence(
    _1_fastq_gz::String,
    _2_fastq_gz::String,
    sample_name::String,
    dna_fasta_gz::String,
    bam::String,
    n_job::Int,
    gb_memory::Int,
)

    start_time = now()

    println("($start_time) Aligning sequence ...")

    dna_fasta_gz_mmi = "$dna_fasta_gz.mmi"

    if !ispath(dna_fasta_gz_mmi)

        print_and_run_cmd(`minimap2 -t $n_job -d $dna_fasta_gz_mmi $dna_fasta_gz`)

    end

    if isfile(bam)

        error("$bam exists.")
    
    end

    mkpath(splitdir(bam)[1])

    print_and_run_cmd(pipeline(
        `minimap2 -x sr -t $n_job -K $(gb_memory)G -R "@RG\tID:$sample_name\tSM:$sample_name" -a $dna_fasta_gz_mmi $_1_fastq_gz $_2_fastq_gz`,
        `samtools sort --threads $n_job -m $(gb_memory)G -n`,
        `samtools fixmate --threads $n_job -m - -`,
        `samtools sort --threads $n_job -m $(gb_memory)G`,
        "$bam.tmp",
    ))

    print_and_run_cmd(`samtools markdup --threads $n_job -s $bam.tmp $bam`)

    rm("$bam.tmp")

    print_and_run_cmd(`samtools index -@ $n_job $bam`)

    print_and_run_cmd(pipeline(
        `samtools flagstat --threads $n_job $bam`,
        "$bam.flagstat",
    ))

    end_time = now()

    run_time = canonicalize(Dates.CompoundPeriod(end_time - start_time))

    println("($end_time) Done in $run_time.")

end
