using Dates

include("print_and_run_cmd.jl")


function align_sequence(
    _1_fastq_gz_file_path::String,
    _2_fastq_gz_file_path::String,
    sample_name::String,
    dna_fasta_gz_file_path::String,
    bam_file_path::String,
    n_job::Int,
    n_gb_memory_per_job::Int,
)

    start_time = now()

    println("($start_time) Aligning sequence...")

    dna_fasta_gz_mmi_file_path = "$dna_fasta_gz_file_path.mmi"

    if !isfile(dna_fasta_gz_mmi_file_path)

        print_and_run_cmd(`minimap2 -t $n_job -d $dna_fasta_gz_mmi_file_path $dna_fasta_gz_file_path`)

    end

    if isfile(bam_file_path)

        error("$bam_file_path exists.")

    end

    mkpath(splitdir(bam_file_path)[1])

    print_and_run_cmd(pipeline(
        `minimap2 -x sr -t $n_job -K $(n_gb_memory_per_job)G -R "@RG\tID:$sample_name\tSM:$sample_name" -a $dna_fasta_gz_mmi_file_path $_1_fastq_gz_file_path $_2_fastq_gz_file_path`,
        `samtools sort --threads $n_job -m $(n_gb_memory_per_job)G -n`,
        `samtools fixmate --threads $n_job -m - -`,
        `samtools sort --threads $n_job -m $(n_gb_memory_per_job)G`,
        "$bam_file_path.tmp",
    ))

    print_and_run_cmd(`samtools markdup --threads $n_job -s $bam_file_path.tmp $bam_file_path`)

    rm("$bam_file_path.tmp")

    print_and_run_cmd(`samtools index -@ $n_job $bam_file_path`)

    print_and_run_cmd(pipeline(
        `samtools flagstat --threads $n_job $bam_file_path`,
        "$bam_file_path.flagstat",
    ))

    end_time = now()

    run_time = canonicalize(Dates.CompoundPeriod(end_time - start_time))

    println("($end_time) Done in $run_time.")

    return nothing

end
